# ghyd.py
# Copyright (c) 2010 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.
"""
:mod:`ghyd` -- calculate free energy of hydration
=================================================

Set up and run free energy perturbation calculations to calculate the
free energy of hydration of a solute in a water box. The protocol
follows the works of D. Mobley (`Free Energy Tutorial`_) and
M. Shirts, and uses Gromacs 4.x.

Required Input:
- topology
- equilibrated structure of the solvated molecule

TODO
----

- run minimization, NVT-equil, NPT-equil prior to production (probably
  use preprocessed topology from ``grompp -pp`` for portability)
  
  See `Free Energy Tutorial`_.

.. _Free Energy Tutorial:
   http://www.dillgroup.ucsf.edu/group/wiki/index.php?title=Free_Energy:_Tutorial
"""
from __future__ import with_statement

import os
from subprocess import call
from pkg_resources import resource_filename

import numpy

import gromacs, gromacs.setup
from gromacs.utilities import asiterable, AttributeDict, in_dir

import logging
logger = logging.getLogger('Ligands.ghyd')

#: `Avogadro's constant (NIST)`_ in mol^-1.
#: .. _Avogadro's constant (NIST): http://physics.nist.gov/cgi-bin/cuu/Value?na
N_AVOGADRO = 6.02214179e23
#: `Boltzmann's constant (NIST)`_ in kJ mol^-1.
#: .. _Boltzmann's constant (NIST): http://physics.nist.gov/cgi-bin/cuu/Value?k
kBOLTZ = 1.3806504e-23 *1e-3 * N_AVOGADRO

def molar_to_nm3(c):
    """Convert a concentration in Molar to nm^-3."""
    return c * N_AVOGADRO * 1e-24

# add equilibration, too?
fep_templates = {
    'production_mdp': resource_filename(__name__, 'templates/fep_opls.mdp'),
    }

class FEPschedule(AttributeDict):
    """Describe mdp parameter choices as key - value pairs."""
    mdp_keywords = ('sc_alpha', 'sc_power', 'sc_sigma',
                    'couple_lambda0', 'couple_lambda1',
                    )

    @property
    def mdp_dict(self):
        """Dict of key-values that can be set in a mdp file."""
        return dict(((k,v) for k,v in self.items() if k in self.mdp_keywords))


class Ghyd(object):
    """Simulations to calculate the hydration free energy."""

    # TODO: make the class saveable/pickle

    schedules = {'coulomb':
                 FEPschedule(name='coulomb',
                             description="dis-charging vdw+q --> vdw",
                             label='Coul',
                             couple_lambda0='vdw-q', couple_lambda1='vdw',
                             sc_alpha=0,      # linear scaling for coulomb
                             ),
                 'vdw':
                 FEPschedule(name='vdw',
                             description="decoupling vdw --> none",
                             label='VDW',
                             couple_lambda0='vdw', couple_lambda1='none',
                             sc_alpha=0.5, sc_power=1.0, sc_sigma=0.3, # recommended values
                             ),
                 }

    def __init__(self, molecule, top, struct, **kwargs):
        """Prepare all input files.
        
        :Arguments:
           *molecule*
               name of the molecule for which the hydration free 
               energy is to be computed (as in the gromacs topology)
           *top*
               topology
           *struct* 
               solvated and equilibrated input structure 
           *dirname*
               directory to work under [.]
           *lambda_coulomb*
               list of lambdas for discharging: q+vdw --> vdw
           *lambda_vdw*
               list of lambdas for decoupling: vdw --> none
           *runtime*
               simulation time per window in ps [100]
           *temperature*
               temperature in Kelvin of the simulation [300.0]
           *templates*
               template or list of templates for queuing system scripts
               (see :data:`gromacs.config.templates` for details) [local.sh]
           *kwargs*
               other undocumented arguments (see source for the moment)
        """
        self.molecule = molecule   # should check that this is in top (?)
        self.top = top
        self.struct = struct
        self.Temperature = kwargs.pop('temperature', 300.0)
        self.templates = kwargs.pop('templates', ['local.sh'])
        self.deffnm = kwargs.setdefault('deffnm', 'md')

        kwargs.setdefault('mdp', fep_templates['production_mdp'])

        kwargs.setdefault('lambda_coulomb', [0, 0.25, 0.5, 0.75, 1.0])
        kwargs.setdefault('lambda_vdw', [0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6,
                                         0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1])
        # should be longer, this is for testing
        kwargs.setdefault('runtime', 100) # ps
        kwargs.setdefault('dirname', os.path.curdir)

        self.mdp = kwargs['mdp']
        self.lambdas = {'coulomb': kwargs['lambda_coulomb'],
                        'vdw': kwargs['lambda_vdw']}
        self.runtime = kwargs['runtime']
        self.dirname = kwargs['dirname']
        self.component_dirs = {'coulomb': os.path.join(self.dirname, 'Coulomb'),
                               'vdw': os.path.join(self.dirname, 'VDW')}

        logger.info("Hydration free energy calculation for molecule %(molecule)r." % vars(self))
        logger.info("Using directories under %(dirname)r: %(component_dirs)r" % vars(self))

        # other variables
        #: Results from the analysis
        self.results = AttributeDict(xvg=AttributeDict(),
                                     dvdl=AttributeDict(),
                                     DeltaA=AttributeDict(),
                                     )
        #: Generated run scripts
        self.scripts = AttributeDict()

    def wdir(self, component, lmbda):
        """Return path to the work directory for *component* and *lmbda*."""
        return os.path.join(self.component_dirs[component], "%04d" % (1000 * lmbda))

    def label(self, component):
        """Simple label for component, e.g. for use in filenames."""
        return self.schedules[component].label

    def tasklabel(self, component, lmbda):
        """Batch submission script name for a single task job."""
        return self.molecule[:3]+'_'+self.schedules[component].label+"%04d" % (1000 * lmbda)

    def arraylabel(self, component):
        """Batch submission script name for a job array."""
        return self.molecule[:3]+'_'+self.schedules[component].label

    def setup(self, **kwargs):
        """Prepare the input files for all Gromacs runs.

        :Keywords:
           *sge*
               (List of) template(s) for batch submission scripts; if not set then 
               the templates are used that were supplied to the constructor.
           *kwargs*
               Most kwargs are passed on to :func:`gromacs.setup.MD` although some
               are set to values that are required for the FEP functionality.
        """
        kwargs['mdrun_opts'] = " ".join([kwargs.pop('mdrun_opts',''), '-dgdl'])  # crucial for FEP!!        
        qsubargs = kwargs.copy()
        qsubargs['dirname'] = self.dirname
        # handle templates separately (necessary for array jobs)
        templates = qsubargs.pop('sge', None) or self.templates
        # make sure that the individual batch scripts are also written
        kwargs.setdefault('sge', templates)
        
        for component, lambdas in self.lambdas.items():
            for l in lambdas:
                self._setup(component, l, **kwargs)
            # array job
            directories = [self.wdir(component, l) for l in lambdas]
            qsubargs['jobname'] = self.arraylabel(component)
            qsubargs['prefix'] = self.label(component)+'_'
            self.scripts[component] = gromacs.qsub.generate_submit_array(templates, directories, **qsubargs)
            logger.info("[%s] Wrote array job scripts %r", component, self.scripts[component])
        logger.info("Finished setting up all individual simulations. Now run them...")

    def _setup(self, component, lmbda, **kwargs):
        """Prepare the input files for an individual  Gromacs runs."""

        # note that all arguments pertinent to the submission scripts should be in kwargs
        # and have been set in setup() so that it is easy to generate array job scripts
        logger.info("Preparing %(component)s for lambda=%(lmbda)g" % vars())

        wdir = self.wdir(component, lmbda)
        kwargs.setdefault('couple-intramol', 'no')
        kwargs.update(self.schedules[component].mdp_dict)  # sets soft core & lambda0/1 state
        kwargs.update(dirname=wdir, struct=self.struct, top=self.top,
                      mdp=self.mdp, # ... add ndx ?
                      mainselection=None,
                      runtime=self.runtime,
                      ref_t=self.Temperature,    # TODO: maybe not working yet, check _setup()
                      gen_temp=self.Temperature, # needed until gromacs.setup() is smarter
                      sgename=self.tasklabel(component,lmbda),
                      free_energy='yes',
                      couple_moltype=self.molecule,
                      init_lambda=lmbda,
                      delta_lambda=0,
                      )
        gromacs.setup.MD(**kwargs)

    def DeltaA0(self, c0=1.0, N=1):
        """Standard state correction to the free energy.

        The total standard hydration free energy is calculated by taking the
        change to standard concentration *c0* from the simulation box volume into
        account:

          DeltaA_std = Delta A_0 - Delta A = -kT ln V0/V = -kT ln 1/V*c0*N_A

        where V is the volume of the simulation cell.
        (This assumes that there is only a single molecule *N* = 1 for which the free
        energy change was computed.)

        .. Note:: This correction is not exact; see Michael Shirts' thesis for
                  the correct one (but the difference should be small).        
        """

        V = gromacs.cbook.get_volume(self.struct)
        V0 = float(N)/molar_to_nm3(c0)
        self.results.DeltaA.standardstate = -kBOLTZ * \
                                            self.Temperature * numpy.log(V0/V)
        logger.debug("Standard state input: V=%g nm^3  V0=%g nm^3 T=%g K",
                     V, V0, self.Temperature)
        logger.info("Standard state correction for c0=%g M (V/V0 = %.2f): "
                    "dA_std = %+.2f kJ/mol", c0, V/V0, self.results.DeltaA.standardstate)
        return self.results.DeltaA.standardstate


    def analyze(self, c0=1.0):
        """Extract dV/dl from output and calculate dG by TI.

        Thermodynamic integration (TI) is performed on the individual component
        window calculation (typically the Coulomb and the VDW part). The
        dV/dlambda graphs are integrated with Simpson's rule (and averaging of
        results if the number of datapoints is odd; see
        :func:`scipy.integrate.simps` for details).

        The total standard hydration free energy is calculated by taking the
        change to standard concentration from the simulation box volume into
        account:

            DeltaA_std = A_0 - A_sim = -kT ln V0/Vsim = -kT ln 1/Vsim*c0*N_A

        (This assumes that there is only a single molecule for which the free
        energy change was computed.)

        Data are stored in :attr:`Ghyd.results`.

        :Keywords:
          *c_standard*
              standard state concentration in mol L^-1 (i.e. M) [1.0]
        """
        from gromacs.formats import XVG
        import scipy.integrate

        def dgdl_xvg(*args):
            return os.path.join(*args + (self.deffnm + '.xvg',))

        logger.info("Finding dgdl xvg files.")
        for component, lambdas in self.lambdas.items():
            xvg_files = [dgdl_xvg(self.wdir(component, l)) for l in lambdas]
            self.results.xvg[component] = (numpy.array(lambdas),
                                           [XVG(xvg) for xvg in xvg_files])            
        self.results.DeltaA.total = 0.0   # total free energy difference
        for component, (lambdas, xvgs) in self.results.xvg.items():
            logger.info("[%s] Computing averages <dV/dl> for %d lambda values.",
                        component, len(lambdas))
            # for TI just get the average dv/dl value (in array column 1; col 0 is the time)
            Y = numpy.array([x.array[1].mean() for x in xvgs])
            DY = numpy.array([x.array[1].std()  for x in xvgs])
            self.results.dvdl[component] = (lambdas, Y, DY)
            self.results.DeltaA[component] = scipy.integrate.simps(Y, x=lambdas) 
            self.results.DeltaA.total += self.results.DeltaA[component]

        # hydration free energy Delta A = -(Delta A_coul + Delta A_vdw)
        self.results.DeltaA.total *= -1

        # standard state
        self.results.DeltaA.total += self.DeltaA0(c0)
        
        # TODO: error estimate (e.g. from boot strapping from the raw data)
        logger.info("Hydration free energy %g kJ/mol", self.results.DeltaA.total)
        return self.results.DeltaA.total


    def plot(self, **kwargs):
        """Plot the TI data with error bars.

        Run :meth:`Ligands.ghyd.Ghyd.analyze` first.

        All *kwargs* are passed on to :func:`pylab.errorbar`.
        """
        from pylab import subplot, errorbar, xlabel, ylabel, legend, xlim, title

        kwargs.setdefault('color', 'black')
        kwargs.setdefault('capsize', 0)

        try:
            if self.results.DeltaA.total is None:
                raise KeyError
        except KeyError:
            logger.info("Data were not analyzed yet -- doing that now... patience.")
            self.analyze()

        dvdl = self.results.dvdl
        nplots = len(dvdl)
        for i, (component, (x,y,dy)) in enumerate(dvdl.items()):
            iplot = i+1
            subplot(1, nplots, iplot)
            label = r"$\Delta A_{\rm{%s}} = %.2f$ kJ/mol" \
                    % (component, self.results.DeltaA[component])
            errorbar(x, y, yerr=dy, label=label, **kwargs)
            xlabel(r'$\lambda$')
            legend(loc='best')
            xlim(-0.05, 1.05)
        subplot(1, nplots, 1)
        ylabel(r'$dV/d\lambda$ in kJ/mol')
        title(r"Free energy difference $\Delta A^{0}$")
        subplot(1, nplots, 2)        
        title(r"for %s: %.2f kJ/mol" %
              (self.molecule, self.results.DeltaA.total))
        

    def qsub(self, script=None):
        """Submit a batch script locally.

        If *script* == ``None`` then take the first script (works well if only
        one template was provided).
        """
        from gromacs.qsub import relpath
        
        def choose_script_from(scripts):
            if script is None:
                s = scripts[0]
            elif script in scripts:
                s = script
            else:
                errmsg = "No script matching %(script)r in %(scripts)r" % vars()
                logger.error(errmsg)
                raise ValueError(errmsg)
            return s

        with in_dir(self.dirname, create=False):
            for component, scripts in self.scripts.items():
                s = relpath(choose_script_from(scripts), self.dirname) # relative to dirname
                cmd = ['qsub', s]
                logger.debug("[%s] submitting locally: %s", " ".join(cmd), component)
                rc = call(cmd)
                if rc != 0:
                    errmsg = "submitting job %(s)r failed with rc=%(rc)d" % vars()
                    logger.error(errmsg)
                    raise OSError(errmsg)

        logger.info("[%r] Submitted jobs locally for %r", self.dirname, self.scripts.keys())
