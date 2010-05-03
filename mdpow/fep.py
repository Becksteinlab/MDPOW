# mdpow: fep.py
# Copyright (c) 2010 Oliver Beckstein
"""
:mod:`mdpow.fep` -- calculate free energy of solvation
======================================================

Set up and run free energy perturbation (FEP) calculations to
calculate the free energy of hydration of a solute in a solvent
box. The protocol follows the works of D. Mobley (`Free Energy
Tutorial`_) and M. Shirts, and uses Gromacs 4.0.x.

Required Input:

 - topology
 - equilibrated structure of the solvated molecule


.. _Free Energy Tutorial:
   http://www.dillgroup.ucsf.edu/group/wiki/index.php?title=Free_Energy:_Tutorial


Example
-------

see :mod:`mdpow`


User reference
--------------

Simulation setup and analysis of all FEP simulations is encapsulated
by a :class:`mdpow.fep.Gsolv` object. For the hydration free energy
there is a special class :class:`~mdpow.fep.Ghyd` and for the
solvation free energy in octanol there is
:class:`~mdpow.fep.Goct`. See the description of
:class:`~mdpow.fep.Gsolv` for methods common to both.

.. autoclass:: Ghyd
.. autoclass:: Goct
.. autoclass:: Gsolv
   :members:

A user really only needs to access classes derived from
:class:`mdpow.fep.Gsolv` such as ; all other classes and functions are
auxiliary and only of interest to developers.


Developer notes
---------------

Additional objects that support :class:`mdpow.fep.Gsolv`.

.. autoclass:: FEPschedule
   :members:
.. autofunction:: molar_to_nm3

.. autodata:: N_AVOGADRO
.. autodata:: kBOLTZ
.. autodata:: fep_templates

.. |NA| replace:: *N*\ :sub:`A`
.. |kB| replace:: *k*\ :sub:`B`
.. _NA NIST value: http://physics.nist.gov/cgi-bin/cuu/Value?na
.. _kB NIST value: http://physics.nist.gov/cgi-bin/cuu/Value?k

TODO
~~~~

- run minimization, NVT-equil, NPT-equil prior to production (probably
  use preprocessed topology from ``grompp -pp`` for portability)
  
  See `Free Energy Tutorial`_.

"""
from __future__ import with_statement

import os
from subprocess import call
import warnings
from pkg_resources import resource_filename
import cPickle

import numpy

import gromacs, gromacs.setup, gromacs.utilities
from gromacs.utilities import asiterable, AttributeDict, in_dir

import logging
logger = logging.getLogger('mdpow.fep')

import config


#: Avogadro's constant |NA| in mol^-1 (`NA NIST value`_).
N_AVOGADRO = 6.02214179e23
#: Boltzmann's constant |kB| in kJ mol^-1 (`kB NIST value`_).
kBOLTZ = 1.3806504e-23 *1e-3 * N_AVOGADRO

def molar_to_nm3(c):
    """Convert a concentration in Molar to nm^-3."""
    return c * N_AVOGADRO * 1e-24

#: Template mdp files for different stages of the FEP protocol. (add equilibration, too?)
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


class Gsolv(object):
    """Simulations to calculate and analyze the solvation free energy.

    Typical work flow::

       G = Gsolv(simulation='drug.simulation')           # continue from :mod:`mdpow.equil`
       G.setup(qscript=['my_template.sge', 'local.sh'])  # my_template.sge is user supplied
       G.qsub()    # run SGE job arrays as generated from my_template.sge
       G.analyze()
       G.plot()

    See :mod:`gromacs.qsub` for notes on how to write templates for
    queuing system scripts (in particular `queuing system templates`_).

    .. _queuing system templates:
       http://sbcb.bioch.ox.ac.uk/oliver/software/GromacsWrapper/html/gromacs/building_blocks.html?highlight=qsub#queuing-system-templates
    """

    solvent_default = "water"

    schedules = {'coulomb':
                 FEPschedule(name='coulomb',
                             description="dis-charging vdw+q --> vdw",
                             label='Coul',
                             couple_lambda0='vdw-q', couple_lambda1='vdw',
                             sc_alpha=0,      # linear scaling for coulomb
                             lambdas=[0, 0.25, 0.5, 0.75, 1.0],  # default values
                             ),
                 'vdw':
                 FEPschedule(name='vdw',
                             description="decoupling vdw --> none",
                             label='VDW',
                             couple_lambda0='vdw', couple_lambda1='none',
                             sc_alpha=0.5, sc_power=1.0, sc_sigma=0.3, # recommended values
                             lambdas=[0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6,  # defaults
                                      0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1],
                             ),
                 }


    def __init__(self, molecule=None, top=None, struct=None, **kwargs):
        """Set up Gsolv from input files or a equilibrium simulation.
        
        :Arguments:
           *molecule*
               name of the molecule for which the hydration free 
               energy is to be computed (as in the gromacs topology)
               [REQUIRED]
           *top*
               topology [REQUIRED]
           *struct* 
               solvated and equilibrated input structure [REQUIRED]
           *ndx*
               index file
           *dirname*
               directory to work under ['FEP/*solvent*']
           *solvent*
               name of the solvent (only used for path names); will not 
               affect the simulations
           *lambda_coulomb*
               list of lambdas for discharging: q+vdw --> vdw
           *lambda_vdw*
               list of lambdas for decoupling: vdw --> none
           *runtime*
               simulation time per window in ps [5000]
           *temperature*
               temperature in Kelvin of the simulation [300.0]
           *qscript*
               template or list of templates for queuing system scripts
               (see :data:`gromacs.config.templates` for details) [local.sh]
           *includes*
               include directories
           *simulation*
               Instead of providing the required arguments, obtain the input
               files from a :class:`mdpow.equil.Simulation` instance.
           *filename*
               Instead of providing the required arguments, load from pickle
               file
           *kwargs*
               other undocumented arguments (see source for the moment)
        """
        required_args = ('molecule', 'top', 'struct')

        filename = kwargs.pop('filename', None)
        simulation = kwargs.pop('simulation', None)
        solvent = kwargs.pop('solvent', self.solvent_default)

        if (None in (molecule, top, struct) and simulation is None) and not filename is None:
            # load from pickle file
            self.load(filename)
            self.filename = filename
            kwargs = {}    # for super
        else:
            if not simulation is None:
                # load data from Simulation instance
                self.molecule = simulation.molecule
                self.top = simulation.files.processed_topology or simulation.files.topology
                self.struct = simulation.files.MD_NPT
                self.ndx = simulation.files.ndx
                if simulation.solvent_type == solvent:
                    self.solvent_type = simulation.solvent_type
                else:
                    errmsg = "Solvent mismatch: simulation was run for %s but Gsolv is set up for %s" % \
                        (simulation.solvent_type, solvent)
                    logger.error(errmsg)
                    raise ValueError(errmsg)
            else:
                self.molecule = molecule   # should check that this is in top (?)
                self.top = top
                self.struct = struct
                self.ndx = kwargs.pop('ndx', None)
                self.solvent_type = solvent

            for attr in required_args:
                if self.__getattribute__(attr) is None:
                    raise ValueError("A value is required for %(attr)r." % vars())

            # fix struct (issue with qscripts being independent from rest of code)
            if not os.path.exists(self.struct):
                # struct is not reliable as it depends on qscript so now we just try everything...
                struct = gromacs.utilities.find_first(self.struct, suffices=['pdb', 'gro'])
                if struct is None:
                    logger.error("Starting structure %r does not exist.", self.struct)
                    raise IOError(errno.ENOENT, "Starting structure not found", self.struct)
                else:
                    logger.info("Found starting structure %r (instead of %r).", struct, self.struct)
                    self.struct = struct
                        
            self.Temperature = kwargs.pop('temperature', 300.0)
            self.qscript = kwargs.pop('qscript', ['local.sh'])
            self.deffnm = kwargs.pop('deffnm', 'md')


            self.mdp = kwargs.pop('mdp', fep_templates['production_mdp'])
            self.lambdas = {
                'coulomb': kwargs.pop('lambda_coulomb', self.schedules['coulomb'].lambdas),
                'vdw':     kwargs.pop('lambda_vdw', self.schedules['vdw'].lambdas),
                }
            self.runtime = kwargs.pop('runtime', 5000.0)   # ps
            self.dirname = kwargs.pop('dirname', os.path.join('FEP', self.solvent_type))
            self.includes = list(asiterable(kwargs.pop('includes',[]))) + [config.includedir]
            self.component_dirs = {'coulomb': os.path.join(self.dirname, 'Coulomb'),
                                   'vdw': os.path.join(self.dirname, 'VDW')}

            self.filename = filename or \
                            os.path.join(self.dirname, self.__class__.__name__ + '.fep')

            # other variables
            #: Results from the analysis
            self.results = AttributeDict(xvg=AttributeDict(),
                                         dvdl=AttributeDict(),
                                         DeltaA=AttributeDict(),
                                         )
            #: Generated run scripts
            self.scripts = AttributeDict()

            # sanity checks
            if os.path.exists(self.dirname):
                wmsg = "Directory %(dirname)r already exists --- will overwrite " \
                       "existing files." % vars(self)
                warnings.warn(wmsg)
                logger.warn(wmsg)

        super(Gsolv, self).__init__(**kwargs)

        logger.info("Solvation free energy calculation for molecule "
                    "%(molecule)s in solvent %(solvent_type)s." % vars(self))
        logger.info("Using directories under %(dirname)r: %(component_dirs)r" % vars(self))


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
           *qscript*
               (List of) template(s) for batch submission scripts; if not set then 
               the templates are used that were supplied to the constructor.
           *kwargs*
               Most kwargs are passed on to :func:`gromacs.setup.MD` although some
               are set to values that are required for the FEP functionality.
               Note: *runtime* is set from the constructor.
        """
        kwargs['mdrun_opts'] = " ".join([kwargs.pop('mdrun_opts',''), '-dgdl'])  # crucial for FEP!!        
        kwargs['includes'] = asiterable(kwargs.pop('includes',[])) + self.includes        
        qsubargs = kwargs.copy()
        qsubargs['dirname'] = self.dirname
        # handle templates separately (necessary for array jobs)
        qscripts = qsubargs.pop('sge', None) or self.qscript
        qscripts.extend(qsubargs.pop('qscript',[]))  # also allow canonical 'templates'        
        # make sure that the individual batch scripts are also written
        kwargs.setdefault('qscript', qscripts)
        
        for component, lambdas in self.lambdas.items():
            for l in lambdas:
                self._setup(component, l, **kwargs)
            # array job
            directories = [self.wdir(component, l) for l in lambdas]
            qsubargs['jobname'] = self.arraylabel(component)
            qsubargs['prefix'] = self.label(component)+'_'
            self.scripts[component] = gromacs.qsub.generate_submit_array(qscripts, directories, **qsubargs)
            logger.info("[%s] Wrote array job scripts %r", component, self.scripts[component])

        self.save(self.filename)
        logger.info("Saved state information to %r; reload later with G = %r.", self.filename, self)
        logger.info("Finished setting up all individual simulations. Now run them...")

    def _setup(self, component, lmbda, **kwargs):
        """Prepare the input files for an individual Gromacs runs."""

        # note that all arguments pertinent to the submission scripts should be in kwargs
        # and have been set in setup() so that it is easy to generate array job scripts
        logger.info("Preparing %(component)s for lambda=%(lmbda)g" % vars())

        wdir = self.wdir(component, lmbda)
        kwargs.setdefault('couple-intramol', 'no')
        kwargs.update(self.schedules[component].mdp_dict)  # sets soft core & lambda0/1 state
        kwargs.update(dirname=wdir, struct=self.struct, top=self.top,
                      mdp=self.mdp,
                      ndx=self.ndx,
                      mainselection=None,
                      runtime=self.runtime,
                      ref_t=self.Temperature,    # TODO: maybe not working yet, check _setup()
                      gen_temp=self.Temperature, # needed until gromacs.setup() is smarter
                      qname=self.tasklabel(component,lmbda),
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


    def collect(self, autosave=True):
        """Collect dV/dl from output"""
        
        from gromacs.formats import XVG
        def dgdl_xvg(*args):
            return os.path.join(*args + (self.deffnm + '.xvg',))

        logger.info("[%(dirname)s] Finding dgdl xvg files." % vars(self))
        for component, lambdas in self.lambdas.items():
            xvg_files = [dgdl_xvg(self.wdir(component, l)) for l in lambdas]
            self.results.xvg[component] = (numpy.array(lambdas),
                                           [XVG(xvg) for xvg in xvg_files])            
        if autosave:
            self.save()


    def analyze(self, c0=1.0, force=False, autosave=True):
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

        Data are stored in :attr:`Gsolv.results`.

        :Keywords:
          *c0*
              standard state concentration in mol L^-1 (i.e. M) [1.0]
          *force*
              reload raw data even though it is already loaded
          *autosave*
              save to the pickle file when results have been computed
        """

        import scipy.integrate

        if force or numpy.any(numpy.array(
            [len(xvgs) for (lambdas,xvgs) in self.results.xvg.values()]) == 0):
            # get data if any of the xvg lists have 0 entries
            self.collect(autosave=False)
        else:
            logger.info("Analyzing stored data.")
        
        self.results.DeltaA.total = 0.0   # total free energy difference
        for component, (lambdas, xvgs) in self.results.xvg.items():
            logger.info("[%s] Computing averages <dV/dl> for %d lambda values.",
                        component, len(lambdas))
            # for TI just get the average dv/dl value (in array column 1; col 0 is the time)
            # (This can take a while if the XVG is now reading the array from disk first time)
            Y = numpy.array([x.array[1].mean() for x in xvgs])
            DY = numpy.array([x.array[1].std()  for x in xvgs])
            self.results.dvdl[component] = (lambdas, Y, DY)
            self.results.DeltaA[component] = scipy.integrate.simps(Y, x=lambdas) 
            self.results.DeltaA.total += self.results.DeltaA[component]

        # hydration free energy Delta A = -(Delta A_coul + Delta A_vdw)
        self.results.DeltaA.total *= -1

        # standard state
        self.results.DeltaA.total += self.DeltaA0(c0)

        if autosave:
            self.save()
        
        # TODO: error estimate (e.g. from boot strapping from the raw data)
        logger.info("Hydration free energy %g kJ/mol", self.results.DeltaA.total)
        return self.results.DeltaA.total    


    def plot(self, **kwargs):
        """Plot the TI data with error bars.

        Run :meth:`mdpow.fep.Gsolv.analyze` first.

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

    def save(self, filename=None):
        """Save instance to a pickle file.

        The default filename is the name of the file that was last loaded from
        or saved to.
        """
        if filename is None:
            filename = self.filename
        else:
            self.filename = filename
        with open(filename, 'wb') as f:
            cPickle.dump(self, f, protocol=cPickle.HIGHEST_PROTOCOL)
        logger.debug("Instance pickled to %(filename)r" % vars())
        
    def load(self, filename=None):
        """Re-instantiate class from pickled file."""
        if filename is None:
            filename = self.filename        
        with open(filename, 'rb') as f:
            instance = cPickle.load(f)
        self.__dict__.update(instance.__dict__)
        logger.debug("Instance loaded from %(filename)r" % vars())

    def __repr__(self):
        return "%s(filename=%r)" % (self.__class__.__name__, self.filename)

class Ghyd(Gsolv):
    """Sets up and analyses MD to obtain the hydration free energy of a solute."""
    solvent_default = "water"

class Goct(Gsolv):
    """Sets up and analyses MD to obtain the solvation free energy of a solute in octanol."""
    solvent_default = "octanol"

