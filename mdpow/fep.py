# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
# mdpow: fep.py
# Copyright (c) 2010 Oliver Beckstein

"""
:mod:`mdpow.fep` -- Calculate free energy of solvation
======================================================

Set up and run free energy perturbation (FEP) calculations to
calculate the free energy of hydration of a solute in a solvent
box. The protocol follows the works of D. Mobley (`Free Energy
Tutorial`_) and M. Shirts, and uses Gromacs 4.0.x.

Required Input:

 - topology
 - equilibrated structure of the solvated molecule

See the docs for :class:`Gsolv` for details on what is calculated.


.. _Free Energy Tutorial:
   http://www.dillgroup.ucsf.edu/group/wiki/index.php?title=Free_Energy:_Tutorial

Differences to published protocols
----------------------------------

Some notes on how the approach encoded here differs from what others
(notabley Mobley) did:

- We use Gromacs 4.x and use the new decoupling feature ::

    couple-intramol = no

  which indicates that "intra-molecular interactions remain". It should (as I
  understand it) allow us to only calculate the free energy changes in solution
  so that we do not have to do an extra calculation in vacuo.
  http://www.mail-archive.com/gmx-users@gromacs.org/msg18803.html

  Mobley does an extra discharging calculation in vacuo and calculates

  .. math::

      \Delta A = \Delta A_{\mathrm{coul}}(\mathrm{vac}) - (\Delta A_{\mathrm{coul}}(\mathrm{sol}) + \Delta A_{\mathrm{vdw}}(\mathrm{sol}))

  (but also annihilates the interactions on the solute, corresponding to
  ``couple-intramol = yes``) whereas we do

  .. math::

      \Delta A = - (\Delta A_{\mathrm{coul}}(\mathrm{sol}) + \Delta A_{\mathrm{vdw}}(\mathrm{sol}))

- simulations are run as NPT (but make sure to use Langevin dynamics
  for temperature control)


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

.. autoclass:: Gsolv
   :members:
.. autoclass:: Ghyd
.. autoclass:: Goct
.. autoclass:: Gcyclo
.. autofunction:: pOW
.. autofunction:: pCW


Developer notes
---------------

A user really only needs to access classes derived from
:class:`mdpow.fep.Gsolv`; all other classes and functions are auxiliary and
only of interest to developers.

Additional objects that support :class:`mdpow.fep.Gsolv`.

.. autoclass:: FEPschedule
   :members:
.. autofunction:: molar_to_nm3
.. autofunction:: kcal_to_kJ
.. autofunction:: kJ_to_kcal

.. data:: N_AVOGADRO

          Avogadro's constant |NA| in mol\ :sup:`-1` (`NA NIST value`_).

.. data:: kBOLTZ

          Boltzmann's constant |kB| in kJ mol\ :sup:`-1` (`kB NIST value`_).

.. |NA| replace:: *N*\ :sub:`A`
.. |kB| replace:: *k*\ :sub:`B`
.. _NA NIST value: http://physics.nist.gov/cgi-bin/cuu/Value?na
.. _kB NIST value: http://physics.nist.gov/cgi-bin/cuu/Value?k
.. |^-1| replace:: \ :sup:`-1`
.. |^-3| replace:: \ :sup:`-3`
.. |^3|  replace:: \ :sup:`3`

TODO
~~~~

- run minimization, NVT-equil, NPT-equil prior to production (probably
  use preprocessed topology from ``grompp -pp`` for portability)

  See `Free Energy Tutorial`_.

"""
from __future__ import absolute_import, with_statement

import os
import errno
import copy
from subprocess import call
import warnings
from pkg_resources import resource_filename

import numpy

import gromacs, gromacs.utilities
try:
    import gromacs.setup
except (ImportError, OSError):
    raise ImportError("Gromacs installation not found, source GMXRC?")
from gromacs.utilities import asiterable, AttributeDict, in_dir, openany
from numkit.observables import QuantityWithError

import logging
logger = logging.getLogger('mdpow.fep')

from . import config
from .restart import Journalled
from . import kBOLTZ, N_AVOGADRO

def molar_to_nm3(c):
    """Convert a concentration in Molar to nm|^-3|."""
    return c * N_AVOGADRO * 1e-24

def bar_to_kJmolnm3(p):
    """Convert pressure in bar to kJ mol|^-1| nm|^-3|.

    1 bar = 1e5 J m|^-3|
    """
    return p * N_AVOGADRO * 1e-25

def kcal_to_kJ(x):
    """Convert a energy in kcal to kJ."""
    return 4.184 * x

def kJ_to_kcal(x):
    """Convert a energy in kJ to kcal."""
    return x / 4.184


class FEPschedule(AttributeDict):
    """Describe mdp parameter choices as key - value pairs.

    The FEP schedule can be loaded from a configuration parser with
    the static method :meth:`FEPschedule.load`.

    See the example runinput file for an example. It contains the sections::

       [FEP_schedule_Coulomb]
       name = Coulomb
       description = dis-charging vdw+q --> vdw
       label = Coul
       couple_lambda0 = vdw-q
       couple_lambda1 = vdw
       # soft core alpha: linear scaling for coulomb
       sc_alpha = 0
       lambdas = 0, 0.25, 0.5, 0.75, 1.0


       [FEP_schedule_VDW]
       name = vdw
       description = decoupling vdw --> none
       label = VDW
       couple_lambda0 = vdw
       couple_lambda1 = none
       # recommended values for soft cores (Mobley, Shirts et al)
       sc_alpha = 0.5
       sc_power = 1.0
       sc_sigma = 0.3
       lambdas = 0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1

    """
    mdp_keywords = dict((('sc_alpha', float),
                         ('sc_power', float),
                         ('sc_sigma', float),
                         ('couple_lambda0', str),
                         ('couple_lambda1', str),
                         ))
    meta_keywords = dict((('name', str), ('description', str), ('label', str)))
    other_keywords = dict((('lambdas', list), ))

    @property
    def mdp_dict(self):
        """Dict of key-values that can be set in a mdp file."""
        return dict(((k,v) for k,v in self.items() if k in self.mdp_keywords))

    @staticmethod
    def load(cfg, section):
        """Initialize a :class:`FEPschedule` from the *section* in the configuration *cfg*"""

        keys = {}
        keys.update(FEPschedule.mdp_keywords)
        keys.update(FEPschedule.meta_keywords)
        keys.update(FEPschedule.other_keywords)

        cfg_get = {float: cfg.getfloat,
                   str: cfg.getstr,    # literal strings, no conversion of None (which we need for the MDP!)
                   list: cfg.getarray  # numpy float array from list
                   }
        def getter(type, section, key):
            try:
                return cfg_get[type](section, key)
            except ConfigParser.NoOptionError:
                return None
        return FEPschedule((key, getter(keytype, section, key)) for key,keytype in keys.items()
                           if getter(keytype, section, key) is not None)

    def __deepcopy__(self, memo):
        x = type(self)()
        for k,v in self.iteritems():
            x[k] = copy.deepcopy(v, memo)
        return x

class Gsolv(Journalled):
    """Simulations to calculate and analyze the solvation free energy.

    :math:`\Delta A` is computed from the decharging and the
    decoupling step. With our choice of ``lambda=0`` being the fully
    interacting and ``lambda=1`` the non-interacting state, it is computed
    as

    .. math::

            \Delta A = -(\Delta A_{\mathrm{coul}} + \Delta A_{\mathrm{vdw}})

    With this protocol, the concentration in the liquid and in the gas
    phase is the same. (Under the assumption of ideal solution/ideal
    gas behaviour this directly relates to the Ben-Naim 1M/1M standard
    state.)

    (We neglect the negligible correction :math:`-kT \ln
    V_x/V_{\mathrm{sim}} = -kT \ln(1 - v_s/V_{\mathrm{sim}})` where
    :math:`V_x` is the volume of the system without the solute but the
    same number of water molecules as in the fully interacting case
    [see Michael Shirts' Thesis, p82].)


    Typical work flow::

       G = Gsolv(simulation='drug.simulation')           # continue from :mod:`mdpow.equil`
       G.setup(qscript=['my_template.sge', 'local.sh'])  # my_template.sge is user supplied
       G.qsub()    # run SGE job arrays as generated from my_template.sge
       G.analyze()
       G.plot()

    See :mod:`gromacs.qsub` for notes on how to write templates for
    queuing system scripts (in particular `queuing system templates`_).

    .. _queuing system templates:
       http://gromacswrapper.readthedocs.io/en/latest/gromacs/blocks/qsub.html#queuing-system-templates
    """

    topdir_default = "FEP"
    dirname_default = os.path.curdir
    solvent_default = "water"
    #: Check list of all methods that can be run as an independent protocol; see also
    #: :meth:`Simulation.get_protocol` and :class:`restart.Journal`
    protocols = ["setup", "fep_run"]

    # TODO: initialize from default cfg
    schedules_default = {'coulomb':
                         FEPschedule(name='coulomb',
                                     description="dis-charging vdw+q --> vdw",
                                     label='Coul',
                                     couple_lambda0='vdw-q', couple_lambda1='vdw',
                                     sc_alpha=0,      # linear scaling for coulomb
                                     lambdas=[0.0, 0.25, 0.5, 0.75, 1.0],  # default values
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

    #: Default Gromacs *MDP* run parameter file for FEP.
    #: (The file is part of the package and is found with :func:`mdpow.config.get_template`.)
    mdp_default = 'fep_opls.mdp'


    def __init__(self, molecule=None, top=None, struct=None, method="BAR", **kwargs):
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
               name of the solvent (only used for path names); will not affect
               the simulations because the input coordinates and topologies are
               determined by either *simulation* or *molecule*, *top*, and
               *struct*. If *solvent* is not provided then it is set to
               :attr:`solvent_default`.
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
           *method*
               "TI" for thermodynamic integration or "BAR" for Bennett acceptance
               ratio; using "BAR" in Gromacs also writes TI data so this is the
               default.
           *mdp*
               MDP file name (if no entry then the package defaults are used)
           *filename*
               Instead of providing the required arguments, load from pickle
               file. If either *simulation* or *molecule*, *top*, and *struct*
               are provided then simply set the attribute :attr:`filename` of
               the default checkpoint file to *filename*.
           *basedir*
               Prepend *basedir* to all filenames; ``None`` disables [.]
           *permissive*
               Set to ``True`` if you want to read past corrupt data in output
               xvg files (see :class:`gromacs.formats.XVG` for details); note that
               *permissive*=``True`` can lead to **wrong results**. Overrides
               the value set in a loaded pickle file. [``False``]
           *stride*
               collect every *stride* data line, see :meth:`Gsolv.collect` [1]
           *kwargs*
               other undocumented arguments (see source for the moment)

        .. Note::

           Only a subset of the FEP parameters can currently be set with
           keywords (lambdas); other parameters such as the soft cores are
           currently fixed. For advanced use: pass a dict with keys "coulomb"
           and "vdw" and values of :class:`FEPschedule` in the keyword
           *schedules*. The *lambdas* will override these settings (which
           typically come from a cfg). This method allows setting all
           parameters.

        """
        required_args = ('molecule', 'top', 'struct')

        # should this be below somewhere?
        if not method in ("TI", "BAR"):
            raise ValueError("method can only be TI or BAR")
        self.method = method

        filename = kwargs.pop('filename', None)
        basedir = kwargs.pop('basedir', os.path.curdir)  # all other dirs relative to basedir
        simulation = kwargs.pop('simulation', None)
        solvent = kwargs.pop('solvent', self.solvent_default)
        if (None in (molecule, top, struct) and simulation is None) and filename is not None:
            # load from pickle file
            self.load(filename)
            self.filename = filename
            kwargs = {}    # for super
        else:
            if simulation is not None:
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

            self.mdp = kwargs.pop('mdp', config.get_template(self.mdp_default))

            # schedules (deepcopy because we might modify)
            self.schedules = copy.deepcopy(self.schedules_default)
            schedules = kwargs.pop('schedules', {})
            self.schedules.update(schedules)
            self.lambdas = {
                'coulomb': kwargs.pop('lambda_coulomb', self.schedules['coulomb'].lambdas),
                'vdw':     kwargs.pop('lambda_vdw', self.schedules['vdw'].lambdas),
                }
            self.runtime = kwargs.pop('runtime', 5000.0)   # ps
            self.dirname = kwargs.pop('dirname', self.dirname_default)
            self.includes = list(asiterable(kwargs.pop('includes',[]))) + [config.includedir]
            self.component_dirs = {'coulomb': os.path.join(self.dirname, 'Coulomb'),
                                   'vdw': os.path.join(self.dirname, 'VDW')}

            # for analysis
            self.stride = kwargs.pop('stride', 1)

            # other variables
            #: Results from the analysis
            self.results = AttributeDict(xvg=AttributeDict(),
                                         dvdl=AttributeDict(),
                                         DeltaA=AttributeDict(),  # contains QuantityWithError
                                         )
            #: Generated run scripts
            self.scripts = AttributeDict()

            # sanity checks
            if os.path.exists(self.dirname):
                wmsg = "Directory %(dirname)r already exists --- will overwrite " \
                       "existing files." % vars(self)
                warnings.warn(wmsg)
                logger.warn(wmsg)

        # overrides pickle file so that we can run from elsewhere
        if not basedir is None:
            self.basedir = os.path.realpath(basedir)
        else:
            self.basedir = None

        try:
            self.filename = os.path.abspath(self.filename)
        except (AttributeError, TypeError):
            # default filename if none was provided
            self.filename = self.frombase(self.dirname, self.__class__.__name__+os.extsep+'fep')

        # override pickle file for this dangerous option: must be set
        # on a case-by-case basis
        self.permissive = kwargs.pop('permissive', False)

        logger.info("Solvation free energy calculation for molecule "
                    "%(molecule)s in solvent %(solvent_type)s.", vars(self))
        logger.info("Base directory is %(basedir)r", vars(self))
        logger.info("Using setup directories under %(dirname)r: %(component_dirs)r", vars(self))
        logger.info("Default checkpoint file is %(filename)r", vars(self))
        logger.debug("Coulomb lambdas = %(coulomb)r", self.lambdas)
        logger.debug("VDW lambdas = %(vdw)r", self.lambdas)

        super(Gsolv, self).__init__(**kwargs)

    def frombase(self, *args):
        """Return path with :attr:`Gsolv.basedir` prefixed."""
        # wrap paths with frombase() and hopefully this allows us fairly
        # flexible use of the class, especially for analysis
        if self.basedir is None:
            return os.path.join(*args)
        return os.path.join(self.basedir, *args)

    def wname(self, component, lmbda):
        """Return name of the window directory itself.

        Typically something like ``VDW/0000``, ``VDW/0500``, ..., ``Coulomb/1000``
        """
        return os.path.join(self.component_dirs[component], "%04d" % (1000 * lmbda))

    def wdir(self, component, lmbda):
        """Return rooted path to the work directory for *component* and *lmbda*.

        (Constructed from :meth:`frombase` and :meth:`wname`.)
        """
        return self.frombase(self.wname(component, lmbda))

    def label(self, component):
        """Simple label for component, e.g. for use in filenames."""
        return self.schedules[component].label

    def tasklabel(self, component, lmbda):
        """Batch submission script name for a single task job."""
        return self.molecule[:3]+'_'+self.schedules[component].label+"%04d" % (1000 * lmbda)

    def arraylabel(self, component):
        """Batch submission script name for a job array."""
        return self.molecule[:3]+'_'+self.schedules[component].label

    def fep_dirs(self):
        """Generator for all simulation sub directories"""
        for component, lambdas in self.lambdas.items():
            for l in lambdas:
                yield self.wdir(component, l)

    def setup(self, **kwargs):
        """Prepare the input files for all Gromacs runs.

        :Keywords:
           *qscript*
               (List of) template(s) for batch submission scripts; if not set then
               the templates are used that were supplied to the constructor.
           *kwargs*
               Most kwargs are passed on to :func:`gromacs.setup.MD` although some
               are set to values that are required for the FEP functionality.

               *mdrun_opts*
                  list of options to :program:`mdrun`; ``-dhdl`` is always added
                  to this list as it is required for the thermodynamic integration
                  calculations
               *includes*
                  list of directories where Gromacs input files can be found

               The following keywords cannot be changed: *dirname*, *jobname*,
               *prefix*, *runtime*, *deffnm*. The last two can be specified when
               constructing :class:`Gsolv`.

        .. SeeAlso:: :func:`gromacs.setup.MD` and
                     :func:`gromacs.qsub.generate_submit_scripts`

        .. versionchanged:: 0.6.0
           Gromacs now uses option ``-dhdl`` instead of ``-dgdl``.
        """
        self.journal.start('setup')

        # -dgdl for FEP output (although that seems to have been changed to -dHdl in Gromacs 4.5.3)
        # NOW use -dhdl
        kwargs['mdrun_opts'] = " ".join([kwargs.pop('mdrun_opts',''), '-dhdl'])
        kwargs['includes'] = asiterable(kwargs.pop('includes',[])) + self.includes
        kwargs['deffnm'] = self.deffnm
        qsubargs = kwargs.copy()
        qsubargs['dirname'] = self.frombase(self.dirname)
        # handle templates separately (necessary for array jobs)
        qscripts = qsubargs.pop('sge', None) or self.qscript
        qscripts.extend(qsubargs.pop('qscript',[]))  # also allow canonical 'templates'
        # make sure that the individual batch scripts are also written
        kwargs.setdefault('qscript', qscripts)

        for component, lambdas in self.lambdas.items():
            if self.method == "TI":
                for l in lambdas:
                    # set up gromacs job for each FEP window in TI
                    params = self._setup(component, l, **kwargs)
            elif self.method == "BAR":
                xlambdas = [None] + list(lambdas) + [None]
                for l in lambdas:
                    foreign_lambdas = [xlambdas[xlambdas.index(l) - 1],
                                       xlambdas[xlambdas.index(l) + 1]]
                    # set up gromacs job for each FEP window in BAR
                    params = self._setup(component, l,
                                         foreign_lambdas=foreign_lambdas, **kwargs)
            else:
                raise ValueError("Unknown method {}".format(self.method))

            # generate queuing system script for array job
            directories = [self.wdir(component, l) for l in lambdas]
            qsubargs['jobname'] = self.arraylabel(component)
            qsubargs['prefix'] = self.label(component)+'_'
            self.scripts[component] = gromacs.qsub.generate_submit_array(qscripts, directories, **qsubargs)
            logger.info("[%s] Wrote array job scripts %r", component, self.scripts[component])

        self.journal.completed('setup')
        self.save(self.filename)
        logger.info("Saved state information to %r; reload later with G = %r.", self.filename, self)
        logger.info("Finished setting up all individual simulations. Now run them...")
        params.pop('struct', None)   # scrub window-specific params
        return params

    def _setup(self, component, lmbda, foreign_lambdas=None, **kwargs):
        """Prepare the input files for an individual Gromacs runs."""

        # note that all arguments pertinent to the submission scripts should be in kwargs
        # and have been set in setup() so that it is easy to generate array job scripts
        logger.info("Preparing %(component)s for lambda=%(lmbda)g" % vars())

        wdir = self.wdir(component, lmbda)
        kwargs.setdefault('couple-intramol', 'no')

        ### XXX Issue 20: if an entry is None then the dict will not be updated:
        ###     I *must* keep "none" as a legal string value
        kwargs.update(self.schedules[component].mdp_dict)  # sets soft core & lambda0/1 state

        # for BAR
        if foreign_lambdas is not None:
            # foreign_lambdas is only passed as an argument if BAR is
            # desired method, defaults to TI otherwise
            if foreign_lambdas[0] is None:
                feplambdas = foreign_lambdas[1]
            elif foreign_lambdas[1] is None:
                feplambdas = foreign_lambdas[0]
            else:
                feplambdas = "{0} {1}".format(foreign_lambdas[0], foreign_lambdas[1])

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
                          fep_lambdas=feplambdas,
                          )

        # for TI
        else:
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
        return gromacs.setup.MD(**kwargs)

    def dgdl_xvg(self, *args):
        """Return filename of the dgdl XVG file.

        Recognizes uncompressed, gzipped (gz), and bzip2ed (bz2)
        files.
g
        :Arguments:
           *args*
               joins the arguments into a path and adds the default
               filename for the dvdl file

        :Returns: Checks if a compressed file exists and returns
                  the appropriate filename.

        :Raises: :exc:`IOError` with error code ENOENT if no file
                 could be found

         """
        EXTENSIONS = ('', os.path.extsep+'bz2', os.path.extsep+'gz')
        root = os.path.join(*args + (self.deffnm + '.xvg',))
        for ext in EXTENSIONS:
            fn = root + ext
            if os.path.exists(fn):
                return fn
        logger.error("Missing dgdl.xvg file %(root)r.", vars())
        raise IOError(errno.ENOENT, "Missing dgdl.xvg file", root)

    def collect(self, stride=None, autosave=True, autocompress=True):
        """Collect dV/dl from output.

        :Keywords:
           *stride*
              read data every *stride* lines, ``None`` uses the class default
           *autosave*
              immediately save the class pickle fep file
           *autocompress*
              compress the xvg file with bzip2 (saves >70% of space)
        """

        from gromacs.formats import XVG

        if not stride is None:
            self.stride = stride

        if autocompress:
            # must be done before adding to results.xvg or we will not find the file later
            self.compress_dgdl_xvg()

        logger.info("[%(dirname)s] Finding dgdl xvg files, reading with "
                    "stride=%(stride)d permissive=%(permissive)r." % vars(self))

        for component, lambdas in self.lambdas.items():
            xvg_files = [self.dgdl_xvg(self.wdir(component, l)) for l in lambdas]
            self.results.xvg[component] = (numpy.array(lambdas),
                                           [XVG(xvg, permissive=self.permissive, stride=self.stride)
                                            for xvg in xvg_files])
        if autosave:
            self.save()

    def compress_dgdl_xvg(self):
        """Compress *all* dgdl xvg files with bzip2.

        .. Note:: After running this method you might want to run
           :meth:`collect` to ensure that the results in :attr:`results.xvg`
           point to the *compressed* files. Otherwise :exc:`IOError` might
           occur which fail to find a `md.xvg` file.

        """
        for component, lambdas in self.lambdas.items():
            xvg_files = [self.dgdl_xvg(self.wdir(component, l)) for l in lambdas]
            for xvg in xvg_files:
                root,ext = os.path.splitext(xvg)
                if ext == os.path.extsep+"xvg":
                    fnbz2 = xvg + os.path.extsep + "bz2"
                    logger.info("[%s] Compressing dgdl file %r with bzip2", self.dirname, xvg)
                    # speed is similar to 'bzip2 -9 FILE' (using a 1 Mio buffer)
                    with open(xvg, 'r', buffering=1048576) as source:
                        with openany(fnbz2, 'w', buffering=1048576) as target:
                            target.writelines(source)
                    if os.path.exists(fnbz2) and os.path.exists(xvg):
                        os.unlink(xvg)
                    if not os.path.exists(fnbz2):
                        logger.error("[%s] Failed to compress %r --- mysterious!", self.dirname, fnbz2)
                    else:
                        logger.info("[%s] Compression complete: %r", self.dirname, fnbz2)


    def contains_corrupted_xvgs(self):
        """Check if any of the source datafiles had reported corrupted lines.

        :Returns: ``True`` if any of the xvg dgdl files were produced with the
                  permissive=True flag and had skipped lines.

        For debugging purposes, the number of corrupted lines are stored in
        :attr:`Gsolv._corrupted` as dicts of dicts with the component as
        primary and the lambda as secondary key.
        """
        from itertools import izip
        def _lencorrupted(xvg):
            try:
                return len(xvg.corrupted_lineno)
            except AttributeError:  # backwards compatible (pre gw 0.1.10 are always ok)
                return 0
            except TypeError:       # len(None): XVG.parse() has not been run yet
                return 0            # ... so we cannot conclude that it does contain bad ones
        corrupted = {}
        self._corrupted = {}        # debugging ...
        for component, (lambdas, xvgs) in self.results.xvg.items():
            corrupted[component] = numpy.any([(_lencorrupted(xvg) > 0) for xvg in xvgs])
            self._corrupted[component] = dict(((l, _lencorrupted(xvg)) for l,xvg in izip(lambdas, xvgs)))
        return numpy.any([x for x in corrupted.values()])

    def analyze(self, force=False, stride=None, autosave=True):
        """Extract dV/dl from output and calculate dG by TI.

        Thermodynamic integration (TI) is performed on the individual
        component window calculation (typically the Coulomb and the
        VDW part, :math:`\Delta A_{\mathrm{coul}}` and :math:`\Delta
        A_{\mathrm{vdW}}`). :math:`\Delta A_{\mathrm{coul}}` is the free
        energy component of discharging the molecule and :math:`\Delta
        A_{\mathrm{vdW}}` of decoupling (switching off LJ interactions
        with the environment). The free energy components must be
        interpreted in this way because we defined ``lambda=0`` as
        interaction switched on and ``lambda=1`` as switched off.

        .. math::
            \Delta A* &= -(\Delta A_{\mathrm{coul}} + \Delta A_{\mathrm{vdw}})\\

        Data are stored in :attr:`Gsolv.results`.

        The dV/dlambda graphs are integrated with the composite Simpson's rule
        (and if the number of datapoints are even, the first interval is
        evaluated with the trapezoidal rule); see :func:`scipy.integrate.simps`
        for details). Note that this implementation of Simpson's rule does not
        require equidistant spacing on the lambda axis.

        For the Coulomb part using Simpson's rule has been shown to produce
        more accurate results than the trapezoidal rule [Jorge2010]_.

        Errors are estimated from the errors of the individual <dV/dlambda>:

         1. The error of the mean <dV/dlambda> is calculated via the decay time
            of the fluctuation around the mean (see
            :attr:`gromacs.formats.XVG.error`).

         2. The error on the integral is calculated analytically via
            propagation of errors through Simpson's rule (with the
            approximation that all spacings are assumed to be equal; taken as
            the average over all spacings as implemented in
            :func:`numkit.integration.simps_error`).

        .. Note:: For the Coulomb part, which typically only contains about 5
           lambdas, it is recommended to have a odd number of lambda values to
           fully benefit from the higher accuracy of the integration scheme.

        .. [Jorge2010] M. Jorge, N.M. Garrido, A.J. Queimada, I.G. Economou,
                       and E.A. Macedo. Effect of the integration method on the
                       accuracy and computational efficiency of free energy
                       calculations using thermodynamic integration. Journal of
                       Chemical Theory and Computation, 6 (4):1018--1027,
                       2010. 10.1021/ct900661c.

        :Keywords:
          *force*
              reload raw data even though it is already loaded
          *stride*
              read data every *stride* lines, ``None`` uses the class default
          *autosave*
              save to the pickle file when results have been computed
        """

        import scipy.integrate
        import numkit.integration

        stride = stride or self.stride

        if force or not self.has_dVdl():
            self.collect(stride=stride, autosave=False)
        else:
            logger.info("Analyzing stored data.")

        # total free energy difference at const P (all simulations are done in NPT)
        GibbsFreeEnergy = QuantityWithError(0,0)

        for component, (lambdas, xvgs) in self.results.xvg.items():
            logger.info("[%s %s] Computing averages <dV/dl> and errors for %d lambda values.",
                        self.molecule, component, len(lambdas))
            # for TI just get the average dv/dl value (in array column 1; col 0 is the time)
            # (This can take a while if the XVG is now reading the array from disk first time)
            # Use XVG class properties: first data in column 0!
            Y = numpy.array([x.mean[0] for x in xvgs])
            stdY = numpy.array([x.std[0]  for x in xvgs])
            DY = numpy.array([x.error[0]  for x in xvgs])   # takes a while: computes correl.time
            tc = numpy.array([x.tc[0]  for x in xvgs])
            self.results.dvdl[component] = {'lambdas':lambdas, 'mean':Y, 'error':DY,
                                            'stddev':stdY, 'tcorrel':tc}
            # Combined Simpson rule integration:
            # even="last" because dV/dl is smoother at the beginning so using trapezoidal
            # integration there makes less of an error (one hopes...)
            a = scipy.integrate.simps(Y, x=lambdas, even='last')
            da = numkit.integration.simps_error(DY, x=lambdas, even='last')
            self.results.DeltaA[component] = QuantityWithError(a, da)
            GibbsFreeEnergy += self.results.DeltaA[component]  # error propagation is automagic!

        # hydration free energy Delta A = -(Delta A_coul + Delta A_vdw)
        GibbsFreeEnergy *= -1
        self.results.DeltaA.Gibbs = GibbsFreeEnergy

        if autosave:
            self.save()

        self.logger_DeltaA0()
        return self.results.DeltaA.Gibbs

    def write_DeltaA0(self, filename, mode='w'):
        """Write free energy components to a file.

        :Arguments:
            *filename*
                 name of the text file
            *mode*
                 'w' for overwrite or 'a' for append ['w']

        Format::
          .                 ----- kJ/mol ------
          molecule solvent  total  coulomb  vdw
        """
        with open(filename, mode) as tab:
            tab.write(self.summary() + '\n')

    def summary(self):
        """Return a string that summarizes the energetics.

        Each energy component is followed by its error.

        Format::

          .                 ------ kJ/mol -----
          molecule solvent  total  coulomb  vdw

        """
        fmt = "%-10s %-14s %+8.2f %8.2f  %+8.2f %8.2f  %+8.2f %8.2f"
        d = self.results.DeltaA
        return fmt % ((self.molecule, self.solvent_type) + \
            d.Gibbs.astuple() +  d.coulomb.astuple() + \
            d.vdw.astuple())

    def logger_DeltaA0(self):
        """Print the free energy contributions (errors in parentheses)."""
        if not 'DeltaA' in self.results or len(self.results.DeltaA) == 0:
            logger.info("No DeltaA free energies computed yet.")
            return

        logger.info("DeltaG0 = -(DeltaG_coul + DeltaG_vdw)")
        for component, energy in self.results.DeltaA.items():
            logger.info("[%s] %s solvation free energy (%s) %g (%.2f) kJ/mol",
                        self.molecule, self.solvent_type.capitalize(), component,
                        energy.value, energy.error)

    def has_dVdl(self):
        """Check if dV/dl data have already been collected.

        :Returns: ``True`` if the dV/dl data have bee aquired
                  (:meth:`Gsolv.collect`) for all FEP simulations.
        """
        try:
            if len(self.results.xvg) == 0:
                return False
        except AttributeError:
            return False
        return numpy.all(numpy.array([len(xvgs) for (lambdas,xvgs) in self.results.xvg.values()]) > 0)

    def plot(self, **kwargs):
        """Plot the TI data with error bars.

        Run :meth:`mdpow.fep.Gsolv.analyze` first.

        All *kwargs* are passed on to :func:`pylab.errorbar`.
        """
        from pylab import subplot, errorbar, xlabel, ylabel, legend, xlim, title

        kwargs.setdefault('color', 'black')
        kwargs.setdefault('capsize', 0)
        kwargs.setdefault('elinewidth', 2)

        try:
            if self.results.DeltaA.Gibbs is None or len(self.results.dvdl) == 0:
                raise KeyError
        except KeyError:
            logger.info("Data were not analyzed yet -- doing that now... patience.")
            self.analyze()

        dvdl = self.results.dvdl
        nplots = len(dvdl)
        for i, component in enumerate(numpy.sort(dvdl.keys())):  # stable plot order
            x,y,dy = [dvdl[component][k] for k in ('lambdas', 'mean', 'error')]
            iplot = i+1
            subplot(1, nplots, iplot)
            energy = self.results.DeltaA[component]
            label = r"$\Delta A^{\rm{%s}}_{\rm{%s}} = %.2f\pm%.2f$ kJ/mol" \
                    % (component, self.solvent_type, energy.value, energy.error)
            errorbar(x, y, yerr=dy, label=label, **kwargs)
            xlabel(r'$\lambda$')
            legend(loc='best')
            xlim(-0.05, 1.05)
        subplot(1, nplots, 1)
        ylabel(r'$dV/d\lambda$ in kJ/mol')
        title(r"Free energy difference $\Delta A^{0}_{\rm{%s}}$" % self.solvent_type)
        subplot(1, nplots, 2)
        title(r"for %s: $%.2f\pm%.2f$ kJ/mol" %
              ((self.molecule,) + self.results.DeltaA.Gibbs.astuple()))


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

    def __repr__(self):
        return "%s(filename=%r)" % (self.__class__.__name__, self.filename)

class Ghyd(Gsolv):
    """Sets up and analyses MD to obtain the hydration free energy of a solute."""
    solvent_default = "water"
    dirname_default = os.path.join(Gsolv.topdir_default, solvent_default)


class Gcyclo(Gsolv):
    solvent_default = "cyclohexane"
    dirname_default = os.path.join(Gsolv.topdir_default, solvent_default)


class Goct(Gsolv):
    """Sets up and analyses MD to obtain the solvation free energy of a solute in octanol.

    The *coulomb* lambda schedule is enhanced compared to water as the initial
    part of the dV/dl curve is quite sensitive. By adding two additional points
    we hope to reduce the overall error on the dis-charging free energy.
    """
    solvent_default = "octanol"
    dirname_default = os.path.join(Gsolv.topdir_default, solvent_default)

    schedules = {'coulomb':
                 FEPschedule(name='coulomb',
                             description="dis-charging vdw+q --> vdw",
                             label='Coul',
                             couple_lambda0='vdw-q', couple_lambda1='vdw',
                             sc_alpha=0,      # linear scaling for coulomb
                             #lambdas=[0, 0.25, 0.5, 0.75, 1.0],  # default
                             lambdas=[0, 0.125, 0.25, 0.375, 0.5, 0.75, 1.0],  # +0.125, 0.375 enhanced
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



def p_transfer(G1, G2, **kwargs):
    """Compute partition coefficient from two :class:`Gsolv` objects.

       The order determines the direction of transfer: G1 --> G2.

       E.g.: ``G1 = Gwat`` and ``G2 = Goct`` will compute a water-octanol
       transfer free energy and a water-octanol partition coefficient.

       transfer free energy from water into octanol::

            DeltaDeltaG0 = DeltaG0_oct - DeltaG0_water

       water-octanol partition coefficient::

            log P_oct/wat =  log [X]_oct/[X]_wat

    .. note::

       The partition coefficient corresponds to a transfer in the *opposite*
       direction compared to the free energy difference computed here. For
       instance, ``log P_oct/water`` is related to the free energy of transfer
       from octanol to water, ``Delta G_water - Delta G_oct``. ``log
       P_oct/water < 0`` means that the water phase is preferred, i.e. the
       probability to find the solute in octanol is *smaller* than the
       probability to find it in the water phase, or in other words, the
       hydration free energy ``Delta G_water`` is smaller (more negative) than
       the octanol solvation free energy ``Delta G_oct``.


    :Arguments:
       *G1*, *G2*
           *G1* and *G2* should be two  :class:`Gsolv` instances,
           order matters.
       *force*
           force rereading of data files even if some data were already stored [False]
       *stride*
           analyze every *stride*-th datapoint in the dV/dlambda files
    :Returns: (transfer free energy, log10 of the water-octanol partition coefficient = log Pow)

    """

    kwargs.setdefault('force', False)

    if G1.molecule != G2.molecule:
        raise ValueError("The two simulations were done for different molecules.")
    if G1.Temperature != G2.Temperature:
        raise ValueError("The two simulations were done at different temperatures.")

    logger.info("[%s] transfer free energy %s --> %s calculation",
                G1.molecule, G1.solvent_type, G2.solvent_type)
    for G in (G1, G2):
        if kwargs['force']:
            G.analyze(**kwargs)
        try:
            G.results.DeltaA.Gibbs
            G.logger_DeltaA0()
        except (KeyError, AttributeError):   # KeyError because results is a AttributeDict
            logger.warn("Must analyze simulation because no hydration free energy data available...")
            G.analyze(**kwargs)

    # x.Gibbs are QuantityWithError so they do error propagation
    transferFE = G2.results.DeltaA.Gibbs - G1.results.DeltaA.Gibbs
    # note minus sign, with our convention of the free energy difference to be
    # opposite to the partition coefficient.
    logPow = -transferFE / (kBOLTZ * G1.Temperature) * numpy.log10(numpy.e)

    molecule = G1.molecule

    # lower case initials, in reverse order of transfer, e.g.
    # water -> octanol:      P_ow
    # water -> cyclohexane:  P_cw
    coefficient = "P_{0}{1}".format(
        G2.solvent_type.lower()[0], G1.solvent_type.lower()[0])

    logger.info("[%s] Values at T = %g K", molecule, G1.Temperature)
    logger.info("[%s] Free energy of transfer %s --> %s: %.3f (%.3f) kJ/mol",
                molecule,
                G1.solvent_type, G2.solvent_type,
                transferFE.value, transferFE.error)
    logger.info("[%s] log %s: %.3f (%.3f)",
                molecule, coefficient, logPow.value, logPow.error)

    return transferFE, logPow

def pOW(G1, G2, **kwargs):
    """Compute water-octanol partition coefficient from two :class:`Gsolv` objects.

    transfer free energy from water into octanol::

            DeltaDeltaG0 = DeltaG0_oct - DeltaG0_water

    octanol/water partition coefficient::

            log P_oct/wat =  log [X]_oct/[X]_wat

    :Arguments:
       *G1*, *G2*
           *G1* and *G2* should be a :class:`Ghyd` and a :class:`Goct` instance,
           but order does not matter
       *force*
           force rereading of data files even if some data were already stored [False]
       *stride*
           analyze every *stride*-th datapoint in the dV/dlambda files

    :Returns: (transfer free energy, log10 of the octanol/water partition coefficient = log Pow)
    """
    if G1.solvent_type == "water" and G2.solvent_type == "octanol":
        args = (G1, G2)
    elif G1.solvent_type == "octanol" and G2.solvent_type == "water":
        args = (G2, G1)
    else:
        msg = "For pOW need water and octanol simulations but instead got {0} and {1}".format(
            G1.solvent_type, G2.solvent_type)
        logger.error(msg)
        raise ValueError(msg)
    return p_transfer(*args, **kwargs)

def pCW(G1, G2, **kwargs):
    """Compute water-cyclohexane partition coefficient from two :class:`Gsolv` objects.

    transfer free energy from water into cyclohexane::

            DeltaDeltaG0 = DeltaG0_cyclohexane - DeltaG0_water

    cyclohexane/water partition coefficient::

            log P_CW =  log [X]_cyclohexane/[X]_water

    :Arguments:
       *G1*, *G2*
           *G1* and *G2* should be a :class:`Ghyd` and a :class:`Gcyclo` instance,
           but order does not matter
       *force*
           force rereading of data files even if some data were already stored [False]
       *stride*
           analyze every *stride*-th datapoint in the dV/dlambda files

    :Returns: (transfer free energy, log10 of the cyclohexane/water partition coefficient = log Pcw)
    """
    if G1.solvent_type == "water" and G2.solvent_type == "cyclohexane":
        args = (G1, G2)
    elif G1.solvent_type == "cyclohexane" and G2.solvent_type == "water":
        args = (G2, G1)
    else:
        msg = "For pCW need water and cyclohexane simulations but instead got {0} and {1}".format(
            G1.solvent_type, G2.solvent_type)
        logger.error(msg)
        raise ValueError(msg)
    return p_transfer(*args, **kwargs)
