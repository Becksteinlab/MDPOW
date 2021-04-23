# equil.py
# Copyright (c) 2010-2011 Oliver Beckstein

"""
:mod:`mdpow.equil` --- Setting up and running equilibrium MD
============================================================

The :mod:`mdpow.equil` module facilitates the setup of equilibrium
molecular dynamics simulations of a compound molecule in a simulation
box of water or other solvent such as octanol.

It requires as input

 - the itp file for the compound
 - a coordinate (structure) file (in pdb or gro format)

By default it uses the *OPLS/AA* forcefield and the *TIP4P* water
model.

.. warning:: Other forcefields than OPLS/AA are currently not
   officially supported; it is not hard to do but requires tedious
   changes to a few paths in template scripts.

.. autoclass:: Simulation
   :members:
.. autoclass:: WaterSimulation
.. autoclass:: OctanolSimulation

.. autodata:: DIST
"""

from __future__ import absolute_import, with_statement

import os, errno
import shutil
import cPickle
import MDAnalysis as mda

try:
    import gromacs.setup, gromacs.cbook
except (ImportError, OSError):
    raise ImportError("Gromacs installation not found, source GMXRC?")
from gromacs.utilities import in_dir, realpath, asiterable, AttributeDict
import gromacs.utilities

from . import config
from . import forcefields
from .restart import Journalled

import logging
logger = logging.getLogger('mdpow.equil')

# ITP <-- forcefields.get_solvent_model(id).itp
# BOX <-- forcefields.get_solvent_model(id).coordinates


# TODO: change to water distance 1.2 in the future (1.0 for
#       compatibility with our SAMPL5 runs)
#: minimum distance between solute and box surface (in nm)
DIST = {'water': 1.0, 'octanol': 1.5, 'cyclohexane': 1.5, 'wetoctanol': 1.5}

class Simulation(Journalled):
    """Simple MD simulation of a single compound molecule in water.

    Typical use ::

       S = Simulation(molecule='DRUG')
       S.topology(itp='drug.itp')
       S.solvate(struct='DRUG-H.pdb')
       S.energy_minimize()
       S.MD_relaxed()
       S.MD()

    .. Note:: The OPLS/AA force field and the TIP4P water molecule is the
              default; changing this is possible but will require provision of
              customized itp, mdp and template top files at various stages.
    """

    #: Keyword arguments to pre-set some file names; they are keys in :attr:`Simulation.files`.
    filekeys = ('topology', 'processed_topology', 'structure', 'solvated', 'ndx',
                'energy_minimized', 'MD_relaxed', 'MD_restrained', 'MD_NPT')
    topdir_default = "Equilibrium"
    dirname_default = os.path.curdir
    solvent_default = 'water'

    #: Coordinate files of the full system in increasing order of advancement of
    #: the protocol; the later the better. The values are keys into :attr:`Simulation.files`.
    coordinate_structures = ('solvated', 'energy_minimized', 'MD_relaxed',
                             'MD_restrained', 'MD_NPT')
    checkpoints = ('solvated','energy_minimized','MD_relaxed','MD_restrained','MD_NPT')


    #: Check list of all methods that can be run as an independent protocol; see also
    #: :meth:`Simulation.get_protocol` and :class:`restart.Journal`
    protocols = ("MD_NPT", "MD_NPT_run",                 # *_run as dummies for the ...
                 "MD_relaxed", "MD_relaxed_run",         # ...checkpointing logic
                 "MD_restrained", "MD_restrained_run",
                 "energy_minimize", "solvate", "topology")

    #: Default Gromacs *MDP* run parameter files for the different stages.
    #: (All are part of the package and are found with :func:`mdpow.config.get_template`.)
    mdp_defaults = {'MD_relaxed': 'NPT_opls.mdp',
                    'MD_restrained': 'NPT_opls.mdp',
                    'MD_NPT': 'NPT_opls.mdp',
                    'energy_minimize': 'em_opls.mdp',
                    }

    def __init__(self, molecule=None, **kwargs):
        """Set up Simulation instance.

        The *molecule* of the compound molecule should be supplied. Existing files
        (which have been generated in previous runs) can also be supplied.

        :Keywords:
          *molecule*
              Identifier for the compound molecule. This is the same as the
              entry in the ``[ molecule ]`` section of the itp file. ["DRUG"]
          *filename*
              If provided and *molecule* is ``None`` then load the instance from
              the pickle file *filename*, which was generated with
              :meth:`~mdpow.equil.Simulation.save`.
          *dirname*
              base directory; all other directories are created under it
          *forcefield*
              'OPLS-AA' or 'CHARMM' or 'AMBER'
          *solvent*
              'water' or 'octanol' or 'cyclohexane' or 'wetoctanol'
          *solventmodel*
              ``None`` chooses the default (e.g, :data:`mdpow.forcefields.DEFAULT_WATER_MODEL`
              for ``solvent == "water"``. Other options are the models defined in
              :data:`mdpow.forcefields.GROMACS_WATER_MODELS`. At the moment, there are no
              alternative parameterizations included for other solvents.
          *mdp*
              dict with keys corresponding to the stages ``energy_minimize``,
              ``MD_restrained``, ``MD_relaxed``,
              ``MD_NPT`` and values *mdp* file names (if no entry then the
              package defaults are used)
          *distance*
               minimum distance between solute and closest box face
          *kwargs*
              advanced keywords for short-circuiting; see
              :data:`mdpow.equil.Simulation.filekeys`.

        """
        self.__cache = {}
        filename = kwargs.pop('filename', None)
        dirname = kwargs.pop('dirname', self.dirname_default)

        forcefield = kwargs.pop('forcefield', 'OPLS-AA')
        solvent = kwargs.pop('solvent', self.solvent_default)
        # mdp files --- should get values from default runinput.cfg
        # None values in the kwarg mdp dict are ignored
        # self.mdp: key = stage, value = path to MDP file

        # 'water' will choose the default ('tip4p'), other choices are
        # 'tip3p', 'spc', 'spce', 'm24', for water; no choices
        # available for 'cyclohexane' and 'octanol'
        solventmodel = kwargs.pop('solventmodel', None)

        mdp_kw = kwargs.pop('mdp', {})
        self.mdp = dict((stage, config.get_template(fn)) for stage,fn in self.mdp_defaults.items())
        self.mdp.update(dict((stage, config.get_template(fn)) for stage,fn in mdp_kw.items() if fn is not None))

        if molecule is None and filename is not None:
            # load from pickle file
            self.load(filename)
            self.filename = filename
            kwargs = {}    # for super
        else:
            self.molecule = molecule or 'DRUG'
            self.dirs = AttributeDict(
                basedir=realpath(dirname),    # .../Equilibrium/<solvent>
                includes=list(asiterable(kwargs.pop('includes',[]))) + [config.includedir],
                )
            # pre-set filenames: keyword == variable name
            self.files = AttributeDict([(k, kwargs.pop(k, None)) for k in self.filekeys])
            self.deffnm = kwargs.pop("deffnm", "md")

            if self.files.topology:
                # assume that a user-supplied topology lives in a 'standard' top dir
                # that includes the necessary itp file(s)
                self.dirs.topology = realpath(os.path.dirname(self.files.topology))
                self.dirs.includes.append(self.dirs.topology)

            self.forcefield = forcefield
            self.solvent_type = solvent
            self.solventmodel_identifier = forcefields.get_solvent_identifier(
                 solvent,
                 model=solventmodel,
                 forcefield=forcefield,
                 )
            if self.solventmodel_identifier is None:
                msg = "No parameters for solvent {0} and solventmodel {1} available.".format(
                    solvent, solventmodel)
                logger.error(msg)
                raise ValueError(msg)
            self.solventmodel = forcefields.get_solvent_model(
                self.solventmodel_identifier,
                forcefield=forcefield,
                )

            distance = kwargs.pop('distance', None)
            distance = distance if distance is not None else DIST[solvent]

            self.solvent = AttributeDict(itp=self.solventmodel.itp,
                                         box=self.solventmodel.coordinates,
                                         distance=distance)

            self.filename = filename or self.solvent_type+'.simulation'

        super(Simulation, self).__init__(**kwargs)

    def BASEDIR(self, *args):
        return os.path.join(self.dirs.basedir, *args)

    def save(self, filename=None):
        """Save instance to a pickle file.

        The default filename is the name of the file that was last loaded from
        or saved to.
        """
        if filename is None:
            if self.filename is None:
                self.filename = filename or self.solvent_type+'.simulation'
                logger.warning("No filename known, saving instance under name %r", self.filename)
            filename = self.filename
        else:
            self.filename = filename
        with open(filename, 'wb') as f:
            cPickle.dump(self, f, protocol=cPickle.HIGHEST_PROTOCOL)
        logger.debug("Instance pickled to %(filename)r" % vars())

    def load(self, filename=None):
        """Re-instantiate class from pickled file."""
        if filename is None:
            if self.filename is None:
                self.filename = self.molecule.lower() + '.pickle'
                logger.warning("No filename known, trying name %r", self.filename)
            filename = self.filename
        with open(filename, 'rb') as f:
            instance = cPickle.load(f)
        self.__dict__.update(instance.__dict__)
        logger.debug("Instance loaded from %(filename)r" % vars())

    def make_paths_relative(self, prefix=os.path.curdir):
        """Hack to be able to copy directories around: prune basedir from paths.

        .. Warning:: This is not guaranteed to work for all paths. In particular,
                     check :attr:`mdpow.equil.Simulation.dirs.includes` and adjust
                     manually if necessary.
        """
        def assinglet(m):
            if len(m) == 1:
                return m[0]
            elif len(m) == 0:
                return None
            return m

        basedir = self.dirs.basedir
        for key, fn in self.files.items():
            try:
                self.files[key] = fn.replace(basedir, prefix)
            except AttributeError:
                pass
        for key, val in self.dirs.items():
            fns = asiterable(val)  # treat them all as lists
            try:
                self.dirs[key] = assinglet([fn.replace(basedir, prefix) for fn in fns])
            except AttributeError:
                pass
        for key, fn in self.mdp.items():
            try:
                self.mdp[key] = fn.replace(basedir, prefix)
            except AttributeError:
                pass
        logger.warn("make_paths_relative(): check/manually adjust %s.dirs.includes = %r !",
                    self.__class__.__name__, self.dirs.includes)

    def topology(self, itp='drug.itp', prm=None, **kwargs):
        """Generate a topology for compound *molecule*.

        :Keywords:
            *itp*
               Gromacs itp file; will be copied to topology dir and
               included in topology
            *prm*
               Gromacs prm file; if given, will be copied to topology
               dir and included in topology
            *dirname*
               name of the topology directory ["top"]
            *kwargs*
               see source for *top_template*, *topol*
        """
        self.journal.start('topology')

        dirname = kwargs.pop('dirname', self.BASEDIR('top'))
        self.dirs.topology = realpath(dirname)

        setting = forcefields.get_ff_paths(self.forcefield)
        template = forcefields.get_top_template(self.solvent_type)

        top_template = config.get_template(kwargs.pop('top_template', template))
        topol = kwargs.pop('topol', os.path.basename(top_template))
        self.top_template = top_template
        itp = os.path.realpath(itp)
        _itp = os.path.basename(itp)

        if prm==None:
            prm_kw = ''
        else:
            prm = os.path.realpath(prm)
            _prm = os.path.basename(prm)
            prm_kw = '#include "{}"'.format(_prm)

        with in_dir(dirname):
            shutil.copy(itp, _itp)
            if prm is not None:
                shutil.copy(prm, _prm)
            gromacs.cbook.edit_txt(top_template,
                                   [('#include +"oplsaa\.ff/forcefield\.itp"',
                                    'oplsaa\.ff/', setting[0]),
                                    ('#include +"compound\.itp"', 'compound\.itp', _itp),
                                    ('#include +"oplsaa\.ff/tip4p\.itp"',
                                     'oplsaa\.ff/tip4p\.itp', setting[0]+self.solvent.itp),
                                    ('#include +"oplsaa\.ff/ions_opls\.itp"',
                                     'oplsaa\.ff/ions_opls\.itp', setting[1]),
                                    ('#include +"compound\.prm"',
                                     '#include +"compound\.prm"', prm_kw),
                                    ('#include +"water\.itp"', 'water\.itp', setting[2]),
                                    ('Compound', 'solvent', self.solvent_type),
                                    ('Compound', 'DRUG', self.molecule),
                                    ('DRUG\s*1', 'DRUG', self.molecule),
                                    ],
                                   newname=topol)
        logger.info('[%(dirname)s] Created topology %(topol)r that includes %(_itp)r', vars())

        # update known files and dirs
        self.files.topology = realpath(dirname, topol)
        if not self.dirs.topology in self.dirs.includes:
            self.dirs.includes.append(self.dirs.topology)

        self.journal.completed('topology')
        return {'dirname': dirname, 'topol': topol}

    @staticmethod
    def _setup_solvate(**kwargs):
        """Solvate structure in a single solvent box."""
        return gromacs.setup.solvate(**kwargs)

    def solvate(self, struct=None, **kwargs):
        """Solvate structure *struct* in a box of solvent.

        The solvent is determined with the *solvent* keyword to the constructor.

        :Keywords:
          *struct*
              pdb or gro coordinate file (if not supplied, the value is used
              that was supplied to the constructor of :class:`~mdpow.equil.Simulation`)
          *distance*
               minimum distance between solute and the closes box face; the default depends
               on the solvent but can be set explicitly here, too.
          *bt*
               any box type understood by :func:`gromacs.editconf` (``-bt``):

               * "triclinic" is a triclinic box,
               * "cubic" is a rectangular box with all sides equal;
               * "dodecahedron" represents a rhombic dodecahedron;
               * "octahedron" is a truncated octahedron.

               The default is "dodecahedron".
          *kwargs*
              All other arguments are passed on to :func:`gromacs.setup.solvate`, but
              set to sensible default values. *top* and *water* are always fixed.
        """
        self.journal.start('solvate')

        self.dirs.solvation = realpath(kwargs.setdefault('dirname', self.BASEDIR('solvation')))
        kwargs['struct'] = self._checknotempty(struct or self.files.structure, 'struct')
        kwargs['top'] = self._checknotempty(self.files.topology, 'top')
        kwargs['water'] = self.solvent.box
        kwargs.setdefault('mainselection', '"%s"' % self.molecule)  # quotes are needed for make_ndx
        kwargs.setdefault('distance', self.solvent.distance)

        boxtype = kwargs.pop('bt', None)
        boxtype = boxtype if boxtype is not None else "dodecahedron"
        if boxtype not in ("dodecahedron", "triclinic", "cubic", "octahedron"):
            msg = "Invalid boxtype '{0}', not suitable for 'gmx editconf'.".format(boxtype)
            logger.error(msg)
            raise ValueError(msg)
        kwargs['bt'] = boxtype

        kwargs['includes'] = asiterable(kwargs.pop('includes',[])) + self.dirs.includes

        params = self._setup_solvate(**kwargs)

        self.files.structure = kwargs['struct']
        self.files.solvated = params['struct']
        self.files.ndx = params['ndx']
        # we can also make a processed topology right now
        self.processed_topology(**kwargs)

        self.journal.completed('solvate')
        return params

    def processed_topology(self, **kwargs):
        """Create a portable topology file from the topology and the solvated system."""
        if self.files.solvated is None or not os.path.exists(self.files.solvated):
            self.solvate(**kwargs)
        kwargs['topol'] = self.files.topology
        kwargs['struct'] = self.files.solvated
        kwargs['includes'] = self.dirs.includes
        self.files.processed_topology = gromacs.cbook.create_portable_topology(**kwargs)
        return self.files.processed_topology

    def energy_minimize(self, **kwargs):
        """Energy minimize the solvated structure on the local machine.

        *kwargs* are passed to :func:`gromacs.setup.energ_minimize` but if
        :meth:`~mdpow.equil.Simulation.solvate` step has been carried out
        previously all the defaults should just work.
        """
        self.journal.start('energy_minimize')

        self.dirs.energy_minimization = realpath(kwargs.setdefault('dirname', self.BASEDIR('em')))
        kwargs['top'] = self.files.topology
        kwargs.setdefault('struct', self.files.solvated)
        kwargs.setdefault('mdp', self.mdp['energy_minimize'])
        kwargs['mainselection'] = None
        kwargs['includes'] = asiterable(kwargs.pop('includes',[])) + self.dirs.includes

        params = gromacs.setup.energy_minimize(**kwargs)

        self.files.energy_minimized = params['struct']

        self.journal.completed('energy_minimize')
        return params

    def _MD(self, protocol, **kwargs):
        """Basic MD driver for this Simulation. Do not call directly."""
        self.journal.start(protocol)

        kwargs.setdefault('dirname', self.BASEDIR(protocol))
        kwargs.setdefault('deffnm', self.deffnm)
        kwargs.setdefault('mdp', config.get_template('NPT_opls.mdp'))
        self.dirs[protocol] = realpath(kwargs['dirname'])
        setupMD = kwargs.pop('MDfunc', gromacs.setup.MD)
        kwargs['top'] = self.files.topology
        kwargs['includes'] = asiterable(kwargs.pop('includes',[])) + self.dirs.includes
        kwargs['ndx'] = self.files.ndx
        kwargs['mainselection'] = None # important for SD (use custom mdp and ndx!, gromacs.setup._MD)
        self._checknotempty(kwargs['struct'], 'struct')
        if not os.path.exists(kwargs['struct']):
            # struct is not reliable as it depends on qscript so now we just try everything...
            struct = gromacs.utilities.find_first(kwargs['struct'], suffices=['pdb', 'gro'])
            if struct is None:
                logger.error("Starting structure %(struct)r does not exist (yet)" % kwargs)
                raise IOError(errno.ENOENT, "Starting structure not found", kwargs['struct'])
            else:
                logger.info("Found starting structure %r (instead of %r).", struct, kwargs['struct'])
                kwargs['struct'] = struct
        # now setup the whole simulation (this is typically gromacs.setup.MD() )
        params =  setupMD(**kwargs)
        # params['struct'] is md.gro but could also be md.pdb --- depends entirely on qscript
        self.files[protocol] = params['struct']
        # Gromacs 4.5.x 'mdrun -c PDB'  fails if it cannot find 'residuetypes.dat'
        # so instead of fuffing with GMXLIB we just dump it into the directory
        try:
            shutil.copy(config.topfiles['residuetypes.dat'], self.dirs[protocol])
        except:
            logger.warn("Failed to copy 'residuetypes.dat': mdrun will likely fail to write a final structure")

        self.journal.completed(protocol)
        return params

    def MD_relaxed(self, **kwargs):
        """Short MD simulation with *timestep* = 0.1 fs to relax strain.

        Energy minimization does not always remove all problems and LINCS
        constraint errors occur. A very short *runtime* = 5 ps MD with very
        short integration time step *dt* tends to solve these problems.

        .. See Also:: :func:`gromacs.setup.MD`

        :Keywords:
          *struct*
             starting coordinates (typically guessed)
          *mdp*
             MDP run parameter file for Gromacs
          *qscript*
             list of queuing system submission scripts; probably a
             good idea to always include the default "local.sh" even
             if you have your own ["local.sh"]
          *qname*
             name of the job as shown in the queuing system
          *startdir*
             **advanced uses**: path of the directory on a remote
             system, which will be hard-coded into the queuing system
             script(s); see :func:`gromacs.setup.MD` and
             :class:`gromacs.manager.Manager`

        """
        # user structure or restrained or solvated
        kwargs.setdefault('struct', self.files.energy_minimized)
        kwargs.setdefault('dt', 0.0001)  # ps
        kwargs.setdefault('runtime', 5)  # ps
        kwargs.setdefault('mdp', self.mdp['MD_relaxed'])
        return self._MD('MD_relaxed', **kwargs)

    def MD_restrained(self, **kwargs):
        """Short MD simulation with position restraints on compound.

        See documentation of :func:`gromacs.setup.MD_restrained` for
        details. The following keywords can not be changed: top, mdp, ndx,
        mainselection

        .. Note:: Position restraints are activated with ``-DPOSRES`` directives
                  for :func:`gromacs.grompp`. Hence this will only work if the
                  compound itp file does indeed contain a ``[ posres ]``
                  section that is protected by a ``#ifdef POSRES`` clause.

        .. See Also:: :func:`gromacs.setup.MD_restrained`

        :Keywords:
          *struct*
              starting coordinates (leave empty for inspired guess of file name)
          *mdp*
              MDP run parameter file for Gromacs
          *qscript*
             list of queuing system submission scripts; probably a
             good idea to always include the default "local.sh" even
             if you have your own ["local.sh"]
          *qname*
             name of the job as shown in the queuing system
          *startdir*
             **advanced uses**: path of the directory on a remote
             system, which will be hard-coded into the queuing system
             script(s); see :func:`gromacs.setup.MD` and
             :class:`gromacs.manager.Manager`

        """
        kwargs.setdefault('struct',
                          self._lastnotempty([self.files.energy_minimized, self.files.MD_relaxed]))
        kwargs.setdefault('mdp', self.mdp['MD_restrained'])
        kwargs['MDfunc'] = gromacs.setup.MD_restrained
        return self._MD('MD_restrained', **kwargs)

    def MD_NPT(self, **kwargs):
        """Short NPT MD simulation.

        See documentation of :func:`gromacs.setup.MD` for details such
        as *runtime* or specific queuing system options. The following
        keywords can not be changed: *top*, *mdp*, *ndx*, *mainselection*.

        .. Note:: If the system crashes (with LINCS errors), try initial
                  equilibration with timestep *dt* = 0.0001 ps (0.1 fs instead
                  of 2 fs) and *runtime* = 5 ps as done in :meth:`~Simulation.MD_relaxed`

        .. See Also:: :func:`gromacs.setup.MD` and :meth:`Simulation.MD_relaxed`

        :Keywords:
          *struct*
               starting conformation; by default, the *struct* is the last frame
               from the position restraints run, or, if this file cannot be
               found (e.g. because :meth:`Simulation.MD_restrained` was not run)
               it falls back to the relaxed and then the solvated system.
          *mdp*
               MDP run parameter file for Gromacs
          *runtime*
               total run time in ps
          *qscript*
               list of queuing system scripts to prepare; available values are
               in :data:`gromacs.config.templates` or you can provide your own
               filename(s) in the current directory (see :mod:`gromacs.qsub` for
               the format of the templates)
          *qname*
             name of the job as shown in the queuing system
          *startdir*
             **advanced uses**: path of the directory on a remote
             system, which will be hard-coded into the queuing system
             script(s); see :func:`gromacs.setup.MD` and
             :class:`gromacs.manager.Manager`

        """
        # user structure or relaxed or restrained or solvated
        kwargs.setdefault('struct', self.get_last_structure())
        kwargs.setdefault('t',self.get_last_checkpoint()) # Pass checkpoint file from md_relaxed
        kwargs.setdefault('mdp', self.mdp['MD_NPT'])
        return self._MD('MD_NPT', **kwargs)

    # for convenience and compatibility
    MD = MD_NPT

    @staticmethod
    def _checknotempty(value, name):
        if value is None or value == "":
            raise ValueError("Parameter %s cannot be empty." % name)
        return value

    @staticmethod
    def _lastnotempty(l):
        """Return the last non-empty value in list *l* (or None :-p)"""
        nonempty = [None] + [x for x in l if not (x is None or x == "" or x == [])]
        return nonempty[-1]

    def get_last_structure(self):
        """Returns the coordinates of the most advanced step in the protocol."""
        return self._lastnotempty([self.files[name] for name in self.coordinate_structures])

    def get_last_checkpoint(self):
        """Returns the checkpoint of the most advanced step in the protocol.
        Relies on md.gro being present from previous simulation, assumes that checkpoint is then present.
        """
        return self._lastnotempty([self.files[name] for name in self.checkpoints]).replace('.gro','.cpt')

class WaterSimulation(Simulation):
    """Equilibrium MD of a solute in a box of water."""
    solvent_default = 'water'
    dirname_default = os.path.join(Simulation.topdir_default, solvent_default)

class CyclohexaneSimulation(Simulation):
    """Equilibrium MD of a solute in a box of cyclohexane."""
    solvent_default = 'cyclohexane'
    dirname_default = os.path.join(Simulation.topdir_default, solvent_default)

class OctanolSimulation(Simulation):
    """Equilibrium MD of a solute in a box of octanol."""
    solvent_default = 'octanol'
    dirname_default = os.path.join(Simulation.topdir_default, solvent_default)

class WetOctanolSimulation(Simulation):
    """Equilibrium MD of a solute in a box of wet octanol."""
    solvent_default = 'wetoctanol'
    dirname_default = os.path.join(Simulation.topdir_default, solvent_default)

    def  _setup_solvate(self, **kwargs):
        sol = gromacs.setup.solvate_sol(**kwargs)
        with in_dir(self.dirs.solvation, create=False):
            u = mda.Universe('solvated.gro')
            octanol = u.select_atoms('resname OcOH')
            n = octanol.n_residues
        with in_dir(self.dirs.topology, create=False):
            gromacs.cbook.edit_txt(self.files.topology,
                                   [('OcOH               1', '1', n)])
        ionkwargs = kwargs
        ionkwargs['struct'] = sol['struct']
        params = gromacs.setup.solvate_ion(**ionkwargs)
        return params
