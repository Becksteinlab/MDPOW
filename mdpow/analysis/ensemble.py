# MDPOW: ensemble.py
# 2021 Alia Lescoulie

import os
import errno
from typing import Optional, List

import numpy as np

import MDAnalysis as mda
from MDAnalysis.lib.log import ProgressBar
from MDAnalysis.exceptions import FileFormatWarning, NoDataError, MissingDataWarning, SelectionError

from gromacs.utilities import in_dir

import logging

logger = logging.getLogger('mdpow._ensemble')


class Ensemble(object):
    """ Collection of related :class:`MDAnalysis.Universe <MDAnalysis.core.groups.universe.Universe>`
    objects.

    Stores systems produced by running mdpow-fep organized
    by solvent, interaction, and lambda.

    Given a mdpow simulation directory will load the MD
    simulation files with the directory structure as keys.

    :Keywords:
    *dirname*
        Molecule Simulation directory. Loads simulation files present in
        lambda directories into the new instance. With this method for
        generating an :class:`~mdpow.analysis.ensemble.Ensemble` the
        lambda directories are explored and
        :meth:`~mdpow.analysis.ensemble.Ensemble._load_universe_from_dir`
        searches for .gro, .gro.b2z, .gro.gz, and .tpr files for topology,
        and .xtc files for trajectory. It will default to using the tpr file
        available.

    *solvents*
        Solvents from directory given to the new instance. Default
        :code:`solvents=('water', 'octanol')`

    *topology_paths*
        Specifies topologies used in loading simulated systems. Given
        with a dictionary with keys-value pair for each solvent and
        its respective topology path.

    *interactions*
        Interactions from directory given to the instance. Default
        :code:`interactions=('Coulomb', 'VDW')`

    *universe_kwargs*
        `Keywords arguments <https://docs.mdanalysis.org/stable/documentation_pages/core/universe>`_
        for loading :class:`MDAnalysis.Universe <MDAnalysis.core.groups.universe.Universe>`
        objects from MDPOW files in :code:`dirname` argument directory when creating an
        :class:`~mdpow.analysis.ensemble.Ensemble` .

    .. rubric:: Examples


    Typical workflow for MDPOW directory::

        ens = Ensemble(dirname='molecule')

    Typical workflow for adding universes individually::

        ens = Ensemble()
        u = mda.Universe(md.gro', 'md.xtc')
        ens.add_system(u)

    Topology paths can be specified when defining the _ensemble
    by giving the paths to each solvent topology in a dictionary
    with the topology_paths argument::

        ens = Ensemble(dirname='molecule', topology_paths={'water': water_path,
                                                           'octanol': octanol_path}

    Interactions can also be specified when initializing the with
    a list using the interactions argument::

        ens = Ensemble(dirname='molecule', interactions=['Coulomb']

    .. versionadded:: 0.8.0
    """

    def __init__(self, dirname=None, solvents=('octanol', 'water'),
                 topology_paths=None, interactions=('Coulomb', 'VDW'),
                 **universe_kwargs):
        self.top_dict = topology_paths
        self._num_systems = 0
        self._ensemble = {}
        self._keys = []
        self._solvents = solvents
        self._interactions = interactions
        self.unv_kwargs = universe_kwargs

        if dirname is None:
            return

        if not os.path.exists(dirname):
            logger.error(f"Directory {dirname} does not exist")
            raise FileNotFoundError(errno.ENOENT, 'Directory does not'
                                                  'exist', dirname)

        self._ensemble_dir = dirname
        self._build_ensemble()

    def __repr__(self):
        return f"<Ensemble Containing {self._num_systems} System>"

    def __len__(self):
        return self._num_systems

    def __getitem__(self, index):
        """Allows dictionary like indexing"""
        return self._ensemble[index]

    @staticmethod
    def _load_universe_from_dir(solv_dir=None, **universe_kwargs) -> Optional[mda.Universe]:
        """Loads system simulation files in directory into an
        :class:`MDAnalysis.Universe <MDAnalysis.core.groups.universe.Universe>`


        If multiple topologies are found it will default to using
        the .tpr file. If more than one trajectory is present they
        will be sorted alphabetically and passed into the
        :class:`MDAnalysis.Universe <MDAnalysis.core.groups.universe.Universe>`
        This method is run automatically by :meth:`~mdpow/analysis/ensemble.Ensemble._build_ensemble`
        when initializing the class using the :code:`dirname` argument.
        """

        def _sort_trajectories(trajectories: list) -> list:
            """Sorts list of trajectory files alphabetically and makes paths
            absolute"""
            sorted(trajectories)
            return [os.path.abspath(t) for t in trajectories]

        def _sort_topologies(topologies: list) -> list:
            """sorts list of trajectory files with .tpr first"""
            tops = []
            logger.info('If more than one topology is present the tpr will be the one used')
            for i in range(len(topologies)):
                f = topologies[i]
                if f.endswith('.tpr'):
                    topologies.pop(i)
                    tops = [f] + topologies
                    break
            return tops

        cur_dir = os.listdir(os.curdir)
        trj = []
        top = []

        if solv_dir is not None:
            # if top is specified in kwargs, saved to list
            top = [solv_dir]

        for file in cur_dir:
            if file.endswith('.xtc'):
                # Saving trajectory directories
                trj.append(file)
            elif (file.endswith('gro') or file.endswith('.tpr') or file.endswith('gro.bz2')
                  or file.endswith('gro.gz')) and solv_dir is None:
                # Saving topology directories
                top.append(file)

        if len(top) == 0 or len(trj) == 0:
            logger.warning('No MD files detected in %s', os.curdir)
            return

        trj = _sort_trajectories(trj)
        if len(top) > 1:
            top = _sort_topologies(top)

        try:
            return mda.Universe(os.path.abspath(top[0]), trj, **universe_kwargs)
        except (ValueError, FileFormatWarning, NoDataError, MissingDataWarning, OSError) as err:
            logger.error(f'{err} raised while loading {top[0]} {trj} in dir {cur_dir}')
            raise NoDataError

    def keys(self):
        """Returns list of system keys"""
        return self._keys

    def _build_ensemble(self):
        """Finds simulation files genderated by MDPOW and attempts to build
        :class:`MDAnalysis.Universe <MDAnalysis.core.groups.universe.Universe>`
        in the lambda directories.

        Run if :code:`dirname` argument is given when initializing the class.
        First enters FEP directory, then traverses solvent and interaction
        directories to search lambda directories for system files."""
        fep_dir = os.path.join(self._ensemble_dir, 'FEP')
        solv_top_path = None
        for solvent in self._solvents:  # Ugly set of loops, may have to find way to clean up
            if self.top_dict is not None:
                solv_top_path = self.top_dict[solvent]
            for dirs in self._interactions:  # Attribute folder names
                int_dir = os.path.join(fep_dir, solvent, dirs)
                with in_dir(int_dir, create=False):  # Entering attribute folders
                    logger.info("Searching %s directory for systems", os.curdir)
                    files = os.listdir(os.curdir)
                    for file in sorted(files):  # Traversing lambda directories
                        if os.path.isdir(file):
                            with in_dir(file, create=False):
                                u = self._load_universe_from_dir(solv_dir=solv_top_path,
                                                                 **self.unv_kwargs)
                                if u is None:
                                    logger.warning(f'No system loaded in {file}')
                                else:
                                    self.add_system((solvent, dirs, file), u)

    def add_system(self, key, universe: mda.Universe):
        """Adds system from universe object for trajectory and topology files

        Existing mda.Universe object or trajectory and topology path. Ensure
        that paths are set to absolute when creating the universe."""
        self._ensemble[key] = universe
        self._keys.append(key)
        self._num_systems += 1

    def pop(self, key):
        """Removes and returns system at specified key.

        Logs if KeyError is raised."""
        system = self._ensemble.pop(key)
        self._num_systems -= 1
        return system

    def select_atoms(self, *args, **kwargs):
        """Returns :class:`~mdpow.analysis.ensemble.EnsembleAtomGroup` containing selections
        from the :class:`~mdpow.analysis.ensemble.Ensemble`

        Uses the same
        `selection commands <https://docs.mdanalysis.org/stable/documentation_pages/selections.html>`_
        as MDAnalysis, and has the same keys as the :class:`~mdpow.analysis.ensemble.Ensemble`"""
        selections = {}
        for key in self.keys():
            try:
                ag = self[key].select_atoms(*args, **kwargs)
            except SelectionError as err:
                logger.error("%r on system %r with selection settings %r %r", err, key, args, kwargs)
                raise
            else:
                selections[key] = ag
        return EnsembleAtomGroup(selections, ensemble=self)

    def select_systems(self, keys=None, solvents=None, interactions=None,
                       lambdas=None, lambda_range=None):
        """
        Select specific subset of systems and returns them in an Ensemble.

        This can be accomplished in two ways, by specific keys, or by
        specifying the desired system attributes solvents, interactions and
        lambdas. All arguments are stored in list form.

        :keywords:

        *keys*
            System keys from :class:`~mdpow.analysis.ensemble.Ensemble`
            to be returned.

        *solvents*
            Solvents from :class:`~mdpow.analysis.ensemble.Ensemble`
            to be returned.

        *interactions*
            Interactions from :class:`~mdpow.analysis.ensemble.Ensemble`
            to be returned

        *lambdas*
            Specific lambdas to be returned

        *lambda_range*
            Range of lambda to be returned

        .. rubric:: Examples

        Specific key workflow example::

            Ens = Ensemble(dirname='Mol')
            w_v_0_25 = Ens.select_systems(keys=[('water', 'VDW', '0000'),
                                                ('water', 'VDW', '0025')]

        For the system attributes workflow there are two ways of selecting lambdas,
        the :param lambdas: keyword saves specific lambdas or the :param lambda_range:
        which saves the lambdas that fall within the given range.

        Specific lambdas example::

            Ens = Ensemble(dirname='Mol')
            w_v_0_25 = Ens.select_systems(solvents=['water'], interactions=['VDW'],
                                          lambdas=['0000', '0025'])
        Range of lambdas example::

            Ens = Ensemble(dirname='Mol')
            w_v = Ens.select_systems(solvents=['water'], interactions=['VDW'],
                                     lambda_range=[0, 1])

        """
        new_ens = Ensemble()
        new_key = []
        if keys is not None:
            # Selection by giving keys
            new_key = keys
        elif solvents is not None:
            # Selection by attributes
            for s in solvents:
                if interactions is not None:
                    for i in interactions:
                        if lambdas is not None:
                            # Selecting specific lambdas
                            for l in lambdas:
                                new_key.append((s, i, l))
                        elif lambda_range is not None:
                            # Selecting range of lambdas
                            for k in self.keys():
                                if lambda_range[0] <= int(k[2]) / 1000 <= lambda_range[1]:
                                    new_key.append((s, i, k[2]))
        for k in new_key:
            logger.info('adding system %r to ensemble', k)
            new_ens.add_system(k, universe=self[k])
        new_ens._ensemble_dir = self._ensemble_dir
        return new_ens


class EnsembleAtomGroup(object):
    """Group for storing selections from :class:`~mdpow.analysis.ensemble.Ensemble`
     objects made using the :meth:`~mdpow.analysis.ensemble.Ensemble.select_atoms` method.

    :class:`~mdpow.analysis.ensemble.EnsembleAtomGroup` is not set up for manual initialization,
    they should be obtained by selecting atoms from an existing object.
    """

    def __init__(self, group_dict: dict, ensemble: Ensemble):
        self._groups = group_dict
        self._ensemble = ensemble
        self._keys = group_dict.keys()

    def __getitem__(self, index):
        return self._groups[index]

    def __eq__(self, other):
        if self.keys() == other.keys():
            return all(self[k] == other[k] for k in self.keys())
        return False

    def __len__(self):
        return len(self.keys())

    def keys(self):
        """List of keys to specific atom groups in the system"""
        return self._keys

    def positions(self, keys=None):
        """Returns the positions of the keys of the selected atoms.

        If no keys are specified positions for all keys are returned"""
        positions = {}
        if not keys is None:
            for k in keys:
                positions[k] = self._groups[k].positions
        else:
            for k in self.keys():
                positions[k] = self[k].positions
        return positions

    def select_atoms(self, *args, **kwargs):
        """Returns :class:`~mdpow.analysis.ensemble.EnsembleAtomGroup` containing selections
        from the :class:`~mdpow.analysis.ensemble.EnsembleAtomGroup`

        Uses the same
        `selection commands <https://docs.mdanalysis.org/stable/documentation_pages/selections.html>`_
        as MDAnalysis, and has the same keys as :class:`~mdpow.analysis.ensemble.EnsembleAtomGroup`"""
        selections = {}
        for key in self.keys():
            try:
                ag = self[key].select_atoms(*args, **kwargs)
            except SelectionError as err:
                logger.error("%r on system %r with selection settings %r %r", err, key, args, kwargs)
                raise
            else:
                selections[key] = ag
        return EnsembleAtomGroup(selections, ensemble=self._ensemble)

    @property
    def ensemble(self):
        """Returns the ensemble of the EnsembleAtomGroup"""
        return self._ensemble


class EnsembleAnalysis(object):
    """Base class for running multi-system analyses

    The class is designed based on the `AnalysisBase
    class in MDAnalysis <https://docs.mdanalysis.org/stable/documentation_pages/analysis/base.html>`_
    and is a template for creating multi-universe multi-frame analyses using the
    :class:`~mdpow.analysis.ensemble.Ensemble` object

    :Keywords:
    *ensemble*
        The :class:`~mdpow.analysis.ensemble.Ensemble` being analyzed in the class


    .. rubric:: Example Analysis

    Dihedral Analysis Demonstration::

        class DihedralAnalysis(mdpow.ensemble.EnsembleAnalysis):
            def __init__(self, DihedralEnsembleGroup):
                super(DihedralAnalysis, self).__init__(DihedralEnsembleGroup.ensemble)

                self._sel = DihedralEnsembleGroup

            def _prepare_ensemble(self):
                self.result_dict = {}
                for s in ['water', 'octanol']:
                    self.result_dict[s] = {'Coulomb': {},
                                           'VDW': {}}
                for key in self._sel.group_keys():
                    self.result_dict[key[0]][key[1]][key[2]] = None

            def _prepare_universe(self):
                self.angle_dict = {'angle': None,
                                   'time': None}
                self.angles = []

            def _single_frame(self):
                angle = calc_dihedrals(self._sel[self._key].positions[0], self._sel[self._key].positions[1],
                                       self._sel[self._key].positions[2], self._sel[self._key].positions[3])
                self.angles.append(angle)

            def _conclude_universe(self):
                self.angle_dict['time'] = self.times
                self.angle_dict['angle'] = self.angles
                self.result_dict[self._key[0]][self._key[1]][self._key[2]] = self.angle_dict

            def _conclude_ensemble(self):
                self.results = pd.DataFrame(data=self.result_dict)

        D = DihedralAnalysis.run(start=0 stop=10, step=1)

    """

    def __init__(self, ensemble=None):
        self._ensemble = ensemble

    def _setup_system(self, key, start=None, stop=None, step=None):
        self._system = self._ensemble[key]
        self._key = key
        self.start = start
        self.stop = stop
        self.step = step
        self._setup_frames(self._system.trajectory)

    def _setup_frames(self, trajectory):
        self._trajectory = trajectory
        start, stop, step = trajectory.check_slice_indices(self.start, self.stop, self.step)
        self.n_frames = len(range(start, stop, step))
        self.frames = np.zeros(self.n_frames, dtype=int)
        self.times = np.zeros(self.n_frames)

    def _single_universe(self):
        """Calculations on a single Universe object.

            Run on each universe in the ensemble during when
            self.run in called.
        """
        pass  # pragma: no cover

    def _single_frame(self):
        """Calculate data from a single frame of trajectory

        Called on each frame for universes in the Ensemble.
        """
        pass  # pragma: no cover

    def _prepare_ensemble(self):
        """For establishing data structures used in running
        analysis on the entire ensemble.

        Data structures will not be overwritten upon moving to
        next system in ensemble.
        """
        pass  # pragma: no cover

    def _prepare_universe(self):
        """For establishing data structures used in running
        analysis on each trajectory in ensemble

        Data structures will be overwritten between upon after
        each trajectory has been run
        """
        pass  # pragma: no cover

    def _conclude_universe(self):
        """Run after each trajectory is finished"""
        pass  # pragma: no cover

    def _conclude_ensemble(self):
        """Run after all trajectories in ensemble are finished"""
        pass  # pragma: no cover

    def run(self, start=None, stop=None, step=None):
        """Runs _single_universe on each system and _single_frame
        on each frame in the system.

        First iterates through keys of ensemble, then runs _setup_system
        which defines the system and trajectory. Then iterates over
        trajectory frames.
        """
        logger.info("Setting up systems")
        self._prepare_ensemble()
        for self._key in ProgressBar(self._ensemble.keys(), verbose=True):
            self._setup_system(self._key, start=start, stop=stop, step=step)
            self._prepare_universe()
            self._single_universe()
            for i, ts in enumerate(ProgressBar(self._trajectory[self.start:self.stop:self.step], verbose=True,
                                               postfix=f'running system {self._key}')):
                self._frame_index = i
                self._ts = ts
                self.frames[i] = ts.frame
                self.times[i] = ts.time
                self._single_frame()
            self._conclude_universe()
            logger.info("Moving to next universe")
        logger.info("Finishing up")
        self._conclude_ensemble()
        return self

    @staticmethod
    def check_groups_from_common_ensemble(groups: List[EnsembleAtomGroup]):
        """Checks if inputted list of
        :class:`~mdpow.analysis.ensemble.EnsembleAtomGroup` originate from
        the same :class:`~mdpow.analysis.ensemble.Ensemble`

        Checks every :class:`~mdpow.analysis.ensemble.EnsembleAtomGroup` in
        list to determine if their :meth:`ensemble` references the same object
        in memory. If two :class:`~mdpow.analysis.ensemble.EnsembleAtomGroup`
        object don't have a common :class:`~mdpow.analysis.ensemble.Ensemble`
        :class:`ValueError` is raised."""
        for i in range((len(groups) - 1) // 2):
            # Checking if EnsembleAtomGroup.ensemble references same object in memory
            if groups[i].ensemble is not groups[-1 - i]:
                msg = '''Dihedral selections from different Ensembles,
                          ensure that all EnsembleAtomGroups are created
                          from the same Ensemble.'''
                logger.error(msg)
                raise ValueError(msg)
