# MDPOW: _ensemble.py
# 2021 Alia Lescoulie

import os

import numpy as np

import MDAnalysis as mda
from MDAnalysis.lib.log import ProgressBar
from MDAnalysis.exceptions import FileFormatWarning, NoDataError, MissingDataWarning, SelectionError

from gromacs.utilities import in_dir

from . import NoDataWarning

import logging

logger = logging.getLogger('mdpow._ensemble')


class Ensemble(object):
    """Collection of related MDAnalysis Universe objects

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

    .. rubric:: Examples


    Typical workflow for MDPOW directory::

        ens = Ensemble(dirname='molecule')

    Typical workflow for adding universes individually::

        ens = Ensemble()
        ens.add_system(topology='md.gro', trajectory='md.xtc')

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
                 topology_paths=None, interactions=('Coulomb', 'VDW')):
        self.top_dict = topology_paths
        self._num_systems = 0
        self._ensemble = {}
        self._keys = []
        self._solvents = solvents
        self._interactions = interactions

        if dirname is None:
            return

        if not os.path.exists(dirname):
            logger.error(f"Directory {dirname} does not exist")
            raise FileNotFoundError

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
    def _sort_trajectories(trajectories: list):
        """sort list of trajectory files alphabetically"""
        return sorted(trajectories)

    @staticmethod
    def _sort_topologies(topologies: list):
        """sorts list of trajectory files with .tpr first"""
        top = []
        logger.info('If more than one topology is present the tpr will be the one used')
        for i in range(len(topologies)):
            file = topologies[i]
            if file.endswith('.tpr'):
                topologies.pop(i)
                top = [file] + topologies
                break
        return top

    def _load_universe_from_dir(self, solvent):
        """Loads system in directory into an MDAnalysis Universe

        logs warning if more than one topology is in directory. If
        more than one trajectory attempts to load both of them
        in a universe if that fail will try to load each individually"""
        cur_dir = os.listdir(os.curdir)
        self.trj = []
        self.top = []

        if not self.top_dict is None:
            # if top is specified in kwargs, saved to list
            self.top = [self.top_dict[solvent]]

        for file in cur_dir:
            if file.endswith('.xtc'):
                # Saving trajectory directories
                self.trj.append(file)
            elif (file.endswith('gro') or file.endswith('.tpr') or file.endswith('gro.b2z') or file.endswith('gro.gz'))\
                    and self.top_dict is None:
                # Saving topology directories
                self.top.append(file)

        if len(self.top) == 0 or len(self.trj) == 0:
            logger.warning('No MD files detected in %s', os.curdir)
            raise NoDataWarning

        self.trj = self._sort_trajectories(self.trj)
        self.top = self._sort_topologies(self.top)

    def keys(self):
        """Returns list of system keys"""
        return self._keys

    def _build_ensemble(self):
        """Goes through given mdpow directory and attempts to build universe in the lambda directories.
        Run if dirname is passed into __init__. First enters FEP directory, then goes through solvent
        and interaction directories to search lambda directories for system files."""
        fep_dir = os.path.join(self._ensemble_dir, 'FEP')
        for solvent in self._solvents:  # Ugly set of loops, may have to find way to clean up
            for dirs in self._interactions:  # Attribute folder names
                int_dir = os.path.join(fep_dir, solvent, dirs)
                if os.path.exists(int_dir):
                    with in_dir(int_dir, create=False):  # Entering attribute folders
                        logger.info("Searching %s directory for systems", os.curdir)
                        files = os.listdir(os.curdir)
                        for file in files:  # Traversing lambda directories
                            if os.path.isdir(file):
                                with in_dir(file, create=False):
                                    try:
                                        self._load_universe_from_dir(solvent)
                                        self.add_system((solvent, dirs, file), topology=self.top[0],
                                                        trajectory=self.trj)
                                    except NoDataWarning:
                                        logger.warning('Failed to load universe in %s', file)
                                        continue
                else:
                    logger.error('%r directory does not exist', int_dir)
                    raise NoDataError

    def add_system(self, key, topology=None, trajectory=None, universe=None):
        """Adds system from universe object for trajectory and topology files

        Takes specified key and either existing mda.Universe object or
        trajectory and topology path."""
        if not universe is None:
            self._ensemble[key] = universe
            self._keys.append(key)
            self._num_systems += 1
            return
        if isinstance(trajectory, list):
            # setting abs paths for trajectory
            trajectory = [os.path.abspath(trj) for trj in trajectory]
        try:
            self._ensemble[key] = mda.Universe(os.path.abspath(topology), trajectory)
        except (ValueError, FileFormatWarning, NoDataError, MissingDataWarning, OSError) as err:
            logger.warning('Error loading systems at %r', key)
            raise NoDataWarning
        else:
            self._keys.append(key)
            self._num_systems += 1
            return

    def pop(self, key):
        """Removes and returns system at specified key.

        Logs if KeyError is raised."""
        system = self._ensemble.pop(key)
        self._num_systems -= 1
        return system

    def select_atoms(self, *args, **kwargs):
        """Returns :class:`~mdpow.analysis.ensemble.EnsembleAtomGroup` for
         Universes from the :class:`~mdpow.analysis.ensemble.Ensemble`

        Uses the same
        `selection commands <https://docs.mdanalysis.org/stable/documentation_pages/selections.html>`
        as MDAnalysis."""
        selections = {}
        for key in self.keys():
            try:
                ag = self[key].select_atoms(*args, **kwargs)
            except SelectionError as err:
                logger.error("%r on system %r with selection settings %r %r", err, key, args, kwargs)
                raise
            else:
                selections[key] = ag
        return EnsembleAtomGroup(selections, ensemble_dir=self._ensemble_dir)

    def select_systems(self, keys=None, solvents=None, interactions=None,
                       lambdas=None, lambda_range=None):
        """Select specific subset of systems and returns them in an Ensemble.

        This can be accomplished in two ways, by specific keys, or by
        specifying the desired system attributes solvents, interactions and
        lambdas. All arguments are stored in list form.

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
        if not keys is None:
            # Selection by giving keys
            new_key = keys
        elif not solvents is None:
            # Selection by attributes
            for s in solvents:
                if not interactions is None:
                    for i in interactions:
                        if not lambdas is None:
                            # Selecting specific lambdas
                            for l in lambdas:
                                new_key.append((s, i, l))
                        elif not lambda_range is None:
                            # Selecting range of lambdas
                            for k in self.keys():
                                if lambda_range[0] <= int(k[2])/1000 <= lambda_range[1]:
                                    new_key.append((s, i, k[2]))
        for k in new_key:
            logger.info('adding system %r to ensemble', k)
            new_ens.add_system(k, universe=self[k])
        new_ens._ensemble_dir = self._ensemble_dir
        return new_ens


class EnsembleAtomGroup(object):
    """Group for storing selections from :class:`~mdpow.analysis.ensemble.Ensemble`
     objects made using the :meth:`~mdpow.analysis.ensemble.Ensemble.select_atoms` method.
    """

    def __init__(self, group_dict: dict, ensemble_dir=None):
        self._groups = group_dict
        self._ens_dir = ensemble_dir
        self._keys = group_dict.keys()

    def __getitem__(self, index):
        return self._groups[index]

    def __eq__(self, other):
        if self.keys() == other.keys():
            for k in self.keys():
                if self[k] != other[k]:
                    return False
            return True
        else:
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
        """Returns :class:`~mdpow.analysis.ensemble.EnsembleAtomGroup` for
         Universes from the :class:`~mdpow.analysis.ensemble.EnsembleAtomGroup`

        Uses the same
        `selection commands <https://docs.mdanalysis.org/stable/documentation_pages/selections.html>`
        as MDAnalysis."""
        selections = {}
        for key in self.keys():
            try:
                ag = self[key].select_atoms(*args, **kwargs)
            except SelectionError as err:
                logger.error("%r on system %r with selection settings %r %r", err, key, args, kwargs)
                raise
            else:
                selections[key] = ag
        return EnsembleAtomGroup(selections, ensemble_dir=self._ens_dir)

    def ensemble(self):
        """Returns the _ensemble of the EnsembleAtomGroup"""
        ens = Ensemble()
        for k in self.keys():
            ens.add_system(k, universe=self[k].universe)
        ens._ensemble_dir = self._ens_dir
        return ens


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

        class DihedralAnalysis(mdpow._ensemble.EnsembleAnalysis):
            def __init__(self, DihedralEnsembleGroup):
                super(DihedralAnalysis, self).__init__(DihedralEnsembleGroup._ensemble())

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

            Run on each universe in the _ensemble during when
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
        analysis on the entire _ensemble.

        Data structures will not be overwritten upon moving to
        next system in _ensemble.
        """
        pass  # pragma: no cover

    def _prepare_universe(self):
        """For establishing data structures used in running
        analysis on each trajectory in _ensemble

        Data structures will be overwritten between upon after
        each trajectory has been run
        """
        pass  # pragma: no cover

    def _conclude_universe(self):
        """Run after each trajectory is finished"""
        pass  # pragma: no cover

    def _conclude_ensemble(self):
        """Run after all trajectories in _ensemble are finished"""
        pass  # pragma: no cover

    def run(self, start=None, stop=None, step=None):
        """Runs _single_universe on each system and _single_frame
        on each frame in the system.

        First iterates through keys of _ensemble, then runs _setup_system
        which defines the system and trajectory. Then iterates over
        trajectory frames.
        """
        logger.info("Setting up systems")
        self._prepare_ensemble()
        for self._key in ProgressBar(self._ensemble.keys(), verbose=True):
            self._setup_system(self._key, start=start, stop=stop, step=step)
            self._prepare_universe()
            self._single_universe()
            for i, ts in enumerate(ProgressBar(self._trajectory[self.start:self.stop:self.step], verbose=True, postfix=f'running system {self._key}')):
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
