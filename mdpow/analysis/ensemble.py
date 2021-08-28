# MDPOW: ensemble.py
# 2021 Alia Lescoulie

import os

import numpy as np

import MDAnalysis as mda
from MDAnalysis.lib.log import ProgressBar
from MDAnalysis.exceptions import FileFormatWarning, NoDataError, MissingDataWarning, SelectionError

from gromacs.utilities import in_dir

from . import NoDataWarning

import logging

logger = logging.getLogger('mdpow.ensemble')


class Ensemble(object):
    """Collection of related MDAnalysis Universe objects

    Stores systems produced by running mdpow-fep organized
    by solvent, interaction, and lambda.

    Given a mdpow simulation directory will load the MD
    simulation files with the directory structure as keys.

    Typical workflow for MDPOW directory::

        ens = Ensemble(dirname='molecule')

    Typical workflow for adding universes individually::

        ens = Ensemble()
        ens.add_system(topology='md.gro', trajectory='md.xtc')

    Topology paths can be specified when defining the ensemble
    by giving the paths to each solvent topology in a dictionary
    with the topology_paths argument::

        ens = Ensemble(dirname='molecule', topology_paths={'water': water_path,
                                                           'octanol': octanol_path}

    """

    def __init__(self, dirname=None, solvents=['octanol', 'water'], topology_paths=None):

        self.top = topology_paths
        self.num_systems = 0
        self.ensemble = {}
        self.system_keys = []
        self.solvents = solvents
        self.interactions = ["Coulomb", "VDW"]

        if dirname is None:
            return

        if not os.path.exists(dirname):
            logger.error(f"Directory {dirname} does not exist")
        self.ensemble_dir = dirname
        self._build_ensemble()

    def __repr__(self):
        return f"<Ensemble Containing {self.num_systems} System>"

    def __len__(self):
        return self.num_systems

    def __getitem__(self, index):
        """Allows dictionary like indexing"""
        try:
            item = self.ensemble[index]
            return item
        except KeyError:
            logger.error(KeyError("Invalid system key"))

    @staticmethod
    def sort_trajectories(trajectories: list):
        """sort list of trajectory files alphabetically"""
        return sorted(trajectories)

    def _load_universe_from_dir(self, solvent):
        """Loads system in directory into an MDAnalysis Universe

        logs warning if more than one topology is in directory. If
        more than one trajectory attempts to load both of them
        in a universe if that fail will try to load each individually"""
        cur_dir = os.listdir(os.curdir)
        trj = []
        if self.top is None:
            top = []
        else:
            # if top is specified in kwargs, saved to list
            top = [self.top[solvent]]
        for file in cur_dir:
            if file.endswith('.xtc'):
                trj.append(file)
            if file.endswith('.tpr') and self.top is None:
                # Ensures .tpr is used first if available
                top = [file] + top
            if (file.endswith('gro') or file.endswith('gro.b2z') or file.endswith('gro.gz')) and self.top is None:
                top.append(file)
        if len(top) > 1:
            logger.warning('More than one topology detected in %s will use first topology',
                           os.curdir)
        if len(top) == 0 or len(trj) == 0:
            logger.warning('No MD files detected in %s', os.curdir)
            raise NoDataWarning

        trj = self.sort_trajectories(trj)

        try:
            if self.top is None:
                system = mda.Universe(os.path.abspath(top[0]), [os.path.abspath(p) for p in trj])
            else:
                system = mda.Universe(top[0], [os.path.abspath(p) for p in trj])
            return system
        except (FileFormatWarning, MissingDataWarning, NoDataError, ValueError):
            logger.error('Multiple incompatible trajectories detected in %s', os.curdir)
            raise NoDataWarning

    def get_keys(self):
        return self.system_keys

    def _build_ensemble(self):
        """Goes through given mdpow directory and attempts to build universe in the lambda directories.
        Run if dirname is passed into __init__. First enters FEP directory, then goes through solvent
        and interaction directories to search lambda directories for system files."""
        fep_dir = os.path.join(self.ensemble_dir, 'FEP')
        for solvent in self.solvents:  # Ugly set of loops, may have to find way to clean up
            for dirs in self.interactions:  # Attribute folder names
                int_dir = os.path.join(fep_dir, solvent, dirs)
                if os.path.exists(int_dir):
                    with in_dir(int_dir, create=False):  # Entering attribute folders
                        logger.info("Searching %s directory for systems", os.curdir)
                        files = os.listdir(os.curdir)
                        for file in files:  # Traversing lambda directories
                            if os.path.isdir(file):
                                with in_dir(file, create=False):
                                    try:
                                        self.ensemble[(solvent, dirs, file)] = self._load_universe_from_dir(solvent)
                                    except NoDataWarning:
                                        logger.warning('Failed to load universe in %s', file)
                                        continue
                                    else:
                                        self.system_keys.append((solvent, dirs, file))
                                        self.num_systems += 1

    def add_system(self, key, topology=None, trajectory=None, universe=None):
        """Adds system from universe object for trajectory and topology files

        Takes specified key and either existing mda.Universe object or
        trajectory and topology path."""
        if not universe is None:
            self.ensemble[key] = universe
        elif topology is None or trajectory is None:
            err = f"{topology} does not exist"
            logger.error(err)
            raise ValueError(err)
        elif not os.path.exists(trajectory) and not os.path.exists(topology):
            err = f"{trajectory} does not exist"
            logger.error(err)
            raise ValueError(err)
        else:
            self.ensemble[key] = mda.Universe(os.path.abspath(topology), os.path.abspath(trajectory))
        self.system_keys.append(key)
        self.num_systems += 1

    def pop_system(self, key):
        """Removes and returns system at specified key.

        Logs if KeyError is raised."""
        try:
            system = self.ensemble.pop(key)
            self.num_systems -= 1
            return system
        except KeyError as err:
            logger.error("%s with key %r" % (err, key))

    def select_atoms(self, selection):
        """Returns atom groups for Universe from objects

        Presumes that solute in systems is the same"""
        selections = {}
        for key in self.system_keys:
            try:
                ag = self[key].select_atoms(selection)
                selections[key] = ag
            except SelectionError as err:
                logger.warning("%s on system %r" % (err, key))
                continue
        return EnsembleAtomGroup(selections, Ensemble_dir=self.ensemble_dir)


class EnsembleAtomGroup(object):
    """Group for storing AtomGroups from Ensemble.select_atoms()"""

    def __init__(self, group_dict=None, Ensemble_dir=None):
        self.groups = group_dict
        self.ens_dir = Ensemble_dir
        if isinstance(group_dict, dict):
            self.keys = group_dict.keys()
        else:
            self.keys = None

    def __getitem__(self, index):
        try:
            item = self.groups[index]
            return item
        except KeyError:
            logger.error(KeyError("Invalid system key"))

    def __eq__(self, other):
        if self.group_keys() == other.group_keys():
            for k in self.group_keys():
                if self[k] != other[k]:
                    return False
            return True
        else:
            return False

    def __len__(self):
        return len(self.group_keys())

    def group_keys(self):
        """List of keys to specific atom groups in the system"""
        return self.keys

    def positions(self, keys=[None]):
        """Returns the positions of the keys of the selected atoms."""
        positions = {}
        if not keys is None:
            for k in keys:
                try:
                    positions[k] = self.groups[k].positions
                except KeyError as err:
                    logger.warning(f"{err} is an invalid key")
                    continue
        else:
            for k in self.group_keys():
                positions[k] = self[k].positions
        return positions

    def select_atoms(self, selection):
        """Returns atom groups for Universe from objects

        Presumes that solute in systems is the same"""
        selections = {}
        for key in self.group_keys():
            try:
                ag = self[key].select_atoms(selection)
                selections[key] = ag
            except SelectionError as err:
                logger.warning("%s on system %r" % (err, key))
                continue
        return EnsembleAtomGroup(selections)

    def ensemble(self):
        """Returns the ensemble of the EnsembleAtomGroup"""
        ens = Ensemble()
        for k in self.group_keys():
            ens.add_system(k, universe=self[k].universe)
        ens.ensemble_dir = self.ens_dir
        return ens


class EnsembleAnalysis(object):
    """Base class for running multi system analysis

    The class is designed based on the `AnalysisBase
    class in MDAnalysis <https://docs.mdanalysis.org/stable/documentation_pages/analysis/base.html>`_
    and is a template for creating multiuniverse
    multiframe analyses using the Ensemble object

    Typical workflow::

        class DihedralAnalysis(mdpow.ensemble.EnsembleAnalysis):
            def __init__(self, DihedralEnsembleGroup):
                super(DihedralAnalysis, self).__init__(DihedralEnsembleGroup.ensemble())

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

        DihedralAnalysis.run(start=0 stop=10, step=1)

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
        pass

    def _single_frame(self):
        """Calculate data from a single frame of trajectory

        Called on each frame for universes in the Ensemble.
        """
        pass

    def _prepare_ensemble(self):
        """For establishing data structures used in running
        analysis on the entire ensemble.

        Data structures will not be overwritten upon moving to
        next system in ensemble.
        """
        pass

    def _prepare_universe(self):
        """For establishing data structures used in running
        analysis on each trajectory in ensemble

        Data structures will be overwritten between upon after
        each trajectory has been run
        """
        pass

    def _conclude_universe(self):
        """Run after each trajectory is finished"""
        pass

    def _conclude_ensemble(self):
        """Run after all trajectories in ensemble are finished"""
        pass

    def run(self, start=None, stop=None, step=None):
        """Runs _single_universe on each system and _single_frame
        on each frame in the system.

        First iterates through keys of _ensemble, then runs _setup_system
        which defines the system and trajectory. Then iterates over
        trajectory frames.
        """
        logger.info("Setting up systems")
        self._prepare_ensemble()
        for self._key in self._ensemble.get_keys():
                self._setup_system(self._key, start=start, stop=stop, step=step)
                self._prepare_universe()
                self._single_universe()
                for i, ts in enumerate(ProgressBar(self._trajectory[self.start:self.stop:self.step], verbose=True)):
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
