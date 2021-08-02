# mdpow: ensemble.py
# 2021 Alia Lescoulie

"""A set of objects for analyzing MDPOW simulations

"""
from __future__ import absolute_import

import os

from gromacs.utilities import in_dir
import MDAnalysis as mda
from MDAnalysis.analysis.dihedrals import calc_dihedrals
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.lib.log import ProgressBar
from MDAnalysis.exceptions import FileFormatWarning, NoDataError, MissingDataWarning, SelectionError
import numpy as np
import pandas as pd

import six

import logging

logger = logging.getLogger('mdpow.analysis')


class NoDataWarning(Exception):
    pass


class Ensemble(object):
    """Collection of related MDAnalysis Universe objects

    Stores systems produced by running mdpow-fep organized
    by solvent, interaction, and lambda.

    Given a mdpow simulation directory will load the MD
    simulation files with the directory structure as keys.

    Typical work flow::

        ens = Ensemble(dirname='molecule')

    """

    def __init__(self, dirname=None, solvents=['octanol', 'water'], **kwargs):
        """

        """
        self.num_systems = 0
        self.ensemble = {}
        self.system_keys = []
        self.solvents = solvents

        if dirname is None:
            return

        if not os.path.exists(dirname):
            logger.error("Directory %s does not exist" % dirname)
        self.ensemble_dir = dirname
        self.build_ensemble()

    def __repr__(self):
        return "<Ensemble Containing %s Universes>", self.num_systems

    def __len__(self):
        return self.num_systems

    def __getitem__(self, index):
        """Allows dictionary like indexing"""
        try:
            item = self.ensemble[index]
            return item
        except KeyError:
            logger.error(KeyError("Invalid system key"))

    def __iter__(self):
        """Returns list of systems with object iterated over"""
        sys_list = []
        for key in self.system_keys:
            sys_list.append(self[key])
        return sys_list

    @staticmethod
    def _load_dir_unv():
        """Loads system in directory into an MDAnalysis Universe

        logs warning if more than one topology is in directory. If
        more than one trajectory attempts to load both of them
        in a universe if that fail will try to load each individually"""
        cur_dir = os.listdir(os.curdir)
        trj = []
        top = []
        for file in cur_dir:
            if file.endswith('.xtc'):
                trj.append(file)
            if file.endswith('.gro'):
                top.append(file)
        if len(top) > 1:
            logger.warning('More than one topology detected in %s will just use first topology' % os.curdir)
        if len(top) == 0 or len(trj) == 0:
            logger.warning('No MD files detected in %s' % os.curdir)
            raise NoDataWarning

        try:
            system = mda.Universe(top[0], trj)
            return system
        except FileFormatWarning or MissingDataWarning or NoDataError or ValueError:
            logger.warning('Mutiple trajectories deteceted in %s attempting'
                           ' to find correct one for topology' % os.curdir)
            while len(trj) > 0:
                try:
                    system = mda.Universe(top[0], top.pop())
                    return system
                except FileFormatWarning or MissingDataWarning or NoDataError or ValueError:
                    continue
            logger.warning('No files in %s are compatible skipping directory' % os.curdir)
            raise NoDataWarning

    def get_keys(self):
        return self.system_keys

    def build_ensemble(self):
        """Goes through given mdpow directory and attempts to build universe in the lambda directories."""
        fep_dir = os.path.join(self.ensemble_dir, 'FEP')
        for solvent in self.solvents:  # Ugly set of loops, may have to find way to clean up
            solvent_dict = {'Coulomb': {},
                            'VDW': {}}
            for dirs in ["Coulomb", "VDW"]:
                int_dir = os.path.join(fep_dir, solvent, dirs)

                if os.path.exists(int_dir):
                    with in_dir(int_dir, create=False):
                        logger.info("Searching %s directory for systems" % os.curdir)
                        files = os.listdir(os.curdir)
                        for file in files:
                            if os.path.isdir(file):
                                with in_dir(file, create=False):
                                    try:
                                        self.ensemble[solvent, dirs, file] = self._load_dir_unv()
                                        self.system_keys.append((solvent, dirs, file))
                                        self.num_systems += 1
                                    except NoDataWarning:
                                        continue

    def add_system(self, key, topology=None, trajectory=None, universe=None):
        if not universe is None:
            self.ensemble[key] = universe
        elif topology is None or trajectory is None:
            logger.error("%s does not exist" % trajectory)
        elif not os.path.exists(trajectory) and not os.path.exists(topology):
            logger.error("%s does not exist" % trajectory)
        else:
            self.ensemble[key] = mda.Universe(topology, trajectory)
        self.system_keys.append(key)
        self.num_systems += 1

    def pop_system(self, key):
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
    """Group for storing AtomGroups from Ensemble.select_atoms"""

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
        if len(self) == len(other):
            for k_s, k_o in self.group_keys(), other.group_keys:
                if self[k_s] != other[k_o]:
                    return False
            return True
        else:
            return False

    def __len__(self):
        return len(self.group_keys())

    def group_keys(self):
        return self.keys

    def positions(self, keys=[None]):
        positions = {}
        if not keys is None:
            for k in keys:
                try:
                    positions[k] = self[k].positions
                except KeyError:
                    logger.warning("%s is an invalid key")
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
        ens = Ensemble(dirname=self.ens_dir)
        for k in self.group_keys():
            ens.add_system(k, universe=self[k].universe)
        return ens


class EnsembleAnalysis(object):
    """Base class for running multi system analysis

    The class is designed based on the AnalysisBase
    class in MDAnalysis https://docs.mdanalysis.org
    and is a template for creating multiuniverse
    multiframe analyses using the Ensemble object
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
        """Runs _single_system on each system and _single_frame
        on each frame in the system.
        """
        logger.info("Setting up systems")
        self._prepare_ensemble()
        with in_dir(os.path.join(self._sel.ens_dir, 'FEP'), create=False):
            for self._key in self._ensemble.get_keys():
                with in_dir(os.path.join(os.curdir, self._key[0], self._key[1], self._key[2]),
                            create=False):
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
