# mdpow: analysis.py
# 2021 Alia Lescoulie

"""

"""
from __future__ import absolute_import

import os

import gromacs
from gromacs.utilities import in_dir
import MDAnalysis as mda
from MDAnalysis.exceptions import FileFormatWarning, NoDataError, MissingDataWarning, SelectionError
import numpy
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

    def add_system(self, key, topology=None, trajectory=None):
        if not os.path.exists(trajectory) and os.path.exists(topology):
            logger.error("%s does not exist" % trajectory)
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
        return selections
