# mdpow: analysis.py
# 2021 Alia Lesocoulie

"""

"""
from __future__ import absolute_import

import os
import glob

import gromacs
from gromacs.utilities import in_dir
import MDAnalysis as mda
from MDAnalysis.exceptions import FileFormatWarning, NoDataError, MissingDataWarning
import numpy
import pandas as pd

import six

import logging
logger = logging.getLogger('mdpow.analysis')
if not hasattr(logger, 'warning'):
    logger.warning = logger.warn()


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

        if not os.path.exists(dirname):
            logger.error("Directory %s does not exist" % dirname)
        self.ensemble_dir = dirname
        self.build_ensemble()

    def __repr__(self):
        print("<Ensemble Containing %s Universes>", self.num_systems)

    def __len__(self):
        return  self.num_systems

    def __getitem__(self, index):
        """Allows dictionary like indexing"""
        try:
            item = self.ensemble[index]
            return item
        except KeyError:
            logger.error(KeyError("Invalid system key"))

    def __iter__(self):
        """Returns list of systems with object iterated over"""
        return [system for system in self[self.system_keys]]

    @staticmethod
    def _load_dir_unv():
        """Loads system in directory into an MDAnalysis Universe

        logs warning if more than one topology is in directory. If
        more than one trajectory attempts to load both of them
        in a universe if that fail will try to load each individually"""
        cur_dir = os.listdir(os.curdir)
        top = []
        trj = []
        for file in cur_dir:
            if '.xvg' in file:
                trj.append(file)
            if '.gro' or '.pdb' in file:
                top.append(file)
        if len(top) > 1:
            logger.warning('More than one topology detected in %s will just use first topology' % os.curdir)
        if len(top) == 0 or len(trj) == 0:
            logger.warning('No MD files detected in %s' % os.curdir)
            raise NoDataWarning

        try:
            system = mda.Universe(top[0], trj)
            return system
        except FileFormatWarning or MissingDataWarning or NoDataError:
            logger.warning('Mutiple trajectories deteceted in %s attempting'
                           ' to find correct one for topology' % os.curdir)
            while len(trj) > 0:
                try:
                    system = mda.Universe(top[0], top.pop())
                    return system
                except FileFormatWarning or MissingDataWarning or NoDataError:
                    continue
            logger.warning('No files in %s are compatible skipping directory' % os.curdir)
            raise NoDataWarning

    def get_keys(self):
        return self.system_keys

    def build_ensemble(self):
        """Goes through given mdpow directory and attempts to build univserse in the lambda directories."""
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
                                with in_dir(file, create=True):
                                    try:
                                        solvent_dict[dirs, file] = self._load_dir_unv()
                                        self.system_keys.append((solvent, dirs, file))
                                    except NoDataWarning:
                                        continue
                            else:
                                pass
                else:
                    pass
            self.ensemble[solvent] = solvent_dict
