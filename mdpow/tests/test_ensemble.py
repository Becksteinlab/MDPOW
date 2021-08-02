from __future__ import absolute_import

from . import tempdir as td

import os.path
import sys

import pytest
import py.path

import yaml
import pybol

from six.moves import cPickle as pickle

from mdpow.ensemble import Ensemble
import mdpow.fep
from gromacs.utilities import in_dir

from pkg_resources import resource_filename
RESOURCES = py.path.local(resource_filename(__name__, 'testing_resources'))
MANIFEST = RESOURCES.join("manifest.yml")

ensemble_keys = [('water', 'Coulomb', '0000'),
                 ('water', 'Coulomb', '0250'),
                 ('water', 'Coulomb', '0500'),
                 ('water', 'Coulomb', '0750'),
                 ('water', 'Coulomb', '1000'),
                 ('water', 'VDW', '0000'),
                 ('water', 'VDW', '0050'),
                 ('water', 'VDW', '0100'),
                 ('water', 'VDW', '0200'),
                 ('water', 'VDW', '0300'),
                 ('water', 'VDW', '0400'),
                 ('water', 'VDW', '0500'),
                 ('water', 'VDW', '0600'),
                 ('water', 'VDW', '0650'),
                 ('water', 'VDW', '0700'),
                 ('water', 'VDW', '0750'),
                 ('water', 'VDW', '0800'),
                 ('water', 'VDW', '0850'),
                 ('water', 'VDW', '0900'),
                 ('water', 'VDW', '0950'),
                 ('water', 'VDW', '1000')]


class TestEnsemble(object):
    def setup(self):
        self.tmpdir = td.TempDir()
        self.old_path = os.getcwd()
        self.resources = os.path.join(
            self.old_path, 'mdpow', 'tests', 'testing_resources')
        m = pybol.Manifest(os.path.join(self.resources, 'manifest.yml'))
        m.assemble('FEP', self.tmpdir.name)

    def teardown(self):
        self.tmpdir.dissolve()

    def test_build_ensemble(self):
        with in_dir(os.path.join(self.resources, 'states'), create=False):
            bnz = Ensemble(dirname='benzene', solvents=['water'])
            diff = set(bnz.get_keys()) ^ set(ensemble_keys)
            assert not diff

    def test_add_remove_systems(self):
        with in_dir(os.path.join(self.resources, 'states'), create=False):
            bnz = Ensemble()
            l_dir = os.path.join(os.curdir, 'benzene', 'FEP', 'water', 'VDW', '0000')
            top_dir = os.path.join(l_dir, 'md.gro')
            trj_dir = os.path.join(l_dir, 'md.xtc')
            bnz.add_system(('water', 'VDW', '0000'), topology=top_dir, trajectory=trj_dir)
            assert bnz.get_keys() == [('water', 'VDW', '0000')]
            assert bnz.num_systems == 1
            assert len(bnz) == 1
            bnz.pop_system(('water', 'VDW', '0000'))
            assert bnz.num_systems == 0
            assert len(bnz) == 0

    def test_selection(self):
        with in_dir(os.path.join(self.resources, 'states'), create=False):
            bnz = Ensemble(dirname='benzene', solvents=['water'])
            solute = bnz.select_atoms('not resname SOL')
            for k in solute.keys():
                assert len(solute[k]) == 12
