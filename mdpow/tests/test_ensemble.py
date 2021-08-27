from __future__ import absolute_import

from . import tempdir as td

import os.path

import py.path

import pybol

import numpy as np

from ..analysis.ensemble import Ensemble, EnsembleAnalysis, EnsembleAtomGroup
from gromacs.utilities import in_dir

from pkg_resources import resource_filename
RESOURCES = py.path.local(resource_filename(__name__, 'testing_resources'))
MANIFEST = RESOURCES.join("manifest.yml")

ensemble_keys = [('water', 'Coulomb', '0000'),
                 ('water', 'Coulomb', '0500'),
                 ('water', 'Coulomb', '1000'),
                 ('water', 'VDW', '0000'),
                 ('water', 'VDW', '0250'),
                 ('water', 'VDW', '0500'),
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
            bnz = Ensemble(dirname='example_FEP', solvents=['water'])
            diff = set(bnz.get_keys()) ^ set(ensemble_keys)
            assert not diff

    def test_kwargs(self):
        with in_dir(os.path.join(self.resources, 'states'), create=False):
            l_dir = os.path.abspath(os.path.join('example_FEP', 'FEP', 'md.gro'))
            bnz = Ensemble(dirname='example_FEP', solvents=['water'], top_paths={'water': l_dir})
            diff = set(bnz.get_keys()) ^ set(ensemble_keys)
            assert not diff

    def test_add_remove_systems(self):
        with in_dir(os.path.join(self.resources, 'states'), create=False):
            bnz = Ensemble()
            l_dir = os.path.join(os.curdir, 'example_FEP', 'FEP', 'water', 'Coulomb', '0000')
            top_dir = os.path.join(l_dir, 'md.gro')
            trj_dir = os.path.join(l_dir, 'md_red.xtc')
            bnz.add_system(('water', 'Coulomb', '0000'), topology=top_dir, trajectory=trj_dir)
            assert bnz.get_keys() == [('water', 'Coulomb', '0000')]
            assert bnz.num_systems == 1
            assert bnz.__repr__() == "<Ensemble Containing 1 System>"
            assert len(bnz) == 1
            bnz.pop_system(('water', 'Coulomb', '0000'))
            assert bnz.num_systems == 0
            assert len(bnz) == 0

    def test_selection(self):
        with in_dir(os.path.join(self.resources, 'states'), create=False):
            bnz = Ensemble(dirname='example_FEP', solvents=['water'])
            solute = bnz.select_atoms('not resname SOL')
            for k in solute.group_keys():
                assert len(solute[k]) == 42

    def test_ensemble_ag_methods(self):
        with in_dir(os.path.join(self.resources, 'states'), create=False):
            BW = Ensemble(dirname='example_FEP', solvents=['water'])
        Sol1 = BW.select_atoms('resname SOL')
        Sol2 = Sol1.select_atoms('resid 1')
        Sol2_pos = Sol2.positions()
        for k in Sol2_pos:
            assert np.shape(Sol2_pos[k]) == (3, 3)
        assert not Sol1 == Sol2
        assert type(Sol2) == EnsembleAtomGroup
        assert Sol2 == Sol1.select_atoms('resid 1')
        assert BW.get_keys() == Sol1.ensemble().get_keys()

    def test_ensemble_analysis(self):
        class TestAnalysis(EnsembleAnalysis):
            def __init__(self, Ensemble):
                super(TestAnalysis, self).__init__(Ensemble)

                self._ens = Ensemble

            def _prepare_ensemble(self):
                self.key_list = []

            def _single_universe(self):
                self.key_list.append(self._key)

            def _single_frame(self):
                assert len(self._system.select_atoms('not resname SOL')) == 42

            def _conclude_universe(self):
                assert self.n_frames == self.stop

        with in_dir(os.path.join(self.resources, 'states'), create=False):
            BNZ = Ensemble(dirname='example_FEP', solvents=['water'])
        TestRun = TestAnalysis(BNZ).run(start=0, step=1, stop=10)
        assert BNZ.get_keys() == TestRun.key_list
