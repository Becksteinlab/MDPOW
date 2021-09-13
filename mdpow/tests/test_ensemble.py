from __future__ import absolute_import

from . import tempdir as td

import os.path

import py.path

import pybol
import pytest

import numpy as np

import MDAnalysis as mda
from MDAnalysis.exceptions import NoDataError, SelectionError

from gromacs.utilities import in_dir

from ..analysis.ensemble import Ensemble, EnsembleAnalysis, EnsembleAtomGroup
from ..analysis.dihedral import DihedralAnalysis

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
        self.m = pybol.Manifest(str(RESOURCES / 'manifest.yml'))
        self.m.assemble('example_FEP', self.tmpdir.name)

    def teardown(self):
        self.tmpdir.dissolve()

    def test_build_ensemble(self):
        # Octanol will be added later
        Sim = Ensemble(dirname=self.tmpdir.name, solvents=['water'])
        diff = set(Sim.keys()) ^ set(ensemble_keys)
        assert not diff

    def test_kwargs(self):
        l_dir = os.path.abspath(os.path.join(self.tmpdir.name, 'FEP', 'md.gro'))
        bnz = Ensemble(dirname=self.tmpdir.name, solvents=['water'], topology_paths={'water': l_dir})
        diff = set(bnz.keys()) ^ set(ensemble_keys)
        assert not diff

    def test_add_remove_systems(self):
        with in_dir(self.tmpdir.name, create=False):
            bnz = Ensemble()
            l_dir = os.path.join(os.curdir, 'FEP', 'water', 'Coulomb', '0000')
            top_dir = os.path.join(l_dir, 'md.gro')
            trj_dir = os.path.join(l_dir, 'md_red.xtc')
            U = mda.Universe(top_dir, trj_dir)
            bnz.add_system(('water', 'Coulomb', '0000'), U)
            assert bnz.keys() == [('water', 'Coulomb', '0000')]
            assert bnz._num_systems == 1
            assert bnz.__repr__() == "<Ensemble Containing 1 System>"
            assert len(bnz) == 1
            bnz.pop(('water', 'Coulomb', '0000'))
            assert bnz._num_systems == 0
            assert len(bnz) == 0

    def test_select_atoms(self):
        Sim = Ensemble(dirname=self.tmpdir.name, solvents=['water'])
        solute = Sim.select_atoms('not resname SOL')
        assert len(solute) == 7
        for k in solute.keys():
            assert len(solute[k]) == 42

    def test_select_systems(self):
        Sim = Ensemble(dirname=self.tmpdir.name, solvents=['water'])
        Sel1 = Sim.select_systems(keys=[('water', 'Coulomb', '0000'),
                                        ('water', 'VDW', '0500')])
        assert Sel1.keys() == [('water', 'Coulomb', '0000'),
                               ('water', 'VDW', '0500')]
        Sel2 = Sim.select_systems(solvents=['water'], interactions=['Coulomb'],
                                  lambdas=['0000', '1000'])
        assert Sel2.keys() == [('water', 'Coulomb', '0000'),
                               ('water', 'Coulomb', '1000')]
        Sel3 = Sim.select_systems(solvents=['water'], interactions=['VDW'],
                                  lambda_range=[0, 1])
        diff = set(Sel3.keys()) ^ set(ensemble_keys[3:])
        assert not diff

    def test_ensemble_ag_methods(self):
        Solv_system = Ensemble(dirname=self.tmpdir.name, solvents=['water'])
        Sol1 = Solv_system.select_atoms('resname SOL')
        Sol2 = Sol1.select_atoms('resid 2')
        Sol2_pos = Sol2.positions()
        assert len(Sol2_pos) > 0
        for k in Sol2_pos:
            assert np.shape(Sol2_pos[k]) == (3, 3)
        assert not Sol1 == Sol2
        assert isinstance(Sol2, EnsembleAtomGroup)
        assert Sol2 == Sol1.select_atoms('resid 2')
        assert ensemble_keys.sort() == Sol1.ensemble.keys().sort()
        Sol1._groups.pop(('water', 'Coulomb', '0000'))
        Sol1._keys = Sol1._groups.keys()
        assert not Sol1 == Sol2
        pos2 = Sol2.positions(keys=[('water', 'Coulomb', '0000')])
        assert np.shape(pos2[('water', 'Coulomb', '0000')]) == (3, 3)

    def test_ensemble_init_exception(self):
        with pytest.raises(FileNotFoundError):
            Ens = Ensemble(dirname='foo')

    def test_ensemble_build_exceptions(self):
        with pytest.raises(NoDataError):
            ens = Ensemble(self.tmpdir.name, solvents=['test_solv'])

    def test_ensemble_selection_error(self):
        ens = Ensemble(dirname=self.tmpdir.name, solvents=['water'])
        sel1 = ens.select_atoms('resid 1')
        with pytest.raises(SelectionError):
            ens.select_atoms('foo')
        with pytest.raises(SelectionError):
            sel1.select_atoms('foo')

    def test_ensemble_analysis(self):
        class TestAnalysis(EnsembleAnalysis):
            def __init__(self, test_ensemble):
                super(TestAnalysis, self).__init__(test_ensemble)

                self._ens = test_ensemble

            def _prepare_ensemble(self):
                self.key_list = []

            def _single_universe(self):
                self.key_list.append(self._key)

            def _single_frame(self):
                assert len(self._system.select_atoms('not resname SOL')) == 42

            def _conclude_universe(self):
                assert self.n_frames == self.stop

        Sim = Ensemble(dirname=self.tmpdir.name, solvents=['water'])
        TestRun = TestAnalysis(Sim).run(start=0, step=1, stop=10)
        assert Sim.keys() == TestRun.key_list

    def test_value_error(self):
        ens = Ensemble(dirname=self.tmpdir.name, solvents=['water'])
        copy_ens = Ensemble()
        copy_ens._ensemble_dir = self.tmpdir.name
        for k in ens.keys():
            copy_ens.add_system(k, ens[k])
        dh1 = ens.select_atoms('name C4 or name C17 or name S2 or name N3')
        dh2 = copy_ens.select_atoms('name C4 or name C17 or name S2 or name N3')
        dh3 = ens.select_atoms('name C4 or name C17 or name S2 or name N3')
        dh4 = ens.select_atoms('name C4 or name C17 or name S2 or name N3')
        with pytest.raises(ValueError):
            dh_run = DihedralAnalysis([dh1, dh2, dh4, dh3]).run(start=0, stop=4, step=1)
