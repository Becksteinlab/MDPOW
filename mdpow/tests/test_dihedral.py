from __future__ import absolute_import

from . import tempdir as td

import os.path

import py.path

import pybol
import pytest

import numpy as np
import pandas as pd

from MDAnalysis.exceptions import SelectionError

from ..analysis.ensemble import Ensemble, EnsembleAnalysis, EnsembleAtomGroup
from ..analysis.dihedral import DihedralAnalysis

from pkg_resources import resource_filename

RESOURCES = py.path.local(resource_filename(__name__, 'testing_resources'))
MANIFEST = RESOURCES.join("manifest.yml")


class TestDihedral(object):
    def setup(self):
        self.tmpdir = td.TempDir()
        self.m = pybol.Manifest(str(RESOURCES / 'manifest.yml'))
        self.m.assemble('example_FEP', self.tmpdir.name)
        self.Ens = Ensemble(dirname=self.tmpdir.name, solvents=['water'])

    def teardown(self):
        self.tmpdir.dissolve()

    def test_dataframe(self):
        dh1 = self.Ens.select_atoms('name C4 or name C17 or name S2 or name N3')
        dh_run = DihedralAnalysis([dh1]).run(start=0, stop=4, step=1)

        results = dh_run.results

        assert results['selection'][0] == 'S2N3C4C17'
        for s in results['solvent']:
            assert s == 'water'
        for i in results['interaction'][:12]:
            assert i == 'Coulomb'

    def test_memory_error(self):
        copy_ens = Ensemble()
        copy_ens._ensemble_dir = self.tmpdir.name
        for k in self.Ens.keys():
            copy_ens.add_system(k, self.Ens[k])
        dh1 = self.Ens.select_atoms('name C4 or name C17 or name S2 or name N3')
        dh2 = copy_ens.select_atoms('name C4 or name C17 or name S2 or name N3')
        dh3 = self.Ens.select_atoms('name C4 or name C17 or name S2 or name N3')
        dh4 = self.Ens.select_atoms('name C4 or name C17 or name S2 or name N3')
        with pytest.raises(MemoryError):
            dh_run = DihedralAnalysis([dh1, dh2, dh4, dh3]).run(start=0, stop=4, step=1)

    def test_selection_error(self):
        dh1 = self.Ens.select_atoms('name C17 or name S2 or name N3')
        with pytest.raises(SelectionError):
            dh_run = DihedralAnalysis([dh1]).run(start=0, stop=4, step=1)
