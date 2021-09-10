from __future__ import absolute_import

from . import tempdir as td

import os.path

import py.path

import pybol
import pytest

import numpy as np
import pandas as pd

import MDAnalysis as mda

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