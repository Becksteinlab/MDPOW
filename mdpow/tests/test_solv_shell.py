from __future__ import absolute_import

import numpy as np

from . import tempdir as td

import py.path

import pybol
import pytest

from numpy.testing import assert_almost_equal
from scipy.stats import variation

from ..analysis.ensemble import Ensemble, EnsembleAnalysis, EnsembleAtomGroup

from ..analysis.solvation import SolvationAnalysis

from pkg_resources import resource_filename

RESOURCES = py.path.local(resource_filename(__name__, 'testing_resources'))
MANIFEST = RESOURCES.join("manifest.yml")


class TestSolvShell(object):
    mean = 2654.0
    std = 2654.465059103246

    def setup(self):
        self.tmpdir = td.TempDir()
        self.m = pybol.Manifest(str(RESOURCES / 'manifest.yml'))
        self.m.assemble('example_FEP', self.tmpdir.name)
        self.ens = Ensemble(dirname=self.tmpdir.name, solvents=['water'])
        self.solute = self.ens.select_atoms('not resname SOL')
        self.solvent = self.ens.select_atoms('resname SOL and name OW')

    def teardown(self):
        self.tmpdir.dissolve()

    def test_dataframe(self):
        solv = SolvationAnalysis(self.solute, self.solvent, [1.2]).run(start=0, stop=4, step=1)

        for d in solv.results['distance']:
            assert d == 1.2
        for s in solv.results['solvent']:
            assert s == 'water'
        for i in solv.results['interaction'][:12]:
            assert i == 'Coulomb'

    def test_selection(self):
        solv = SolvationAnalysis(self.solute, self.solvent, [2, 10]).run(start=0, stop=4, step=1)
        mean = np.mean(solv.results['N_solvent'])
        std = np.std(solv.results['N_solvent'])
        assert_almost_equal(mean, self.mean, 6)
        assert_almost_equal(std, self.std, 6)
