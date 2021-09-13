from __future__ import absolute_import

from . import tempdir as td

import py.path

import pybol
import pytest

from numpy.testing import assert_almost_equal

from ..analysis.ensemble import Ensemble, EnsembleAnalysis, EnsembleAtomGroup

from ..analysis.solvation import SolvationAnalysis

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
        solv = SolvationAnalysis(self.Ens, [1.2]).run(start=0, stop=4, step=1)

        for d in solv.results['distance']:
            assert d == 1.2
        for s in solv.results['solvent']:
            assert s == 'water'
        for i in solv.results['interaction'][:12]:
            assert i == 'Coulomb'
