import numpy as np
import pandas as pd

from . import tempdir as td

import py.path

import pybol
import pytest

from numpy.testing import assert_almost_equal

from ..analysis.ensemble import Ensemble

from ..analysis.solvation import SolvationAnalysis

from pkg_resources import resource_filename

RESOURCES = py.path.local(resource_filename(__name__, 'testing_resources'))
MANIFEST = RESOURCES.join("manifest.yml")


class TestSolvShell(object):
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
        assert isinstance(solv.results, pd.DataFrame)

        for d in solv.results['distance']:
            assert d == 1.2
        for s in solv.results['solvent']:
            assert s == 'water'
        for i in solv.results['interaction'][:12]:
            assert i == 'Coulomb'

    @pytest.fixture(scope='class')
    def solvation_analysis_list_results(self):
        self.setup()  # Won't have solute and solvent without this
        return SolvationAnalysis(self.solute, self.solvent, [2, 10]).run(start=0, stop=4, step=1)

    @pytest.mark.parametrize("d,ref_mean,ref_std", [(2, 1.10714285,2.07604166),
                                                    (10, 5306.89285714, 129.16720594)])
    def test_selection(self, solvation_analysis_list_results, d, ref_mean, ref_std):
        results = solvation_analysis_list_results.results
        mean = np.mean(results.loc[results['distance'] == d]['N_solvent'])
        std = np.std(results.loc[results['distance'] == d]['N_solvent'])

        assert mean == pytest.approx(ref_mean)
        assert std == pytest.approx(ref_std)
