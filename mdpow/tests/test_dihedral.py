from __future__ import absolute_import

from . import tempdir as td

import py.path

import pybol
import pytest

from numpy.testing import assert_almost_equal
from scipy.stats import circvar, circmean

from MDAnalysis.exceptions import SelectionError

from ..analysis.ensemble import Ensemble, EnsembleAnalysis, EnsembleAtomGroup
from ..analysis.dihedral import DihedralAnalysis

from pkg_resources import resource_filename

RESOURCES = py.path.local(resource_filename(__name__, 'testing_resources'))
MANIFEST = RESOURCES.join("manifest.yml")


class TestDihedral(object):
    DG48910_mean = -172.9849512527183
    DG491011_mean = -3.7300197060478695
    DG48910_var = 1490.6576365537262
    DG491011_var = 128.3805265432388

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

        assert results['selection'][0] == 'S2-N3-C4-C17'
        for s in results['solvent']:
            assert s == 'water'
        for i in results['interaction'][:12]:
            assert i == 'Coulomb'

    def test_selection_error(self):
        dh1 = self.Ens.select_atoms('name C17 or name S2 or name N3')
        with pytest.raises(SelectionError):
            dh_run = DihedralAnalysis([dh1]).run(start=0, stop=4, step=1)

    def test_results_recursive1(self):
        dh1 = self.Ens.select_atoms('name C11 or name C10 or name C9 or name C4')
        dh2 = self.Ens.select_atoms('name C11 or name C10 or name C9 or name C4')

        dh_run1 = DihedralAnalysis([dh1]).run(start=0, stop=4, step=1)
        dh_run2 = DihedralAnalysis([dh2]).run(start=0, stop=4, step=1)
        assert len(dh_run1.results['dihedral']) == len(dh_run2.results['dihedral'])
        for i in range(len(dh_run1.results['dihedral'])):
            assert dh_run1.results['dihedral'][i] == dh_run2.results['dihedral'][i]

    def test_results_recursive2(self):
        dh1 = self.Ens.select_atoms('name C11 or name C10 or name C9 or name C4')
        dh2 = self.Ens.select_atoms('name C8 or name C4 or name C9 or name C10')

        dh_run = DihedralAnalysis([dh1, dh2]).run(start=0, stop=4, step=1)

        dh1_result = dh_run.results.loc[dh_run.results['selection'] == 'C4-C9-C10-C11']['dihedral']
        dh2_result = dh_run.results.loc[dh_run.results['selection'] == 'C4-C8-C9-C10']['dihedral']

        dh1_mean = circmean(dh1_result, high=180, low=-180)
        dh2_mean = circmean(dh2_result, high=180, low=-180)
        dh1_var = circvar(dh1_result, high=180, low=-180)
        dh2_var = circvar(dh2_result, high=180, low=-180)

        assert_almost_equal(self.DG48910_mean, dh1_mean, 6)
        assert_almost_equal(self.DG48910_var, dh1_var, 6)
        assert_almost_equal(self.DG491011_mean, dh2_mean, 6)
        assert_almost_equal(self.DG491011_var, dh2_var, 6)
