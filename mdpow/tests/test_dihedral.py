from __future__ import absolute_import

from . import tempdir as td

import sys

import py.path

import pybol
import pytest

from numpy.testing import assert_almost_equal
from scipy.stats import circvar, circmean

from MDAnalysis.exceptions import SelectionError

from ..analysis.ensemble import Ensemble, EnsembleAnalysis, EnsembleAtomGroup
from ..analysis.dihedral import DihedralAnalysis

from pkg_resources import resource_filename

RESOURCES = py.path.local(resource_filename(__name__, "testing_resources"))
MANIFEST = RESOURCES.join("manifest.yml")


class TestDihedral(object):
    DG48910_mean = -172.9849512527183
    DG491011_mean = 177.74725233051953
    DG48910_var = 0.20311120667628546
    DG491011_var = 0.006976126708773456

    def setup_method(self):
        self.tmpdir = td.TempDir()
        self.m = pybol.Manifest(str(RESOURCES / "manifest.yml"))
        self.m.assemble("example_FEP", self.tmpdir.name)
        self.Ens = Ensemble(dirname=self.tmpdir.name, solvents=["water"])

    def teardown_method(self):
        self.tmpdir.dissolve()

    def test_dataframe(self):
        dh1 = self.Ens.select_atoms("name C4", "name C17", "name S2", "name N3")
        dh_run = DihedralAnalysis([dh1]).run(start=0, stop=4, step=1)

        results = dh_run.results

        assert results["selection"][0] == "C4-C17-S2-N3"
        for s in results["solvent"]:
            assert s == "water"
        for i in results["interaction"][:12]:
            assert i == "Coulomb"

    def test_selection_error(self):
        dh1 = self.Ens.select_atoms("name C17", "name S2", "name N3")
        with pytest.raises(SelectionError):
            dh_run = DihedralAnalysis([dh1]).run(start=0, stop=4, step=1)

    def test_results_recursive1(self):
        dh1 = self.Ens.select_atoms("name C11", "name C10", "name C9", "name C4")
        dh2 = self.Ens.select_atoms("name C11", "name C10", "name C9", "name C4")

        dh_run1 = DihedralAnalysis([dh1]).run(start=0, stop=4, step=1)
        dh_run2 = DihedralAnalysis([dh2]).run(start=0, stop=4, step=1)
        assert len(dh_run1.results["dihedral"]) == len(dh_run2.results["dihedral"])
        for i in range(len(dh_run1.results["dihedral"])):
            assert dh_run1.results["dihedral"][i] == dh_run2.results["dihedral"][i]

    @pytest.mark.skipif(
        sys.version_info < (3, 8), reason="scipy circvar gives wrong answers"
    )
    def test_results_recursive2(self):
        dh1 = self.Ens.select_atoms("name C11", "name C10", "name C9", "name C4")
        dh2 = self.Ens.select_atoms("name C8", "name C4", "name C9", "name C10")

        dh_run = DihedralAnalysis([dh1, dh2]).run(start=0, stop=4, step=1)

        dh1_result = dh_run.results.loc[dh_run.results["selection"] == "C11-C10-C9-C4"][
            "dihedral"
        ]
        dh2_result = dh_run.results.loc[dh_run.results["selection"] == "C8-C4-C9-C10"][
            "dihedral"
        ]

        dh1_mean = circmean(dh1_result, high=180, low=-180)
        dh2_mean = circmean(dh2_result, high=180, low=-180)
        dh1_var = circvar(dh1_result, high=180, low=-180)
        dh2_var = circvar(dh2_result, high=180, low=-180)

        assert_almost_equal(dh1_mean, self.DG48910_mean, 6)
        assert_almost_equal(dh1_var, self.DG48910_var, 6)
        assert_almost_equal(dh2_mean, self.DG491011_mean, 6)
        assert_almost_equal(dh2_var, self.DG491011_var, 6)

    def test_ValueError_different_ensemble(self):
        other = Ensemble(dirname=self.tmpdir.name, solvents=["water"])
        dh1 = self.Ens.select_atoms("name C11", "name C10", "name C9", "name C4")
        dh2 = other.select_atoms("name C8", "name C4", "name C9", "name C10")
        with pytest.raises(
            ValueError, match="Dihedral selections from different Ensembles, "
        ):
            DihedralAnalysis([dh1, dh2])

    def test_single_universe(self):
        dh = self.Ens.select_atoms("name C4", "name C17", "name S2", "name N3")
        with pytest.raises(NotImplementedError):
            DihedralAnalysis([dh])._single_universe()
