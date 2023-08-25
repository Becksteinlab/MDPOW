import re
import os
import sys
import yaml
import pathlib
import logging

import pybol
import pytest
import numpy as np
from numpy.testing import assert_equal
import pandas as pd
import MDAnalysis as mda

from . import RESOURCES, MANIFEST, STATES
from pkg_resources import resource_filename
from mdpow.workflows import base


@pytest.fixture
def molname_workflows_directory(tmp_path):
    m = pybol.Manifest(str(MANIFEST))
    m.assemble("workflows", tmp_path)
    return tmp_path


@pytest.fixture
def universe(request):
    masses, names = request.param
    # build minimal test universe
    u = mda.Universe.empty(n_atoms=len(names))
    u.add_TopologyAttr("names", names)
    u.add_TopologyAttr("masses", masses)
    return u


@pytest.mark.parametrize(
    "universe,elements",
    [
        [
            (
                np.array([12.011, 14.007, 0, 12.011, 35.45, 12.011]),
                np.array(["C", "Nx", "DUMMY", "C0S", "Cl123", "C0U"]),
            ),
            np.array(["C", "N", "DUMMY", "C", "CL", "C"]),
        ],
        [
            (
                np.array([12.011, 14.007, 0, 35.45]),
                np.array(["C", "Nx", "DUMMY", "Cl123"]),
            ),
            np.array(["C", "N", "DUMMY", "CL"]),
        ],
        [
            (
                np.array([15.999, 0, 40.08, 40.08, 40.08, 24.305, 132.9]),
                np.array(["OW", "MW", "C0", "CAL", "CA2+", "MG2+", "CES"]),
            ),
            np.array(["O", "DUMMY", "CA", "CA", "CA", "MG", "CS"]),
        ],
        [
            (np.array([16, 1e-6, 40.085, 133]), np.array(["OW", "MW", "CA2+", "CES"])),
            np.array(["O", "DUMMY", "CA", "CS"]),
        ],
    ],
    indirect=["universe"],
)
def test_guess_elements(universe, elements):
    u = universe
    guessed_elements = base.guess_elements(u.atoms)

    assert_equal(guessed_elements, elements)


class TestWorkflowsBase(object):
    @pytest.fixture
    def SM_tmp_dir(self, molname_workflows_directory):
        dirname = molname_workflows_directory
        return dirname

    @pytest.fixture
    def csv_input_data(self):
        csv_path = STATES["workflows"] / "project_paths.csv"
        csv_df = pd.read_csv(csv_path).reset_index(drop=True)
        return csv_path, csv_df

    @pytest.fixture
    def test_df_data(self):
        test_dict = {"molecule": ["SM25", "SM26"], "resname": ["SM25", "SM26"]}
        test_df = pd.DataFrame(test_dict).reset_index(drop=True)
        return test_df

    @pytest.fixture
    def project_paths_data(self, SM_tmp_dir):
        project_paths = base.project_paths(parent_directory=SM_tmp_dir)
        return project_paths

    def test_project_paths(self, test_df_data, project_paths_data):
        test_df = test_df_data
        project_paths = project_paths_data

        assert project_paths["molecule"][0] == test_df["molecule"][0]
        assert project_paths["molecule"][1] == test_df["molecule"][1]
        assert project_paths["resname"][0] == test_df["resname"][0]
        assert project_paths["resname"][1] == test_df["resname"][1]

    def test_project_paths_csv_input(self, csv_input_data):
        csv_path, csv_df = csv_input_data
        project_paths = base.project_paths(csv=csv_path)

        pd.testing.assert_frame_equal(project_paths, csv_df)

    def test_dihedral_analysis_figdir_requirement(self, project_paths_data, caplog):
        caplog.clear()
        caplog.set_level(logging.ERROR, logger="mdpow.workflows.base")

        project_paths = project_paths_data
        # change resname to match topology (every SAMPL7 resname is 'UNK')
        # only necessary for this dataset, not necessary for normal use
        project_paths["resname"] = "UNK"

        with pytest.raises(
            AssertionError,
            match="figdir MUST be set, even though it is a kwarg. Will be changed with #244",
        ):
            base.automated_project_analysis(
                project_paths, solvents=("water",), ensemble_analysis="DihedralAnalysis"
            )

            assert "all analyses completed" in caplog.text, (
                "automated_dihedral_analysis "
                "did not iteratively run to completion for the provided project"
            )

    def test_automated_project_analysis_KeyError(self, project_paths_data, caplog):
        caplog.clear()
        caplog.set_level(logging.ERROR, logger="mdpow.workflows.base")

        project_paths = project_paths_data
        # change resname to match topology (every SAMPL7 resname is 'UNK')
        # only necessary for this dataset, not necessary for normal use
        project_paths["resname"] = "UNK"

        # test error output when raised
        with pytest.raises(
            KeyError,
            match="Invalid ensemble_analysis 'DarthVaderAnalysis'. "
            "An EnsembleAnalysis type that corresponds to an existing "
            "automated workflow module must be input as a kwarg. ex: "
            "ensemble_analysis='DihedralAnalysis'",
        ):
            base.automated_project_analysis(
                project_paths,
                ensemble_analysis="DarthVaderAnalysis",
                solvents=("water",),
            )

        # test logger error recording
        assert "'DarthVaderAnalysis' is an invalid selection" in caplog.text, (
            "did not catch incorrect "
            "key specification for workflows.registry that results in KeyError"
        )
