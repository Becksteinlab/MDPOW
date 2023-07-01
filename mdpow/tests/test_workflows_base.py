import re
import os
import sys
import yaml
import pybol
import pytest
import pathlib
import logging

import pandas as pd

from mdpow.workflows import base

from pkg_resources import resource_filename

from . import RESOURCES, MANIFEST, STATES


@pytest.fixture(scope='function')
def molname_workflows_directory(tmp_path):
    m = pybol.Manifest(str(MANIFEST))
    m.assemble('workflows', tmp_path)
    return tmp_path

class TestWorkflowsBase(object):

    @pytest.fixture(scope='function')
    def SM_tmp_dir(self, molname_workflows_directory):
        dirname = molname_workflows_directory
        return dirname

    @pytest.fixture(scope='function')
    def csv_input_data(self):
        csv_path = STATES['workflows'] / 'project_paths.csv'
        csv_df = pd.read_csv(csv_path).reset_index(drop=True)
        return csv_path, csv_df

    @pytest.fixture(scope='function')
    def test_df_data(self):
        test_dict = {'molecule' : ['SM25', 'SM26'],
                    'resname' : ['SM25', 'SM26']}
        test_df = pd.DataFrame(test_dict).reset_index(drop=True)
        return test_df

    @pytest.fixture(scope='function')
    def project_paths_data(self, SM_tmp_dir):
        project_paths = base.project_paths(parent_directory=SM_tmp_dir)
        return project_paths

    def test_project_paths(self, test_df_data, project_paths_data):
        test_df = test_df_data
        project_paths = project_paths_data

        assert project_paths['molecule'][0] == test_df['molecule'][0]
        assert project_paths['molecule'][1] == test_df['molecule'][1]
        assert project_paths['resname'][0] == test_df['resname'][0]
        assert project_paths['resname'][1] == test_df['resname'][1]

    def test_project_paths_csv_input(self, csv_input_data):
        csv_path, csv_df = csv_input_data
        project_paths = base.project_paths(csv=csv_path)

        pd.testing.assert_frame_equal(project_paths, csv_df)

    def test_dihedral_analysis_figdir_requirement(self, project_paths_data, caplog):
        caplog.clear()
        caplog.set_level(logging.ERROR, logger='mdpow.workflows.base')
        
        project_paths = project_paths_data
        # change resname to match topology (every SAMPL7 resname is 'UNK')
        # only necessary for this dataset, not necessary for normal use
        project_paths['resname'] = 'UNK'

        with pytest.raises(AssertionError,
                           match="figdir MUST be set, even though it is a kwarg. Will be changed with #244"):
        
            base.automated_project_analysis(project_paths, solvents=('water',),
                                            ensemble_analysis='DihedralAnalysis')

            assert 'all analyses completed' in caplog.text, ('automated_dihedral_analysis '
                   'did not iteratively run to completion for the provided project')

    def test_automated_project_analysis_KeyError(self, project_paths_data, caplog):
        caplog.clear()
        caplog.set_level(logging.ERROR, logger='mdpow.workflows.base')

        project_paths = project_paths_data
        # change resname to match topology (every SAMPL7 resname is 'UNK')
        # only necessary for this dataset, not necessary for normal use
        project_paths['resname'] = 'UNK'

        # test error output when raised
        with pytest.raises(KeyError,
                           match="Invalid ensemble_analysis 'DarthVaderAnalysis'. "
                                 "An EnsembleAnalysis type that corresponds to an existing "
                                 "automated workflow module must be input as a kwarg. ex: "
                                 "ensemble_analysis='DihedralAnalysis'"):
            base.automated_project_analysis(project_paths, ensemble_analysis='DarthVaderAnalysis', solvents=('water',))

        # test logger error recording
        assert "'DarthVaderAnalysis' is an invalid selection" in caplog.text, ('did not catch incorrect '
               'key specification for workflows.registry that results in KeyError')
