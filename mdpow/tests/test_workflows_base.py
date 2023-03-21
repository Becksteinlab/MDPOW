import re
import os
import sys
import yaml
import pybol
import pytest
import pathlib
import logging

import pandas as pd

from . import RESOURCES
from . import STATES

import py.path

from ..workflows import base

from pkg_resources import resource_filename

RESOURCES = pathlib.PurePath(resource_filename(__name__, 'testing_resources'))
MANIFEST = RESOURCES / 'manifest.yml'

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
        csv_path = STATES['workflows'] / 'dir_paths.csv'
        csv_df = pd.read_csv(csv_path).reset_index(drop=True)
        return csv_path, csv_df

    @pytest.fixture(scope='function')
    def test_df_data(self):
        test_dict = {'molecule' : ['SM25', 'SM26'],
                    'resname' : ['SM25', 'SM26']}
        test_df = pd.DataFrame(test_dict).reset_index(drop=True)
        return test_df

    @pytest.fixture(scope='function')
    def dir_paths_data(self, SM_tmp_dir):
        directory_paths = base.directory_paths(parent_directory=SM_tmp_dir)
        return directory_paths

    def test_directory_paths(self, test_df_data, dir_paths_data):
        test_df = test_df_data
        directory_paths = dir_paths_data

        assert directory_paths['molecule'][0] == test_df['molecule'][0]
        assert directory_paths['molecule'][1] == test_df['molecule'][1]
        assert directory_paths['resname'][0] == test_df['resname'][0]
        assert directory_paths['resname'][1] == test_df['resname'][1]

    def test_directory_paths_csv_input(self, csv_input_data):
        csv_path, csv_df = csv_input_data
        directory_paths = base.directory_paths(csv=csv_path)

        pd.testing.assert_frame_equal(directory_paths, csv_df)

    def test_directory_iteration(self, dir_paths_data, caplog):
        directory_paths = dir_paths_data
        # change resname to match topology (every SAMPL7 resname is 'UNK')
        # only necessary for this dataset, not necessary for normal use
        directory_paths['resname'] = 'UNK'

        base.directory_iteration(directory_paths, solvents=('water',),
                                 ensemble_analysis='DihedralAnalysis')

        assert 'all analyses completed' in caplog.text, ('automated_dihedral_analysis '
               'did not iteratively run to completion for the provided project')

    def test_directory_iteration_KeyError(self, dir_paths_data, caplog):
        caplog.clear()
        caplog.set_level(logging.ERROR, logger='mdpow.workflows.base')

        directory_paths = dir_paths_data
        # change resname to match topology (every SAMPL7 resname is 'UNK')
        # only necessary for this dataset, not necessary for normal use
        directory_paths['resname'] = 'UNK'

        # test error output when raised
        with pytest.raises(KeyError,
                           match="Invalid ensemble_analysis 'DarthVaderAnalysis'. "
                                 "An EnsembleAnalysis type that corresponds to an existing "
                                 "automated workflow module must be input as a kwarg. ex: "
                                 "ensemble_analysis='DihedralAnalysis'"):
            base.directory_iteration(directory_paths, ensemble_analysis='DarthVaderAnalysis', solvents=('water',))

        # test logger error recording
        assert "'DarthVaderAnalysis' is an invalid selection" in caplog.text, ('did not catch incorrect '
               'key specification for workflows.registry that results in KeyError')

    def test_directory_iteration_TypeError(self, dir_paths_data, caplog):
        # this will currently test for input of analysis type
        # that does not have a corresponding automation module
        # this test and the corresponding error catcher might
        # not be necessary once all automation workflows
        # modules in the registry are completed
        caplog.clear()
        caplog.set_level(logging.ERROR, logger='mdpow.workflows.base')

        directory_paths = dir_paths_data
        # change resname to match topology (every SAMPL7 resname is 'UNK')
        # only necessary for this dataset, not necessary for normal use
        directory_paths['resname'] = 'UNK'

        # test error output when raised
        with pytest.raises(TypeError,
                           match="Invalid ensemble_analysis SolvationAnalysis. An EnsembleAnalysis "
                                 "type that corresponds to an existing automated workflow module must "
                                 "be input as a kwarg. ex: ensemble_analysis='DihedralAnalysis'"):
            base.directory_iteration(directory_paths, ensemble_analysis='SolvationAnalysis', solvents=('water',))

        # test logger error recording
        assert 'workflow module for SolvationAnalysis does not exist yet' in caplog.text, ('did not catch incorrect '
               'key specification for workflows.registry that results in TypeError (NoneType)')
