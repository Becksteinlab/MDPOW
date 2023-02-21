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

import py.path

from ..workflows import base

from pkg_resources import resource_filename

RESOURCES = pathlib.PurePath(resource_filename(__name__, 'testing_resources'))
MANIFEST = RESOURCES / "manifest.yml"

@pytest.fixture(scope="function")
def molname_workflows_directory(tmp_path):
    m = pybol.Manifest(str(MANIFEST))
    m.assemble('workflows', tmp_path)
    return tmp_path

class TestWorkflowsBase(object):

    @pytest.fixture(scope="function")
    def SM_tmp_dir(self, molname_workflows_directory):
        dirname = molname_workflows_directory
        return dirname

    test_dict = {'molecule' : ['SM25', 'SM26'],
                    'resname' : ['SM25', 'SM26']
                   }

    test_df = pd.DataFrame(test_dict).reset_index(drop=True)

    csv_path = 'mdpow/tests/testing_resources/states/workflows/dir_paths.csv'
    csv_df = pd.read_csv(csv_path).reset_index(drop=True)
        
    def test_directory_paths(self, SM_tmp_dir):
        directory_paths = base.directory_paths(parent_directory=SM_tmp_dir)

        assert directory_paths['molecule'][0] == self.test_df['molecule'][0]
        assert directory_paths['molecule'][1] == self.test_df['molecule'][1]
        assert directory_paths['resname'][0] == self.test_df['resname'][0]
        assert directory_paths['resname'][1] == self.test_df['resname'][1]

    def test_directory_paths_csv_input(self):
        directory_paths = base.directory_paths(csv=self.csv_path)
        
        pd.testing.assert_frame_equal(directory_paths, self.csv_df)
        # is additional assertion required here? can this ever fail?

    def test_directory_iteration(self, SM_tmp_dir, caplog):
        directory_paths = base.directory_paths(parent_directory=SM_tmp_dir)
        # change resname to match topology (every SAMPL7 resname is 'UNK')
        # only necessary for this dataset, not necessary for normal use
        directory_paths['resname'] = 'UNK'

        base.directory_iteration(directory_paths, solvents=('water',),
                                 ensemble_analysis='DihedralAnalysis')

        assert 'all analyses completed' in caplog.text, 'automated_dihedral_analysis did not iteratively run to completion for the provided project'
