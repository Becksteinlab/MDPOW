import os
import sys
import yaml
import pybol
import pytest
import pathlib
import logging

import scipy
import numpy as np
import pandas as pd

import rdkit
from rdkit import Chem

import seaborn

from numpy.testing import assert_almost_equal
from scipy.stats import circvar, circmean

from . import RESOURCES

import py.path

#from ..workflows import dihedrals
from ..workflows import base

from pkg_resources import resource_filename

RESOURCES = pathlib.PurePath(resource_filename(__name__, 'testing_resources'))
MANIFEST = RESOURCES / "manifest.yml"

@pytest.fixture(scope="function")
def molname_workflows_directory(tmp_path):
    m = pybol.Manifest(str(MANIFEST))
    m.assemble('workflows', tmp_path)
    return tmp_path

@pytest.fixture(scope="function")
def csv_tmp_dir(tmp_path, csv='tmp_csv'):
    m = pybol.Manifest(str(MANIFEST))
    m.assemble('workflows', tmp_path)
    return tmp_path / csv

class TestWorkflowsBase(object):

    @pytest.fixture(scope="function")
    def SM_tmp_dir(self, molname_workflows_directory):
        dirname = molname_workflows_directory
        return dirname

    resname = 'UNK'
        
    def test_directory_paths(self, molname_workflows_directory, SM_tmp_dir, csv_tmp_dir):
        # tests directory_paths function and resname-molname conversion/substitution
        parent_directory = molname_workflows_directory
        csv_save_dir=csv_tmp_dir
        directory_paths = base.directory_paths(parent_directory=parent_directory, csv_save_dir=csv_save_dir)
        assert directory_paths['molecule'][0] == 'SM25'
        assert directory_paths['molecule'][1] == 'SM26'
        
    def test_csv_ directory_iteration(self, molname_workflows_directory, SM_tmp_dir, csv_tmp_dir):
        parent_directory = molname_workflows_directory
        csv_save_dir=csv_tmp_dir
        directory_paths = base.directory_paths(csv=f'{csv_save_dir}/dir_paths.csv')

        # change resname to match topology (every SAMPL7 resname is 'UNK')
        # only necessary for this dataset, not necessary for normal use
        directory_paths['resname'] = 'UNK'

        base.directory_iteration(directory_paths, df_save_dir=SM_tmp_dir, solvents=('water',),
                                 ensemble_analysis='DihedralAnalysis')

        assert (SM_tmp_dir / 'SM25' / 'SM25_full_df.csv.bz2').exists()
        assert (SM_tmp_dir / 'SM26' / 'SM26_full_df.csv.bz2').exists()
