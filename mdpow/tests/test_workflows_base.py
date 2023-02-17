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

#@pytest.fixture(scope="function")
#def csv_tmp_dir(tmp_path, csv='tmp_csv'):
#    m = pybol.Manifest(str(MANIFEST))
#    m.assemble('workflows', tmp_path)
#    return tmp_path / csv

class TestWorkflowsBase(object):

    @pytest.fixture(scope="function")
    def SM_tmp_dir(self, molname_workflows_directory):
        dirname = molname_workflows_directory
        return dirname

    resname = 'UNK'
        
    def test_directory_paths(self, SM_tmp_dir):
        directory_paths = base.directory_paths(parent_directory=SM_tmp_dir)
        assert directory_paths['molecule'][0] == 'SM25'
        assert directory_paths['molecule'][1] == 'SM26'
        
    def test_directory_iteration(self, SM_tmp_dir, caplog):
        directory_paths = base.directory_paths(parent_directory=SM_tmp_dir)
        # change resname to match topology (every SAMPL7 resname is 'UNK')
        # only necessary for this dataset, not necessary for normal use
        directory_paths['resname'] = 'UNK'

        base.directory_iteration(directory_paths, solvents=('water',),
                                 ensemble_analysis='DihedralAnalysis')

        assert 'all analyses completed' in caplog.text, 'automated_dihedral_analysis did not iteratively run to completion for the provided project'
