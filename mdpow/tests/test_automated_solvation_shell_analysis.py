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

from . import RESOURCES

import py.path

from ..workflows import solvations

from pkg_resources import resource_filename

# ^review and update these as necessary, currently copied from test_ada

RESOURCES = pathlib.PurePath(resource_filename(__name__, 'testing_resources'))
MANIFEST = RESOURCES / "manifest.yml"

@pytest.fixture(scope="function")
def molname_workflows_directory(tmp_path, molname='SM25'):
    m = pybol.Manifest(str(MANIFEST))
    m.assemble('workflows', tmp_path)
    return tmp_path / molname

class TestAutomatedSolvationShellAnalysis(object):

    @pytest.fixture(scope="function")
    def SM25_tmp_dir(self, molname_workflows_directory):
        dirname = molname_workflows_directory
        return dirname