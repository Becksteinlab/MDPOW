'''from __future__ import absolute_import
# review and remove remaning unnecessary code 
import os
import sys
import yaml
import pybol
import pytest
import pathlib

import scipy
import numpy as np
import pandas as pd

import rdkit
from rdkit import Chem

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
def molname_workflows_directory(tmp_path, molname='SM25'):
    m = pybol.Manifest(str(RESOURCES / 'manifest.yml'))
    m.assemble('workflows', tmp_path)
    return tmp_path / molname

class TestAutomatedDihedralAnalysis(object):

    @pytest.fixture(scope="function")
    def SM25_tmp_dir(self, molname_workflows_directory):
        dirname = molname_workflows_directory
        return dirname

    @pytest.fixture(scope="function")
    def gen_data(self, SM25_tmp_dir):
        atom_group_indices = dihedrals.dihedral_indices(dirname=SM25_tmp_dir, resname=self.resname)
        # test user input of SMARTS string
        atom_group_indices_alt = dihedrals.dihedral_indices(dirname=SM25_tmp_dir, resname=self.resname,
                                                            SMARTS='[!$(*#*)&!D1]-!@[!$(*#*)&!D1]')
        df = dihedrals.dihedral_groups_ensemble(atom_group_indices=atom_group_indices, dirname=SM25_tmp_dir,
                                                solvents=('water',))
        df_aug = dihedrals.periodic_angle(df)
        return atom_group_indices, atom_group_indices_alt, df, df_aug
        # gen_data[0]=atom_group_indices, gen_data[1]=atom_group_indices_alt,
        # gen_data[2]=df, gen_data[3]=df_aug

    resname = 'UNK'
        
    def test_base(self, molname_workflows_directory, SM25_tmp_dir):
        # tests directory_paths function and resname-molname conversion
        parent_directory = molname_workflows_directory
        directory_paths = base.directory_paths(parent_directory=parent_directory)
        assert (directory_paths['molecule'] == 'SM25').any() == True
        
    def test_directory_iteration(self, molname_workflows_directory, SM25_tmp_dir):
        parent_directory = molname_workflows_directory
        directory_paths = base.directory_paths(parent_directory=parent_directory)
        directory_paths = directory_paths
        # change resname to match topology (every SAMPL7 resname is 'UNK')
        # only necessary for this dataset, not necessary for normal use
        for f in directory_paths:
            directory_paths.resname = 'UNK'
        base.directory_iteration(directory_paths, df_save_dir=SM25_tmp_dir, solvents=('water',),
                                 ensemble_analysis='DihedralAnalysis', padding=45, width=0.9,
                                 SMARTS='[!#1]~[!$(*#*)&!D1]-!@[!$(*#*)&!D1]~[!#1]')
        assert SM25_tmp_dir / 'SM25' / 'SM25_full_df.bz2'
        assert SM25_tmp_dir / 'SM26' / 'SM26_full_df.bz2'
        '''