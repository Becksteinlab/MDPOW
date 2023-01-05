from __future__ import absolute_import

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

#from ..analysis.ensemble import Ensemble, EnsembleAnalysis, EnsembleAtomGroup
#from ..analysis.dihedral import DihedralAnalysis
from ..analysis.workflows import dihedrals as ada
from ..analysis.workflows import base

from pkg_resources import resource_filename

RESOURCES = pathlib.PurePath(resource_filename(__name__, 'testing_resources'))
MANIFEST = RESOURCES / "manifest.yml"

@pytest.fixture(scope="function")
def molname_FEP_directory(tmp_path, molname='SM25'):
    m = pybol.Manifest(str(RESOURCES / 'manifest.yml'))
    m.assemble('FEP', tmp_path)
    return tmp_path / molname

class TestAutomatedDihedralAnalysis(object):

    @pytest.fixture(scope="function")
    def SM25_tmp_dir(self, molname_FEP_directory):
        dirname = molname_FEP_directory
        return dirname

    @pytest.fixture(scope="function")
    def gen_data(self, SM25_tmp_dir):
        bonds = ada.dihedral_indices(dirname=SM25_tmp_dir, resname=self.resname)
        # test user input of SMARTS string, 'user_SMARTS'
        bonds_alt = ada.dihedral_indices(dirname=SM25_tmp_dir, resname=self.resname,
                                         user_SMARTS='[!$(*#*)&!D1]-!@[!$(*#*)&!D1]')
        df = ada.dihedral_groups_ensemble(bonds=bonds, dirname=SM25_tmp_dir, solvents=('water',))
        df_aug = ada.periodic_angle(df)
        return bonds, bonds_alt, df, df_aug
        # gen_data[0]=bonds, gen_data[1]=bonds_alt, gen_data[2]=df, gen_data[3]=df_aug

    resname = 'UNK'

    # bond indices with default SMARTS string
    check_bonds = ((0, 1, 2, 3),(0, 1, 12, 13),(1, 2, 3, 11),(1, 2, 3, 10),
                   (1, 2, 3, 4),(1, 12, 13, 14),(2, 3, 4, 5),(2, 3, 4, 9),
                   (2, 1, 12, 13),(3, 2, 1, 12),(5, 4, 3, 11),(5, 4, 3, 10),
                   (9, 4, 3, 11),(9, 4, 3, 10),(12, 13, 14, 15),(12, 13, 14, 19))
    
    # bond indices with alternate SMARTS string from user input
    check_bonds_alt = ((1, 2), (1, 12), (2, 3), (3, 4), (12, 13), (13, 14))

    check_groups = [np.array(['O1', 'C2', 'N3', 'S4'], dtype=object),
                    np.array(['O1', 'C2', 'C13', 'C14'], dtype=object),
                    np.array(['C2', 'N3', 'S4', 'O12'], dtype=object),
                    np.array(['C2', 'N3', 'S4', 'O11'], dtype=object),
                    np.array(['C2', 'N3', 'S4', 'C5'], dtype=object),
                    np.array(['C2', 'C13', 'C14', 'C15'], dtype=object),
                    np.array(['N3', 'S4', 'C5', 'C6'], dtype=object),
                    np.array(['N3', 'S4', 'C5', 'C10'], dtype=object),
                    np.array(['N3', 'C2', 'C13', 'C14'], dtype=object),
                    np.array(['S4', 'N3', 'C2', 'C13'], dtype=object),
                    np.array(['C6', 'C5', 'S4', 'O12'], dtype=object),
                    np.array(['C6', 'C5', 'S4', 'O11'], dtype=object),
                    np.array(['C10', 'C5', 'S4', 'O12'], dtype=object),
                    np.array(['C10', 'C5', 'S4', 'O11'], dtype=object),
                    np.array(['C13', 'C14', 'C15', 'C16'], dtype=object),
                    np.array(['C13', 'C14', 'C15', 'C20'], dtype=object)]

    # pre-padding dihedral group statistics
    DG_O1C2N3S4_mean = -0.13089566578887002
    DG_O1C2N3S4_var = 0.03473127346296412

    # pre-padding dihedral group statistics
    DG_C13141520_mean = 90.04076626959608
    DG_C13141520_var = 0.8814931303534448

    # post-padding dihedral group statistics
    # 'A' represents 'augmented'
    ADG_C13141520_mean = 93.50126701923381
    ADG_C13141520_var = 0.8815675813248334

    def test_dihedral_indices(self, gen_data):
        bonds = gen_data[0]
        assert bonds == self.check_bonds
        
    def test_user_SMARTS(self, gen_data):
        bonds = gen_data[1]
        assert bonds == self.check_bonds_alt

    '''def test_dihedral_indices_error_exception(self, SM25_tmp_dir):
        with pytest.raises(AttributeError) as excinfo:
            ada.dihedral_indices(dirname=SM25_tmp_dir, resname=self.resname)
        assert excinfo.value is AttributeError'''
    # needs an example topology without hydrogen

    def test_dihedral_groups(self, SM25_tmp_dir):
        groups = ada.dihedral_groups(dirname=SM25_tmp_dir, resname=self.resname)
        i = 0
        while i < len(groups):
            assert groups[i].all() == self.check_groups[i].all()
            i+=1

    @pytest.mark.skipif(sys.version_info < (3, 8), reason="scipy circvar gives wrong answers")
    def test_dihedral_groups_ensemble(self, gen_data):

        df = gen_data[2]

        dh1_result = df.loc[df['selection'] == 'O1-C2-N3-S4']['dihedral']
        dh1_mean = circmean(dh1_result, high=180, low=-180)
        dh1_var = circvar(dh1_result, high=180, low=-180)

        dh2_result = df.loc[df['selection'] == 'C13-C14-C15-C20']['dihedral']
        dh2_mean = circmean(dh2_result, high=180, low=-180)
        dh2_var = circvar(dh2_result, high=180, low=-180)

        # default relative tolerance is 1e-6 for pytest.approx
        dh1_mean == pytest.approx(self.DG_O1C2N3S4_mean)
        dh1_var == pytest.approx(self.DG_O1C2N3S4_var)

        dh2_mean == pytest.approx(self.DG_C13141520_mean)
        dh2_var == pytest.approx(self.DG_C13141520_var)
        
    def test_save_df(self, gen_data, SM25_tmp_dir):
        ada.save_df(df=gen_data[2], df_save_dir=SM25_tmp_dir, molname='SM25')
        assert SM25_tmp_dir / 'SM25' / 'SM25_full_df.bz2'
                                        

    @pytest.mark.skipif(sys.version_info < (3, 8), reason="scipy circvar gives wrong answers") 
    def test_periodic_angle(self, gen_data):

        df_aug = gen_data[3]

        aug_dh2_result = df_aug.loc[df_aug['selection'] == 'C13-C14-C15-C20']['dihedral']

        aug_dh2_mean = circmean(aug_dh2_result, high=180, low=-180)
        aug_dh2_var = circvar(aug_dh2_result, high=180, low=-180)

        aug_dh2_mean == pytest.approx(self.ADG_C13141520_mean)
        aug_dh2_var == pytest.approx(self.ADG_C13141520_var)
        
    def test_base(self, SM25_tmp_dir):
        # tests directory_paths function and resname-molname conversion
        directory_paths = base.directory_paths(parent_directory=RESOURCES)
        assert (directory_paths['molecule'] == 'SM25').any() == True

    def test_save_fig(self, SM25_tmp_dir):
        ada.automated_dihedral_analysis(dirname=SM25_tmp_dir, figdir=SM25_tmp_dir,
                                        resname=self.resname, molname='SM25',
                                        solvents=('water',))
        assert SM25_tmp_dir / 'SM25' / 'SM25_C10-C5-S4-O11_violins.pdf'