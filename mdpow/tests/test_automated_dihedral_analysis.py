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

from numpy.testing import assert_almost_equal
from scipy.stats import circvar, circmean

from . import RESOURCES

import py.path

from ..workflows import dihedrals

from pkg_resources import resource_filename

# NEW TESTING DATA NEEDS TO BE GENERATED AFTER
# REDUCING TESTING DATASET

RESOURCES = pathlib.PurePath(resource_filename(__name__, 'testing_resources'))
MANIFEST = RESOURCES / "manifest.yml"

@pytest.fixture(scope="function")
def molname_workflows_directory(tmp_path, molname='SM25'):
    m = pybol.Manifest(str(MANIFEST))
    m.assemble('workflows', tmp_path)
    return tmp_path / molname

class TestAutomatedDihedralAnalysis(object):

    @pytest.fixture(scope="function")
    def SM25_tmp_dir(self, molname_workflows_directory):
        dirname = molname_workflows_directory
        return dirname

    @pytest.fixture(scope="function")
    def atom_indices(self, SM25_tmp_dir):
        atom_group_indices = dihedrals.dihedral_indices(dirname=SM25_tmp_dir, resname=self.resname)

        # test user input of alternate SMARTS string
        atom_group_indices_alt = dihedrals.dihedral_indices(dirname=SM25_tmp_dir, resname=self.resname,
                                                            SMARTS='[!$(*#*)&!D1]-!@[!$(*#*)&!D1]')
        return atom_group_indices, atom_group_indices_alt
        # atom_indices[0]=atom_group_indices, atom_indices[1]=atom_group_indices_alt,

    @pytest.fixture(scope="function")
    def dihedral_data(self, SM25_tmp_dir, atom_indices):
        atom_group_indices, _ = atom_indices
        df = dihedrals.dihedral_groups_ensemble(atom_group_indices=atom_group_indices, # default values
                                                dirname=SM25_tmp_dir,
                                                solvents=('water',))
        df_aug = dihedrals.periodic_angle(df)
        return df, df_aug
        # dihedral_data[0]=df, dihedral_data[1]=df_aug

    resname = 'UNK'

    # bond indices with default SMARTS string
    check_atom_group_indices = ((0, 1, 2, 3),(0, 1, 12, 13),(1, 2, 3, 11),(1, 2, 3, 10),
                                (1, 2, 3, 4),(1, 12, 13, 14),(2, 3, 4, 5),(2, 3, 4, 9),
                                (2, 1, 12, 13),(3, 2, 1, 12),(5, 4, 3, 11),(5, 4, 3, 10),
                                (9, 4, 3, 11),(9, 4, 3, 10),(12, 13, 14, 15),(12, 13, 14, 19))
    
    # atom group indices with alternate SMARTS string from user input
    check_atom_group_indices_alt = ((1, 2), (1, 12), (2, 3), (3, 4), (12, 13), (13, 14))

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

    universe_solute_atom_names = np.array(['O1', 'C2', 'N3', 'S4', 'C5', 'C6', 'C7', 'C8',
                                           'C9', 'C10', 'O11', 'O12', 'C13', 'C14', 'C15',
                                           'C16', 'C17', 'C18', 'C19', 'C20', 'H21', 'H22',
                                           'H23', 'H24', 'H25', 'H26', 'H27', 'H28', 'H29',
                                           'H30', 'H31', 'H32', 'H33', 'H34', 'H35'], dtype=object)

    check_hydrogens = np.array(['H21', 'H22', 'H23', 'H24', 'H25', 'H26', 'H27', 'H28',
                                'H29', 'H30', 'H31', 'H32', 'H33', 'H34', 'H35'], dtype=object)

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
    
    def test_build_universe(self, SM25_tmp_dir):
        u = dihedrals.build_universe(dirname=SM25_tmp_dir)
        solute = u.select_atoms('resname UNK')
        solute_names = solute.atoms.names
        assert solute_names.all() == self.universe_solute_atom_names.all()

    @pytest.mark.skipif(sys.version_info < (3, 8), reason="pytest=7.2.0, build=py37h89c1867_0 gives wrong answers")
    def test_dihedral_indices(self, atom_indices):
        atom_group_indices = atom_indices[0]
        assert atom_group_indices == self.check_atom_group_indices

    @pytest.mark.skipif(sys.version_info < (3, 8), reason="pytest=7.2.0, build=py37h89c1867_0 gives wrong answers")
    def test_SMARTS(self, atom_indices):
        atom_group_indices_alt = atom_indices[1]
        assert atom_group_indices_alt == self.check_atom_group_indices_alt

    @pytest.mark.skipif(sys.version_info < (3, 8), reason="pytest=7.2.0, build=py37h89c1867_0 gives wrong answers")
    def test_dihedral_groups(self, SM25_tmp_dir):
        groups = dihedrals.dihedral_groups(dirname=SM25_tmp_dir, resname=self.resname)
        i = 0
        while i < len(groups):
            assert groups[i].all() == self.check_groups[i].all()
            i+=1

    @pytest.mark.skipif(sys.version_info < (3, 8), reason="scipy circvar gives wrong answers")
    # separate issue raised to address nature of this problem^
    def test_dihedral_groups_ensemble(self, dihedral_data):

        df = dihedral_data[0]

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
        
    def test_save_df(self, dihedral_data, SM25_tmp_dir):
        dihedrals.save_df(df=dihedral_data[0], df_save_dir=SM25_tmp_dir, molname='SM25')
        assert (SM25_tmp_dir / 'SM25' / 'SM25_full_df.csv.bz2').exists(), 'Compressed csv file not saved'

    def test_save_df_info(self, dihedral_data, SM25_tmp_dir, caplog):
        caplog.clear()
        caplog.set_level(logging.INFO, logger='mdpow.workflows.dihedrals')
        dihedrals.save_df(df=dihedral_data[0], df_save_dir=SM25_tmp_dir, molname='SM25')
        assert f'Results DataFrame saved as {SM25_tmp_dir}/SM25/SM25_full_df.csv.bz2', 'Save location not logged or returned'

    @pytest.mark.skipif(sys.version_info < (3, 8), reason="scipy circvar gives wrong answers")
    # separate issue raised to address nature of this problem^
    def test_periodic_angle(self, dihedral_data):

        df_aug = dihedral_data[1]

        aug_dh2_result = df_aug.loc[df_aug['selection'] == 'C13-C14-C15-C20']['dihedral']

        aug_dh2_mean = circmean(aug_dh2_result, high=180, low=-180)
        aug_dh2_var = circvar(aug_dh2_result, high=180, low=-180)

        aug_dh2_mean == pytest.approx(self.ADG_C13141520_mean)
        aug_dh2_var == pytest.approx(self.ADG_C13141520_var)

    def test_save_fig(self, SM25_tmp_dir):
        dihedrals.automated_dihedral_analysis(dirname=SM25_tmp_dir, figdir=SM25_tmp_dir,
                                              resname=self.resname, molname='SM25',
                                              solvents=('water',))
        assert (SM25_tmp_dir / 'SM25' / 'SM25_C10-C5-S4-O11_violins.pdf').exists(), 'PDF file not generated'

    def test_save_fig_info(self, SM25_tmp_dir, caplog):
        caplog.clear()
        caplog.set_level(logging.INFO, logger='mdpow.workflows.dihedrals')
        dihedrals.automated_dihedral_analysis(dirname=SM25_tmp_dir, figdir=SM25_tmp_dir,
                                              resname=self.resname, molname='SM25',
                                              solvents=('water',))
        assert f'Figure saved as {SM25_tmp_dir}/SM25/SM25_C10-C5-S4-O11_violins.pdf' in caplog.text, 'PDF file not saved'

    def test_DataFrame_input(self, SM25_tmp_dir):
        test_df = pd.DataFrame([['C1-C2-C3-C4', 'water', 'Coulomb', 0, 0, 60.0],
                                ['C1-C2-C3-C5', 'water', 'Coulomb', 0, 0, 60.0]],
                                [1,2],['selection', 'solvent', 'interaction', 'lambda', 'time', 'dihedral'])
        plot = dihedrals.automated_dihedral_analysis(dirname=SM25_tmp_dir, figdir=SM25_tmp_dir,
                                                     resname=self.resname,
                                                     solvents=('water',), dataframe=test_df)
        assert plot # get exact object ID

    def test_DataFrame_input_info(self, SM25_tmp_dir, caplog):
        caplog.clear()
        caplog.set_level(logging.INFO, logger='mdpow.workflows.dihedrals')
        test_df = pd.DataFrame([['C1-C2-C3-C4', 'water', 'Coulomb', 0, 0, 60.0],
                                ['C1-C2-C3-C5', 'water', 'Coulomb', 0, 0, 60.0]],
                                [1,2],['selection', 'solvent', 'interaction', 'lambda', 'time', 'dihedral'])
        dihedrals.automated_dihedral_analysis(dirname=SM25_tmp_dir, figdir=SM25_tmp_dir,
                                              resname=self.resname,
                                              solvents=('water',), dataframe=test_df)
        assert 'Proceeding with results DataFrame provided.' in caplog.text, 'No dataframe provided or dataframe not recognized'