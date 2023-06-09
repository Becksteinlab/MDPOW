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

from ..workflows import dihedrals

from pkg_resources import resource_filename

RESOURCES = pathlib.PurePath(resource_filename(__name__, 'testing_resources'))
MANIFEST = RESOURCES / "manifest.yml"

resname = "UNK"
molname = "SM25"

@pytest.fixture
def molname_workflows_directory(tmp_path, molname='SM25'):
    m = pybol.Manifest(str(MANIFEST))
    m.assemble('workflows', tmp_path)
    return tmp_path / molname

class TestAutomatedDihedralAnalysis(object):

    @pytest.fixture
    def SM25_tmp_dir(self, molname_workflows_directory):
        dirname = molname_workflows_directory
        return dirname

    @pytest.fixture
    def mol_sol_data(self, SM25_tmp_dir):
        u = dihedrals.build_universe(dirname=SM25_tmp_dir)
        mol, solute = dihedrals.rdkit_conversion(u=u, resname=resname)
        return mol, solute

    @pytest.fixture
    def atom_indices(self, mol_sol_data):
        mol, _ = mol_sol_data
        atom_group_indices = dihedrals.get_atom_indices(mol=mol)

        # testing optional user input of alternate SMARTS string
        # for automated dihedral atom group selection
        atom_group_indices_alt = dihedrals.get_atom_indices(mol=mol, SMARTS='[!$(*#*)&!D1]-!@[!$(*#*)&!D1]')
        return atom_group_indices, atom_group_indices_alt
        # fixture output, tuple:
        # atom_indices[0]=atom_group_indices
        # atom_indices[1]=atom_group_indices_alt

    @pytest.fixture
    def bond_indices(self, mol_sol_data, atom_indices):
        mol, _ = mol_sol_data
        aix, _ = atom_indices
        bond_indices = dihedrals.get_bond_indices(mol=mol, atom_indices=aix)
        return bond_indices

    @pytest.fixture
    def dihedral_groups(self, mol_sol_data, atom_indices):
        _, solute = mol_sol_data
        aix, _ = atom_indices
        dihedral_groups = dihedrals.get_dihedral_groups(solute=solute, atom_indices=aix)
        return dihedral_groups

    @pytest.fixture
    def dihedral_data(self, SM25_tmp_dir, atom_indices):
        atom_group_indices, _ = atom_indices
        df = dihedrals.dihedral_groups_ensemble(atom_indices=atom_group_indices,
                                                dirname=SM25_tmp_dir,
                                                solvents=('water',))
        df_aug = dihedrals.periodic_angle_padding(df)
        return df, df_aug
        # fixture output, tuple:
        # dihedral_data[0]=df
        # dihedral_data[1]=df_aug

    # tuple-tuples of dihedral atom group indices
    # collected using mdpow.workflows.dihedrals.SMARTS_DEFAULT
    check_atom_group_indices = ((0, 1, 2, 3),(0, 1, 12, 13),(1, 2, 3, 11),(1, 2, 3, 10),
                                (1, 2, 3, 4),(1, 12, 13, 14),(2, 3, 4, 5),(2, 3, 4, 9),
                                (2, 1, 12, 13),(3, 2, 1, 12),(5, 4, 3, 11),(5, 4, 3, 10),
                                (9, 4, 3, 11),(9, 4, 3, 10),(12, 13, 14, 15),(12, 13, 14, 19))

    check_bond_indices = ((8, 4, 2), (8, 9, 14), (4, 2, 0), (4, 2, 1),
                          (4, 2, 3), (9, 14, 21), (2, 3, 6), (2, 3, 7),
                          (4, 9, 14), (2, 4, 9), (6, 3, 0), (6, 3, 1),
                          (7, 3, 0), (7, 3, 1), (14, 21, 25), (14, 21, 26))

    check_ab_pairs = {'O1-C2-N3-S4': ((0, 1, 2, 3), (8, 4, 2)),
                      'O1-C2-C13-C14': ((0, 1, 12, 13), (8, 9, 14)),
                      'C2-N3-S4-O12': ((1, 2, 3, 11), (4, 2, 0)),
                      'C2-N3-S4-O11': ((1, 2, 3, 10), (4, 2, 1)),
                      'C2-N3-S4-C5': ((1, 2, 3, 4), (4, 2, 3)),
                      'C2-C13-C14-C15': ((1, 12, 13, 14), (9, 14, 21)),
                      'N3-S4-C5-C6': ((2, 3, 4, 5), (2, 3, 6)),
                      'N3-S4-C5-C10': ((2, 3, 4, 9), (2, 3, 7)),
                      'N3-C2-C13-C14': ((2, 1, 12, 13), (4, 9, 14)),
                      'S4-N3-C2-C13': ((3, 2, 1, 12), (2, 4, 9)),
                      'C6-C5-S4-O12': ((5, 4, 3, 11), (6, 3, 0)),
                      'C6-C5-S4-O11': ((5, 4, 3, 10), (6, 3, 1)),
                      'C10-C5-S4-O12': ((9, 4, 3, 11), (7, 3, 0)),
                      'C10-C5-S4-O11': ((9, 4, 3, 10), (7, 3, 1)),
                      'C13-C14-C15-C16': ((12, 13, 14, 15), (14, 21, 25)),
                      'C13-C14-C15-C20': ((12, 13, 14, 19), (14, 21, 26))}

    # tuple-tuples of dihedral atom group indices
    # collected using alternate SMARTS input (explicitly defined)
    # see: fixture - atom_indices().atom_group_indices_alt
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

    # pre 'angle padding' - scipy.stats for
    # dihedral atom group: O1-C2-N3-S4
    DG_O1C2N3S4_mean = -0.5933808752787115
    DG_O1C2N3S4_var = 0.031351024457919485

    # pre 'angle padding' - scipy.stats for
    # dihedral atom group: C13-C14-C15-C20
    DG_C13141520_mean = 89.22382649857468
    DG_C13141520_var = 0.8753980937068645

    # post 'angle padding' - scipy.stats for
    # dihedral atom group: C13-C14-C15-C20
    # 'A' = 'augmented', referencing those
    # results included in 'df_aug'
    ADG_C13141520_mean = 91.71943996962284
    ADG_C13141520_var = 0.8773028474908289
    
    def test_build_universe(self, SM25_tmp_dir):
        u = dihedrals.build_universe(dirname=SM25_tmp_dir)
        solute = u.select_atoms('resname UNK')
        solute_names = solute.atoms.names
        assert solute_names.all() == self.universe_solute_atom_names.all()

    # the following 'reason' affects every downstream function that relies
    # on the atom indices returned for dihedral atom group selections
    # issue raised (#239) to identify and resolve exact package/version responsible
    @pytest.mark.skipif(sys.version_info < (3, 8), reason='pytest=7.2.0, build=py37h89c1867_0, '
                       'returns incorrect atom_indices for dihedral atom group selections')
    def test_dihedral_indices(self, atom_indices):
        atom_group_indices, _ = atom_indices
        assert atom_group_indices == self.check_atom_group_indices

    # the following 'reason' affects every downstream function that relies
    # on the atom indices returned for dihedral atom group selections
    # issue raised (#239) to identify and resolve exact package/version responsible
    @pytest.mark.skipif(sys.version_info < (3, 8), reason='pytest=7.2.0, build=py37h89c1867_0, '
                       'returns incorrect atom_indices for dihedral atom group selections')
    def test_bond_indices(self, bond_indices):
        bix = bond_indices
        assert bix == self.check_bond_indices

    # the following 'reason' affects every downstream function that relies
    # on the atom indices returned for dihedral atom group selections
    # issue raised (#239) to identify and resolve exact package/version responsible
    @pytest.mark.skipif(sys.version_info < (3, 8), reason='pytest=7.2.0, build=py37h89c1867_0, '
                       'returns incorrect atom_indices for dihedral atom group selections')
    def test_SMARTS(self, atom_indices):
        _, atom_group_indices_alt = atom_indices
        assert atom_group_indices_alt == self.check_atom_group_indices_alt

    # the following 'reason' affects every downstream function that relies
    # on the atom indices returned for dihedral atom group selections
    # issue raised (#239) to identify and resolve exact package/version responsible
    @pytest.mark.skipif(sys.version_info < (3, 8), reason='pytest=7.2.0, build=py37h89c1867_0, '
                       'returns incorrect atom_indices for dihedral atom group selections')
    def test_dihedral_groups(self, dihedral_groups):
        groups = dihedral_groups
        i = 0
        while i < len(groups):
            assert groups[i].all() == self.check_groups[i].all()
            i+=1

    # the following 'reason' affects every downstream function that relies
    # on the atom indices returned for dihedral atom group selections
    # issue raised (#239) to identify and resolve exact package/version responsible
    @pytest.mark.skipif(sys.version_info < (3, 8), reason='pytest=7.2.0, build=py37h89c1867_0, '
                       'returns incorrect atom_indices for dihedral atom group selections')
    def test_ab_pairs(self, atom_indices, bond_indices, dihedral_groups):
        aix, _ = atom_indices
        bix = bond_indices
        groups = dihedral_groups
        ab_pairs = dihedrals.get_paired_indices(atom_indices=aix, bond_indices=bix,
                                                dihedral_groups=groups)
        assert ab_pairs == self.check_ab_pairs

    # the following 'reason' affects every downstream function that relies
    # on the atom indices returned for dihedral atom group selections
    # issue raised (#239) to identify and resolve exact package/version responsible
    @pytest.mark.skipif(sys.version_info < (3, 8), reason='pytest=7.2.0, build=py37h89c1867_0, '
                       'returns incorrect atom_indices for dihedral atom group selections')
    def test_dihedral_groups_ensemble(self, dihedral_data):

        df, _ = dihedral_data

        dh1_result = df.loc[df['selection'] == 'O1-C2-N3-S4']['dihedral']
        dh1_mean = circmean(dh1_result, high=180, low=-180)
        dh1_var = circvar(dh1_result, high=180, low=-180)

        dh2_result = df.loc[df['selection'] == 'C13-C14-C15-C20']['dihedral']
        dh2_mean = circmean(dh2_result, high=180, low=-180)
        dh2_var = circvar(dh2_result, high=180, low=-180)

        # default relative tolerance for pytest.approx is 1e-6
        dh1_mean == pytest.approx(self.DG_O1C2N3S4_mean)
        dh1_var == pytest.approx(self.DG_O1C2N3S4_var)

        dh2_mean == pytest.approx(self.DG_C13141520_mean)
        dh2_var == pytest.approx(self.DG_C13141520_var)

    def test_save_df(self, dihedral_data, SM25_tmp_dir):
        df, _ = dihedral_data
        dihedrals.save_df(df=df, df_save_dir=SM25_tmp_dir, resname='UNK', molname='SM25')
        assert (SM25_tmp_dir / 'SM25' / 'SM25_full_df.csv.bz2').exists(), 'Compressed csv file not saved'

    def test_save_df_info(self, dihedral_data, SM25_tmp_dir, caplog):
        df, _ = dihedral_data
        caplog.clear()
        caplog.set_level(logging.INFO, logger='mdpow.workflows.dihedrals')
        dihedrals.save_df(df=df, df_save_dir=SM25_tmp_dir, resname='UNK', molname='SM25')
        assert f'Results DataFrame saved as {SM25_tmp_dir}/SM25/SM25_full_df.csv.bz2' in caplog.text, 'Save location not logged or returned'

    # the following 'reason' affects every downstream function that relies
    # on the atom indices returned for dihedral atom group selections
    # issue raised (#239) to identify and resolve exact package/version responsible
    @pytest.mark.skipif(sys.version_info < (3, 8), reason='pytest=7.2.0, build=py37h89c1867_0, '
                       'returns incorrect atom_indices for dihedral atom group selections')
    def test_periodic_angle_padding(self, dihedral_data):

        _, df_aug = dihedral_data

        aug_dh2_result = df_aug.loc[df_aug['selection'] == 'C13-C14-C15-C20']['dihedral']

        aug_dh2_mean = circmean(aug_dh2_result, high=180, low=-180)
        aug_dh2_var = circvar(aug_dh2_result, high=180, low=-180)

        aug_dh2_mean == pytest.approx(self.ADG_C13141520_mean)
        aug_dh2_var == pytest.approx(self.ADG_C13141520_var)

    # the following 'reason' affects every downstream function that relies
    # on the atom indices returned for dihedral atom group selections
    # issue raised (#239) to identify and resolve exact package/version responsible
    @pytest.mark.skipif(sys.version_info < (3, 8), reason='pytest=7.2.0, build=py37h89c1867_0, '
                       'returns incorrect atom_indices for dihedral atom group selections')
    def test_save_fig(self, SM25_tmp_dir):
        dihedrals.automated_dihedral_analysis(dirname=SM25_tmp_dir, figdir=SM25_tmp_dir,
                                              resname=resname, molname='SM25',
                                              solvents=('water',))
        assert (SM25_tmp_dir / 'SM25' / 'SM25_C10-C5-S4-O11_violins.pdf').exists(), 'PDF file not generated'

    # the following 'reason' affects every downstream function that relies
    # on the atom indices returned for dihedral atom group selections
    # issue raised (#239) to identify and resolve exact package/version responsible
    @pytest.mark.skipif(sys.version_info < (3, 8), reason='pytest=7.2.0, build=py37h89c1867_0, '
                       'returns incorrect atom_indices for dihedral atom group selections')
    def test_save_fig_info(self, SM25_tmp_dir, caplog):
        caplog.clear()
        caplog.set_level(logging.INFO, logger='mdpow.workflows.dihedrals')
        dihedrals.automated_dihedral_analysis(dirname=SM25_tmp_dir, figdir=SM25_tmp_dir,
                                              resname=resname, molname='SM25',
                                              solvents=('water',))
        assert f'Figure saved as {SM25_tmp_dir}/SM25/SM25_C10-C5-S4-O11_violins.pdf' in caplog.text, 'PDF file not saved'

    def test_DataFrame_input(self, SM25_tmp_dir, dihedral_data):
        df, _ = dihedral_data
        dihedrals.automated_dihedral_analysis(dirname=SM25_tmp_dir, figdir=SM25_tmp_dir,
                                              resname=resname, molname=molname,
                                              solvents=('water',), dataframe=df)
        assert (SM25_tmp_dir / 'SM25' / 'SM25_C10-C5-S4-O11_violins.pdf').exists(), 'PDF file not generated'

    def test_DataFrame_input_info(self, SM25_tmp_dir, dihedral_data, caplog):
        caplog.clear()
        caplog.set_level(logging.INFO, logger='mdpow.workflows.dihedrals')
        df, _ = dihedral_data
        dihedrals.automated_dihedral_analysis(dirname=SM25_tmp_dir, figdir=SM25_tmp_dir,
                                              resname=resname, molname=molname,
                                              solvents=('water',), dataframe=df)
        assert 'Proceeding with results DataFrame provided.' in caplog.text, 'No dataframe provided or dataframe not recognized'
