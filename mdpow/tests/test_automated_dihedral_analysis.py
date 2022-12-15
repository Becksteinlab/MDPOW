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

from numpy.testing import assert_almost_equal
from scipy.stats import circvar, circmean

from . import RESOURCES

import py.path

#from ..analysis.ensemble import Ensemble, EnsembleAnalysis, EnsembleAtomGroup
#from ..analysis.dihedral import DihedralAnalysis
from ..analysis.workflows import dihedrals as ada

from pkg_resources import resource_filename

RESOURCES = py.path.local(resource_filename(__name__, 'testing_resources'))
MANIFEST = RESOURCES.join("manifest.yml")

@pytest.fixture(scope="function")
def molname_FEP_directory(tmp_path, molname='SM25'):
    tmpdir = tmp_path
    m = pybol.Manifest(str(RESOURCES / 'manifest.yml'))
    m.assemble('FEP', tmpdir)
    return tmpdir / molname

class TestAutomatedDihedralAnalysis(object):

    @pytest.fixture(scope="function")
    def SM25_tmp_dir(self, molname_FEP_directory):
        datadir = molname_FEP_directory
        return datadir

    resname = 'UNK'

    check_bonds = ((0, 1, 2, 3),(0, 1, 12, 13),(1, 2, 3, 11),(1, 2, 3, 10),
                   (1, 2, 3, 4),(1, 12, 13, 14),(2, 3, 4, 5),(2, 3, 4, 9),
                   (2, 1, 12, 13),(3, 2, 1, 12),(5, 4, 3, 11),(5, 4, 3, 10),
                   (9, 4, 3, 11),(9, 4, 3, 10),(12, 13, 14, 15),(12, 13, 14, 19))

    check_groups = ('[array([\'O1\', \'C2\', \'N3\', \'S4\'], dtype=object), '
                    'array([\'O1\', \'C2\', \'C13\', \'C14\'], dtype=object), '
                    'array([\'C2\', \'N3\', \'S4\', \'O12\'], dtype=object), '
                    'array([\'C2\', \'N3\', \'S4\', \'O11\'], dtype=object), '
                    'array([\'C2\', \'N3\', \'S4\', \'C5\'], dtype=object), '
                    'array([\'C2\', \'C13\', \'C14\', \'C15\'], dtype=object), '
                    'array([\'N3\', \'S4\', \'C5\', \'C6\'], dtype=object), '
                    'array([\'N3\', \'S4\', \'C5\', \'C10\'], dtype=object), '
                    'array([\'N3\', \'C2\', \'C13\', \'C14\'], dtype=object), '
                    'array([\'S4\', \'N3\', \'C2\', \'C13\'], dtype=object), '
                    'array([\'C6\', \'C5\', \'S4\', \'O12\'], dtype=object), '
                    'array([\'C6\', \'C5\', \'S4\', \'O11\'], dtype=object), '
                    'array([\'C10\', \'C5\', \'S4\', \'O12\'], dtype=object), '
                    'array([\'C10\', \'C5\', \'S4\', \'O11\'], dtype=object), '
                    'array([\'C13\', \'C14\', \'C15\', \'C16\'], dtype=object), '
                    'array([\'C13\', \'C14\', \'C15\', \'C20\'], dtype=object)]')

    DG_O1C2N3S4_mean = -0.17689226523782509
    DG_O1C2N3S4_var = 0.03409665763502634

    DG_C13141520_mean = 90.0308806327842
    DG_C13141520_var = 0.8807194366657267

    ADG_C13141520_mean = 93.22450745932701
    ADG_C13141520_var = 0.8807224235045962

    def test_dihedral_indices(self, SM25_tmp_dir):
        bonds = ada.dihedral_indices(datadir=SM25_tmp_dir, resname=self.resname)
        assert bonds == self.check_bonds

    '''def test_dihedral_indices_error_exception(self, SM25_tmp_dir):
        with pytest.raises(AttributeError) as excinfo:
            ada.dihedral_indices(datadir=SM25_tmp_dir, resname=self.resname)
        assert excinfo.value is AttributeError'''
    # needs an example topology without hydrogen

    def test_dihedral_groups(self, SM25_tmp_dir):
        groups = ada.dihedral_groups(datadir=SM25_tmp_dir, resname=self.resname)
        assert str(groups) == self.check_groups

    @pytest.mark.skipif(sys.version_info < (3, 8), reason="scipy circvar gives wrong answers")
    def test_dihedral_groups_ensemble(self, SM25_tmp_dir):
        bonds = ada.dihedral_indices(datadir=SM25_tmp_dir, resname=self.resname)
        df = ada.dihedral_groups_ensemble(bonds=bonds, datadir=SM25_tmp_dir, solvents=('water',))

        dh1_result = df.loc[df['selection'] == 'O1-C2-N3-S4']['dihedral']
        dh1_mean = circmean(dh1_result, high=180, low=-180)
        dh1_var = circvar(dh1_result, high=180, low=-180)

        dh2_result = df.loc[df['selection'] == 'C13-C14-C15-C20']['dihedral']
        dh2_mean = circmean(dh2_result, high=180, low=-180)
        dh2_var = circvar(dh2_result, high=180, low=-180)


        assert_almost_equal(dh1_mean, self.DG_O1C2N3S4_mean, 6)
        assert_almost_equal(dh1_var, self.DG_O1C2N3S4_var, 6)

        assert_almost_equal(dh2_mean, self.DG_C13141520_mean, 6)
        assert_almost_equal(dh2_var, self.DG_C13141520_var, 6)

    @pytest.mark.skipif(sys.version_info < (3, 8), reason="scipy circvar gives wrong answers") 
    def test_periodic_angle(self, SM25_tmp_dir):
        bonds = ada.dihedral_indices(datadir=SM25_tmp_dir, resname=self.resname)
        df = ada.dihedral_groups_ensemble(bonds=bonds, datadir=SM25_tmp_dir, solvents=('water',))

        dh2_result = df.loc[df['selection'] == 'C13-C14-C15-C20']['dihedral']
        dh2_mean = circmean(dh2_result, high=180, low=-180)
        dh2_var = circvar(dh2_result, high=180, low=-180)

        df_aug = ada.periodic_angle(df)

        aug_dh2_result = df_aug.loc[df_aug['selection'] == 'C13-C14-C15-C20']['dihedral']
        aug_dh2_mean = circmean(aug_dh2_result, high=180, low=-180)
        aug_dh2_var = circvar(aug_dh2_result, high=180, low=-180)

        assert_almost_equal(aug_dh2_mean, self.ADG_C13141520_mean, 6)
        assert_almost_equal(aug_dh2_var, self.ADG_C13141520_var, 6)

        assert dh2_mean != aug_dh2_mean
