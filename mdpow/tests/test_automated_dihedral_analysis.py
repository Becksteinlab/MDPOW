import pytest

from __future__ import absolute_import

from . import tempdir as td

import sys

import py.path

import pybol
import pytest

import numpy as np
import pandas as pd
import os
import pathlib

from rdkit import Chem

import MDAnalysis as mda
from MDAnalysis.topology.guessers import guess_atom_element

from ..analysis.ensemble import Ensemble, EnsembleAnalysis, EnsembleAtomGroup
from ..analysis.dihedral import DihedralAnalysis
from ..analysis.automation.dihedral import automated_dihedral_analysis as ada
from pkg_resources import resource_filename

RESOURCES = py.path.local(resource_filename(__name__, 'testing_resources'))
MANIFEST = RESOURCES.join("manifest.yml")

class TestDihedralAutomation(object):
    
    @pytest.fixture
    def setup(self):
        self.tmpdir = td.TempDir()
        self.m = pybol.Manifest(str(RESOURCES / 'manifest.yml'))
        self.m.assemble('example_FEP/SM25/FEP/', self.tmpdir.name)
        self.Ens = Ensemble(dirname=self.tmpdir.name, solvents=['water'])
        
    def test_dihedral_indices(self):

        topology = self.tmpdir / 'md.tpr'
        trajectory = self.tmpdir / 'md.xtc'
        u = mda.Universe(str(topology), str(trajectory))
        #assert u.atoms == #atoms
        try:
            solute = u.select_atoms(f'resname {resname}')
            mol = solute.convert_to('RDKIT')
        except AttributeError:
            u_aug = add_hydrogens(u)
            solute = u_aug.select_atoms(f'resname {resname}')
            mol = solute.convert_to('RDKIT')
        with.pytest.raises(AttributeError) as excifno:
            #assert excinfo.type == 'AttributeError'

        pattern = Chem.MolFromSmarts('[!#1]~[!$(*#*)&!D1]-!@[!$(*#*)&!D1]~[!#1]')
        bonds = mol.GetSubstructMatches(pattern)
        
        
    def test_dihedral_groups(self):
        
    def test_dihedral_groups_ensemble(self):
        
    def test_periodic_angle(self):
        
    


