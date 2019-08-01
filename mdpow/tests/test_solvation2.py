import os
import shutil

import mdpow.equil
from gromacs.utilities import in_dir
import gromacs

import pytest

@pytest.fixture
def setup(tmpdir):
    newdir = tmpdir.mkdir('resources')
    old_path = os.getcwd()
    resources = os.path.join(
        old_path, 'mdpow', 'tests', 'testing_resources')
    files = ['benzene.pdb','benzene.itp']
    for f in files:
        orig = os.path.join(resources, 'molecules', 'benzene', f)
        shutil.copy(orig, newdir.dirname)
    return newdir.dirname

def solvation(setup, solvent):
    sims = {"water" : mdpow.equil.WaterSimulation,
            "octanol" : mdpow.equil.OctanolSimulation,
            "cyclohexane" : mdpow.equil.CyclohexaneSimulation,
            "wetoctanol" : mdpow.equil.WetOctanolSimulation,
            }
    with in_dir(setup, create=False):
        try:
            S = sims[solvent](molecule='BNZ')
            S.topology(itp='benzene.itp')
            S.solvate(struct='benzene.pdb')
        except Exception:
            raise AssertionError('Solvation failed.')

def test_solvation_water(setup):
    solvation(setup, "water")

def test_solvation_octanol(setup):
    solvation(setup, "octanol")

def test_solvation_cyclohexane(setup):
    solvation(setup, "cyclohexane")

@pytest.mark.xfail("gromacs.release().startswith('2019')")
def test_solvation_wetoctanol(setup):
    solvation(setup, "wetoctanol")
