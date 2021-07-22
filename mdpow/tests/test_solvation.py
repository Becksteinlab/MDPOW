import os
import shutil

import mdpow.equil
from gromacs.utilities import in_dir
import gromacs

import pytest

sims = {"water" : mdpow.equil.WaterSimulation,
        "octanol" : mdpow.equil.OctanolSimulation,
        "cyclohexane" : mdpow.equil.CyclohexaneSimulation,
        "wetoctanol" : mdpow.equil.WetOctanolSimulation,
        }

test_file = {"OPLS-AA": 'benzene.itp',
             "CHARMM": 'benzene_charmm.itp',
             "AMBER": 'benzene_amber.itp',
             }

@pytest.fixture
def setup(tmpdir):
    newdir = tmpdir.mkdir('resources')
    old_path = os.getcwd()
    resources = os.path.join(
        old_path, 'mdpow', 'tests', 'testing_resources')
    files = ['benzene.pdb', 'benzene.itp',
             'benzene_charmm.itp', 'benzene_amber.itp']
    for f in files:
        orig = os.path.join(resources, 'molecules', 'benzene', f)
        shutil.copy(orig, newdir.dirname)
    return newdir.dirname

def solvation(setup, solvent, ff='OPLS-AA'):
    itp = test_file[ff]
    with in_dir(setup, create=False):
        try:
            S = sims[solvent](molecule='BNZ', forcefield=ff)
            S.topology(itp=itp)
            S.solvate(struct='benzene.pdb')
        except Exception:
            raise AssertionError('Solvation failed.')

@pytest.mark.parametrize("ff", ['OPLS-AA', 'CHARMM', 'AMBER'])
def test_solvation_water(setup, ff):
    solvation(setup, "water", ff)

@pytest.mark.parametrize("ff", ['OPLS-AA', 'CHARMM', 'AMBER'])
def test_solvation_octanol(setup, ff):
    solvation(setup, "octanol", ff)

def test_solvation_cyclohexane(setup):
    solvation(setup, "cyclohexane")

@pytest.mark.xfail(gromacs.release.startswith('4')
                   or gromacs.release.startswith('5')
                   or gromacs.release.startswith('2016'),
                   reason="GROMACS < 2018 cannot easily work with mixed solvents "
                   "(see issue #111)")
@pytest.mark.parametrize("ff", ['OPLS-AA', 'CHARMM', 'AMBER'])
def test_solvation_wetoctanol(setup, ff):
    solvation(setup, "wetoctanol", ff)
