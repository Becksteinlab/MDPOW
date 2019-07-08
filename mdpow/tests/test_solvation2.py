import os
import shutil

import mdpow.equil
from gromacs.utilities import in_dir

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

    sims = {"water" : mdpow.equil.WaterSimulation,
            "octanol" : mdpow.equil.OctanolSimulation,
            "cyclohexane" : mdpow.equil.CyclohexaneSimulation,
            #"wetoctanol" : mdpow.equil.WetOctanolSimulation,
    }


@pytest.mark.parametrize("sim", [
mdpow.equil.WaterSimulation,
mdpow.equil.OctanolSimulation,
mdpow.equil.CyclohexaneSimulation,
mdpow.equil.WetOctanolSimulation
])
def test_solvation(setup, sim):
    with in_dir(setup, create=False):
        try:
            S = sim(molecule='BNZ')
            S.topology(itp='benzene.itp')
            S.solvate(struct='benzene.pdb')
        except Exception:
            raise AssertionError('Solvation failed.')
