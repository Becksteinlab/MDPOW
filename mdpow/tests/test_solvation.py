import mdpow.equil
import tempdir as td
import os

class TestSolvation(object):
    sim_filename = ".simulation"
    sims = {"water" : mdpow.equil.WaterSimulation,
            "octanol" : mdpow.equil.OctanolSimulation,
            "cyclohexane" : mdpow.equil.CyclohexaneSimulation}

    def setup(self):
        self.tmpdir = td.TempDir()
        self.old_path = os.getcwd()
        print(self.tmpdir.name)
        print(self.old_path)
        assert 0

    def _test_solvation(structure,solvent):
        raise NotImplementedError

    def test_solvation_water(structure):
        return _test_solvation(structure,'water')

    def test_solvation_octanol(structure):
        return _test_solvation(structure,'octanol')

    def test_solvation_cyclohexane(structure):
        return _test_solvation(structure,'cyclohexane')

