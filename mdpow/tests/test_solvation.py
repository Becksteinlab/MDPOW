import mdpow.equil
import tempdir as td
import os
import shutil

class TestSolvation(object):
    sim_filename = ".simulation"
    sims = {"water" : mdpow.equil.WaterSimulation,
            "octanol" : mdpow.equil.OctanolSimulation,
            "cyclohexane" : mdpow.equil.CyclohexaneSimulation,
           }

    def setup(self):
        self.tmpdir = td.TempDir()
        self.old_path = os.getcwd()
        self.resources = self.old_path + "/mdpow/tests/testing_resources"
        self.solvation_paths = self.resources + "/stages/solvation/water/benzene"
        [shutil.copy('{0}/molecules/benzene/{1}'.format(self.resources, x), self.tmpdir.name) for x in ['benzene.pdb','benzene.itp']]

    def _test_solvation(self, solvent):
        os.chdir(self.tmpdir.name)
        try:
            S = self.sims[solvent](molecule='BNZ')
            S.topology(itp='benzene.itp')
            S.solvate(struct='benzene.pdb')
            assert 1
        finally:
            os.chdir(self.old_path)

    def test_solvation_water(self):
        return self._test_solvation('water')

    def test_solvation_octanol(self):
        return self._test_solvation('octanol')

    def test_solvation_cyclohexane(self):
        return self._test_solvation('cyclohexane')
