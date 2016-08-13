import mdpow.equil
import tempdir as td
import os
import shutil

class TestSolvation(object):

    sims = {"water" : mdpow.equil.WaterSimulation,
            "octanol" : mdpow.equil.OctanolSimulation,
            "cyclohexane" : mdpow.equil.CyclohexaneSimulation,
    }

    def setup(self):
        self.tmpdir = td.TempDir()
        self.old_path = os.getcwd()
        self.resources = self.old_path + "/mdpow/tests/testing_resources"
        self.solvation_paths = self.resources + "/stages/solvation/water/benzene"
        
        files = ['benzene.pdb','benzene.itp']
        for f in files:
            orig = '{0}/molecules/benzene/{1}'.format(self.resources,f)
            shutil.copy(orig, self.tmpdir.name)

    def _test_solvation(self, solvent):
        os.chdir(self.tmpdir.name)
        try:
            if isinstance(solvent, list):
                for sol in solvent:
                    S = self.sims[sol](molecule='BNZ')
                    S.topology(itp='benzene.itp')
                    S.solvate(struct='benzene.pdb')
            else:    
                S = self.sims[solvent](molecule='BNZ')
                S.topology(itp='benzene.itp')
                S.solvate(struct='benzene.pdb')
            assert 1
        except Exception:
            assert 0
        finally:
            os.chdir(self.old_path)

    def test_solvation_water(self):
        return self._test_solvation('water')

    def test_solvation_octanol(self):
        return self._test_solvation('octanol')

    def test_solvation_cyclohexane(self):
        return self._test_solvation('cyclohexane')

    def test_all_solvents(self):
        return self._test_solvation([k for k in self.sims])
