import os
import shutil

from gromacs.utilities import in_dir
import tempdir as td

import mdpow.equil

class TestSolvation(object):
    sims = {"water" : mdpow.equil.WaterSimulation,
            "octanol" : mdpow.equil.OctanolSimulation,
            "cyclohexane" : mdpow.equil.CyclohexaneSimulation,
    }

    def setup(self):
        self.tmpdir = td.TempDir()
        self.old_path = os.getcwd()
        self.resources = os.path.join(
            self.old_path, 'mdpow', 'tests', 'testing_resources')
        self.solvation_paths = os.path.join(
            self.resources, 'stages', 'solvation', 'water', 'benzene')

        # TODO replace by using manifest
        """
        m = manifest.Manifest('manifest.yml')
        m.assemble('base')
        """

        files = ['benzene.pdb','benzene.itp']
        for f in files:
            orig = os.path.join(self.resources, 'molecules', 'benzene', f)
            shutil.copy(orig, self.tmpdir.name)

    def teardown(self):
        self.tmpdir.dissolve()

    def _test_solvation(self, solvent):
        with in_dir(self.tmpdir.name, create=False):
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
            except Exception:
                raise AssertionError('Solvation failed.')

    def test_solvation_water(self):
        return self._test_solvation('water')

    def test_solvation_octanol(self):
        return self._test_solvation('octanol')

    def test_solvation_cyclohexane(self):
        return self._test_solvation('cyclohexane')

    def test_all_solvents(self):
        return self._test_solvation([k for k in self.sims])
