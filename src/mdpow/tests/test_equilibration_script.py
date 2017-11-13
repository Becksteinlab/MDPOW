import os.path
import tempdir as td
import pybol

from gromacs.utilities import in_dir

from mdpow.run import equilibrium_simulation
from mdpow.config import get_configuration

import pkg_resources

TEST_RESOURCES = pkg_resources.resource_filename(
    __name__, 'testing_resources')

class TestEquilibriumScript(object):
    def setup(self):
        self.tmpdir = td.TempDir()
        self.resources = TEST_RESOURCES
        m = pybol.Manifest(os.path.join(self.resources, 'manifest.yml'))
        m.assemble('base', self.tmpdir.name)

    def teardown(self):
        self.tmpdir.dissolve()

    def _run_equil(self, solvent, dirname):
        cfg = get_configuration('runinput.yml')
        self.S = equilibrium_simulation(cfg, solvent, dirname=dirname)

    def test_basic_run(self):
        with in_dir(self.tmpdir.name, create=False):
            try:
                self._run_equil('water','benzene/')
                self._new_structures()
            except Exception as err:
                raise AssertionError('Equilibration simulations failed with exception:\n{0}'.format(str(err)))

    def _new_structures(self):
        assert os.path.exists(
            os.path.join(self.tmpdir.name,
                         'benzene', 'Equilibrium', 'water', 'em', 'em.pdb'))
        assert os.path.exists(
            os.path.join(self.tmpdir.name,
                         'benzene', 'Equilibrium', 'water', 'solvation', 'solvated.gro'))
        assert os.path.exists(
            os.path.join(self.tmpdir.name,
                         'benzene', 'Equilibrium', 'water', 'MD_NPT', 'md.gro'))

