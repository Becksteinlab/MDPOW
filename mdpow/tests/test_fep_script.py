import tempdir as td
import os
import pybol
from mdpow.equil import Simulation
from gromacs.utilities import in_dir

from mdpow.run import fep_simulation
from mdpow.config import get_configuration

class TestFEPScript(object):
    def setup(self):
        self.tmpdir = td.TempDir()
        self.old_path = os.getcwd()
        self.resources = os.path.join(
            self.old_path, 'mdpow', 'tests', 'testing_resources')
        self.m = pybol.Manifest(os.path.join(self.resources,'manifest.yml'))
        self.m.assemble('md_npt',self.tmpdir.name)

        S = Simulation(filename=os.path.join(
            self.tmpdir.name, 'benzene', 'water.simulation'))
        S.make_paths_relative(prefix=os.path.join(
            self.tmpdir.name,'benzene', 'Equilibrium', 'water'))
        S.dirs.includes = os.path.join(self.tmpdir.name, 'top')
        S.save()

    def teardown(self):
        self.tmpdir.dissolve()

    def _run_fep(self, solvent, dirname):
        cfg = get_configuration('runinput.yml')
        self.S = fep_simulation(cfg, solvent, dirname=dirname)

    def test_default_run(self):
        with in_dir(self.tmpdir.name, create=False):
            try:
                self._run_fep('water', 'benzene/')
            except:
                raise AssertionError('FEP simulations failed with exception:\n{0}'.format(str(err)))

            assert os.path.exists(os.path.join(self.tmpdir.name,
                                               'benzene', 'FEP', 'water', 'VDW', '0000', 'md.edr'))
