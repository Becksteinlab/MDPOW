import tempdir as td
import os
import pybol
import numpy as np
from gromacs.utilities import in_dir

import mdpow.fep
import mdpow.equil

import pkg_resources

TEST_RESOURCES = pkg_resources.resource_filename(
    __name__, 'testing_resources')

class Test_Gsolv_manual(object):

    def setup(self):
        self.tmpdir = td.TempDir()
        self.resources = TEST_RESOURCES
        self.m = pybol.Manifest(os.path.join(self.resources,'manifest.yml'))
        self.m.assemble('md_npt',self.tmpdir.name)
        simulation_filename = os.path.join(self.tmpdir.name,'benzene',
                                  'water.simulation')
        self.S = mdpow.equil.Simulation(filename = simulation_filename)

        self.S.make_paths_relative(prefix=os.path.join(
           self.tmpdir.name,'benzene', 'Equilibrium', 'water'))
        self.S.dirs.includes = os.path.join(self.tmpdir.name, 'top')
        self.S.save()

    def teardown(self):
        self.tmpdir.dissolve()

    def _setup(self, **kwargs):
        with in_dir(self.tmpdir.name, create=False):
            self.Gsolv = mdpow.fep.Gsolv(simulation=self.S, molecule='BNZ',
                    **kwargs)
            self.Gsolv.setup(maxwarn=1)

    def test_default_setup(self):
        self._setup()

    def test_list_foreign_lambdas(self):
        lambda_coulomb = [0,0.5,1.0]
        lambda_vdw = [0,0.2,1.0]
        self._setup(lambda_coulomb=lambda_coulomb, lambda_vdw=lambda_vdw)

    def test_array_foreign_lambdas(self):
        lambda_coulomb = np.array([0,0.5,1.0])
        lambda_vdw = np.array([0,0.2,1.0])
        self._setup(lambda_coulomb=lambda_coulomb, lambda_vdw=lambda_vdw)
