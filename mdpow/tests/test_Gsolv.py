import tempdir as td
import os
import pybol
import numpy as np
import gromacs
from gromacs.utilities import in_dir
from mdpow.fep import *
from mdpow.equil import *


class Test_Gsolv_manual(object):

    def setup(self):
        self.tmpdir = td.TempDir()
        self.old_path = os.getcwd()
        self.resources = os.path.join(self.old_path, 'testing_resources')
        self.m = pybol.Manifest(os.path.join(self.resources, 'manifest.yml'))
        self.m.assemble('md_npt', self.tmpdir.name)
        simulation_filename = os.path.join(self.tmpdir.name, 'benzene', 'water.simulation')
        self.S = Simulation(filename=simulation_filename)

        self.S.make_paths_relative(prefix=os.path.join(
           self.tmpdir.name, 'benzene', 'Equilibrium', 'water'))
        self.S.dirs.includes = os.path.join(self.tmpdir.name, 'top')
        self.S.save()

    def teardown(self):
        self.tmpdir.dissolve()

    def _setup(self, **kwargs):
        with in_dir(self.tmpdir.name, create=False):
            self.Gsolv = Gsolv(simulation=self.S, molecule='BNZ',
                               mdp=os.path.join(self.old_path, '../templates', 'bar_opls.mdp'), **kwargs)
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
