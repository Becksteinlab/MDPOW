from . import tempdir as td

import os
import pybol
import numpy as np
from gromacs.utilities import in_dir

from mdpow import fep, equil, config

from . import RESOURCES

class Test_Gsolv_manual(object):

    def setup(self):
        self.tmpdir = td.TempDir()
        self.m = pybol.Manifest(str(RESOURCES / 'manifest.yml'))
        self.m.assemble('md_npt',self.tmpdir.name)
        simulation_filename = os.path.join(self.tmpdir.name,'benzene',
                                  'water.simulation')
        self.S = equil.Simulation(filename = simulation_filename)

        self.S.make_paths_relative(prefix=os.path.join(
           self.tmpdir.name,'benzene', 'Equilibrium', 'water'))
        self.S.dirs.includes = os.path.join(self.tmpdir.name, 'top')
        self.S.save()

    def teardown(self):
        self.tmpdir.dissolve()

    def _setup(self, **kwargs):
        mdp = config.get_template("bar_opls.mdp")
        with in_dir(self.tmpdir.name, create=False):
            self.Gsolv = fep.Gsolv(simulation=self.S, molecule='BNZ',
                                   mdp=mdp ,**kwargs)
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
