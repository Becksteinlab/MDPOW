from __future__ import absolute_import

from . import tempdir as td

import py.path

import pybol
import pytest

from ..analysis.ensemble import Ensemble, EnsembleAnalysis, EnsembleAtomGroup
from ..analysis.dihedral import DihedralAnalysis

from pkg_resources import resource_filename

RESOURCES = py.path.local(resource_filename(__name__, 'testing_resources'))
MANIFEST = RESOURCES.join("manifest.yml")

# https://docs.pytest.org/en/7.2.x/how-to/assert.html#assertraises

class TestRunExcept(object):

    def setup(self):
        self.tmpdir = td.TempDir()
        self.m = pybol.Manifest(str(RESOURCES / 'manifest.yml'))
        self.m.assemble('example_FEP', self.tmpdir.name)
        self.Ens = Ensemble(dirname=self.tmpdir.name, solvents=['water'])

    def teardown(self):
        self.tmpdir.dissolve()

    @pytest.mark.xfail(raises=NotImplementedError)
    def test_except(self):
        dh = self.Ens.select_atoms('name C4 or name C17 or name S2 or name N3')
        dh_run = DihedralAnalysis([dh]).run(start=0, stop=4, step=1)
        dh_run

    def test_single_universe(self):
        dh = self.Ens.select_atoms('name C4 or name C17 or name S2 or name N3')
        with pytest.raises(NotImplementedError):
            DihedralAnalysis([dh])._single_universe()
