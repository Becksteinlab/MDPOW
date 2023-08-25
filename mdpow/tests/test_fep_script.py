import os
import pytest
import pybol

from . import tempdir as td

import gromacs

import mdpow
from mdpow.equil import Simulation
from mdpow.run import fep_simulation
from mdpow.config import get_configuration

from . import RESOURCES


class TestFEPScript(object):
    def setup_method(self):
        self.tmpdir = td.TempDir()
        self.old_path = os.getcwd()
        self.resources = RESOURCES
        self.m = pybol.Manifest(str(self.resources / "manifest.yml"))
        self.m.assemble("md_npt", self.tmpdir.name)

        S = Simulation(
            filename=os.path.join(self.tmpdir.name, "benzene", "water.simulation")
        )
        S.make_paths_relative(
            prefix=os.path.join(self.tmpdir.name, "benzene", "Equilibrium", "water")
        )
        S.dirs.includes = os.path.join(self.tmpdir.name, "top")
        S.save()

    def teardown_method(self):
        self.tmpdir.dissolve()

    def _run_fep(self, solvent, dirname):
        cfg = get_configuration("runinput.yml")
        if gromacs.release.startswith("4"):
            # For GROMACS 4.6.5 explicitly enable the group neighbor
            # scheme by creating a copy of the MDP file in the current
            # directory with MDP cutoff-scheme option changed. The local
            # MDP file will be picked up in preference to the default one
            # in the templates.
            fep_mdp_name = cfg.get("FEP", "mdp")
            mdp = mdpow.config.get_template(fep_mdp_name)
            gromacs.cbook.edit_mdp(
                mdp,
                new_mdp=os.path.join(os.getcwd(), fep_mdp_name),
                cutoff_scheme="group",
            )
        self.savefilename = fep_simulation(
            cfg, solvent, dirname=dirname, exit_on_error=False
        )

    def test_default_run(self, capsys):
        with gromacs.utilities.in_dir(self.tmpdir.name, create=False):
            try:
                self._run_fep("water", "benzene/")
            except Exception as err:
                # check if log file contains
                #   'There are 2 perturbed non-bonded pair interactions beyond the pair-list cutoff'
                # This error can happen and still counts as passing. (It is stochastic and is due
                # to a poorly designed test system together with GROMACS having become better at detecting
                # internal problems.)
                captured = capsys.readouterr()
                for line in captured:
                    if (
                        "perturbed non-bonded pair interactions beyond the pair-list cutoff'"
                        in line
                    ):
                        pytest.xfail(
                            "Stochastic test failure (perturbed non-bonded beyond cutoff). "
                            "Still works as expected, see #175."
                        )
                        break
                else:
                    raise AssertionError(
                        "FEP simulations failed with exception:\n{0}".format(str(err))
                    )

            assert os.path.exists(
                os.path.join(
                    self.tmpdir.name, "benzene", "FEP", "water", "VDW", "0000", "md.edr"
                )
            )
