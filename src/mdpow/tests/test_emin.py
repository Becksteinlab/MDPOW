import mdpow.equil
import tempdir as td
import os

import pkg_resources

TEST_RESOURCES = pkg_resources.resource_filename(
    __name__, 'testing_resources')

class TestEnergyMinimization(object):

    def setup(self):
        self.tmpdir = td.TempDir()
        self.resources = TEST_RESOURCES
        # TODO Instantiate Simulation object from existing files in resources

    def _run_emin(self):
        """To be used as a helper function in tests."""
        pass

    def test_new_struct_exists(self):
        """Tests if the new structure exists where it is expected to."""
        pass

    def test_new_struct_diff(self):
        """Tests if the new structure differs from the initial structure."""
        pass

