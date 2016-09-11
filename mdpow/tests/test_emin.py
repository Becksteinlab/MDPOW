import mdpow.equil
import tempdir as td
import os

class TestEnergyMinimization(object):
    
    def setup(self):
        self.tmpdir = td.TempDir()
        self.old_path = os.getcwd()
        self.resources = os.path.join(self.old_path, 'mdpow', 'tests', 'testing_resources')
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

