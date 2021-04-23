import os.path
import pytest
import py.path

import yaml
import pybol

from numpy.testing import assert_array_almost_equal

from six.moves import cPickle as pickle

import mdpow.fep

from pkg_resources import resource_filename
RESOURCES = py.path.local(resource_filename(__name__, 'testing_resources'))
MANIFEST = RESOURCES.join("manifest.yml")

def fix_manifest(topdir):
    """Create a temporary manifest with a custom `path`.

    Fix manifest in a local temporary copy in existing dir topdir
    where the `path` is an absolute path to our "states" directory. We
    use `pkg_resources.resource_filename` to anchor the path.

    Arguments
    ---------
    topdir : py.path.local
        existing temporary directory (as provided by, for instance,
        `pytest.tmpdir`)

    Returns
    -------
    new_manifest : py.path.local
        Path to the new manifest.

    Example
    -------
    Use as ::

        new_manifest = fix_manifest(tmpdir)
        m = pybol.Manifest(new_manifest.strpath)

    """
    manifest = yaml.load(MANIFEST.open())
    # simple heuristic: last element of the recorded manifest::path is the name
    # of the states directory, typically 'states' (from .../testing_resources/states)
    manifest['path'] = RESOURCES.join(os.path.basename(manifest['path'])).strpath
    new_manifest = topdir.join("local_manifest.yml")
    yaml.dump(manifest, stream=new_manifest.open("w"))
    return new_manifest


# session scope if read-only use

@pytest.fixture(scope="function")
def fep_benzene_directory(tmpdir_factory):
    topdir = tmpdir_factory.mktemp('analysis')
    m = pybol.Manifest(fix_manifest(topdir).strpath)
    m.assemble('FEP', topdir.strpath)
    return topdir.join("benzene")

class TestAnalyze(object):
    def get_Gsolv(self, pth):
        gsolv = pth.join("FEP", "water", "Gsolv.fep")
        G = pickle.load(gsolv.open())
        # patch paths
        G.basedir = pth.strpath
        G.filename = gsolv.strpath
        return G

    @staticmethod
    def assert_DeltaA(G):
        DeltaA = G.results.DeltaA
        assert_array_almost_equal(DeltaA.Gibbs.astuple(),
                                  (-3.7217472974883794, 2.3144288928034911),
                                  decimal=6)
        assert_array_almost_equal(DeltaA.coulomb.astuple(),
                                  (8.3346255170099575, 0.73620918517131495),
                                  decimal=6)
        assert_array_almost_equal(DeltaA.vdw.astuple(),
                                  (-4.6128782195215781, 2.1942144688960972),
                                  decimal=6)


    def test_convert_edr(self, fep_benzene_directory):
        G = self.get_Gsolv(fep_benzene_directory)
        try:
            G.analyze(force=True, autosave=False)
        except IOError as err:
            raise AssertionError("Failed to auto-convert edr to xvg: {0}: {1}".format(
                err.strerror, err.filename))
        self.assert_DeltaA(G)


    def test_TI(self, fep_benzene_directory):
        G = self.get_Gsolv(fep_benzene_directory)
        # ensure conversion EDR to XVG.bz2; if the fixture is session scoped
        # then other workers will pick up these files. Make sure that only one
        # runs convert because there is no file locking, if in doubt, make
        # fep_benzene_directory locally scoped
        G.convert_edr()
        try:
            G.analyze(force=True, autosave=False)
        except IOError as err:
            raise AssertionError("Failed to convert edr to xvg: {0}: {1}".format(
                err.strerror, err.filename))
        self.assert_DeltaA(G)

