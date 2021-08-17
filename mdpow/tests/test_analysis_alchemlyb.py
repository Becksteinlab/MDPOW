import os.path

import pytest

import yaml
import pybol

from numpy.testing import assert_array_almost_equal
import pandas

import pickle

import mdpow.fep
from . import RESOURCES

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
    manifest = yaml.safe_load(MANIFEST.open())
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
        # Loading only works with magic code in restart.Journalled.load()
        # so that Pickles produced with Python 2 can also be read with 3:
        G = mdpow.fep.Ghyd(filename=gsolv.strpath)
        # patch paths
        G.basedir = pth.strpath
        G.filename = gsolv.strpath
        return G

    @pytest.mark.parametrize('method, Gibbs, coulomb, vdw', [
                ('TI',
                 (-3.901068,  0.550272),
                 (8.417035, 0.22289),
                 (-4.515967,  0.50311)),
                ('BAR',
                 (-4.091241, 0.385413),
                 (8.339705, 0.166802),
                 (-4.248463, 0.347449)),
                ('MBAR',
                 (-6.793117,  0.475149),
                 (8.241836, 0.219235),
                 (-1.448719,  0.421548))
                ])

    @pytest.mark.xfail(pandas.__version__.startswith("1.3.0"),
                       reason="bug in pandas 1.3.0 see alchemistry/alchemlyb#147")
    def test_estimator_alchemlyb(self, fep_benzene_directory, method,
                                 Gibbs, coulomb, vdw):
        G = self.get_Gsolv(fep_benzene_directory)
        G.method = method
        G.start = 0
        G.stop = None
        # ensure conversion EDR to XVG.bz2; if the fixture is session scoped
        # then other workers will pick up these files. Make sure that only one
        # runs convert because there is no file locking, if in doubt, make
        # fep_benzene_directory locally scoped
        G.convert_edr()
        try:
            G.analyze_alchemlyb(force=True, autosave=False, SI=False)
        except IOError as err:
            raise AssertionError("Failed to convert edr to xvg: {0}: {1}".format(
                err.strerror, err.filename))
        DeltaA = G.results.DeltaA
        assert_array_almost_equal(DeltaA.Gibbs.astuple(), Gibbs,
                                  decimal=5)  # with more recent versions of pandas/alchemlyb/numpy the original values are only reproduced to 5 decimals, see PR #166")
        assert_array_almost_equal(DeltaA.coulomb.astuple(), coulomb,
                                  decimal=5)  # with more recent versions of pandas/alchemlyb/numpy the original values are only reproduced to 5 decimals, see PR #166")
        assert_array_almost_equal(DeltaA.vdw.astuple(), vdw,
                                  decimal=5)  # with more recent versions of pandas/alchemlyb/numpy the original values are only reproduced to 5 decimals, see PR #166")

    @pytest.mark.xfail(pandas.__version__.startswith("1.3.0"),
                       reason="bug in pandas 1.3.0 see alchemistry/alchemlyb#147")
    def test_SI(self, fep_benzene_directory):
        G = self.get_Gsolv(fep_benzene_directory)
        G.method = 'TI'
        G.start = 0
        G.stop = None
        G.convert_edr()
        try:
            G.analyze_alchemlyb(force=True, SI=True, autosave=False)
        except IOError as err:
            raise AssertionError("Failed to convert edr to xvg: {0}: {1}".format(
                err.strerror, err.filename))
        DeltaA = G.results.DeltaA
        assert_array_almost_equal(DeltaA.Gibbs.astuple(), (-2.908885,  2.175976),
                                  decimal=5)  # with more recent versions of pandas/alchemlyb/numpy the original values are only reproduced to 5 decimals, see PR #166")
        assert_array_almost_equal(DeltaA.coulomb.astuple(), (7.755779, 0.531481),
                                  decimal=5)  # with more recent versions of pandas/alchemlyb/numpy the original values are only reproduced to 5 decimals, see PR #166")
        assert_array_almost_equal(DeltaA.vdw.astuple(), (-4.846894,  2.110071),
                                  decimal=5)  # with more recent versions of pandas/alchemlyb/numpy the original values are only reproduced to 5 decimals, see PR #166")

    @pytest.mark.xfail(pandas.__version__.startswith("1.3.0"),
                       reason="bug in pandas 1.3.0 see alchemistry/alchemlyb#147")
    def test_start_stop_stride(self, fep_benzene_directory):
        G = self.get_Gsolv(fep_benzene_directory)
        G.method = 'TI'
        G.start = 10
        G.stride = 2
        G.stop = 200
        G.convert_edr()
        try:
            G.analyze_alchemlyb(force=True, autosave=False, SI=False)
        except IOError as err:
            raise AssertionError("Failed to convert edr to xvg: {0}: {1}".format(
                err.strerror, err.filename))
        DeltaA = G.results.DeltaA
        assert_array_almost_equal(DeltaA.Gibbs.astuple(), (-3.318109,  0.905128),
                                  decimal=5)  # with more recent versions of pandas/alchemlyb/numpy the original values are only reproduced to 5 decimals, see PR #166")
        assert_array_almost_equal(DeltaA.coulomb.astuple(), (8.146806, 0.348866),
                                  decimal=5)  # with more recent versions of pandas/alchemlyb/numpy the original values are only reproduced to 5 decimals, see PR #166")
        assert_array_almost_equal(DeltaA.vdw.astuple(), (-4.828696,  0.835195),
                                  decimal=5)  # with more recent versions of pandas/alchemlyb/numpy the original values are only reproduced to 5 decimals, see PR #166")
