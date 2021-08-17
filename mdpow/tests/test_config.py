import pytest
import os.path
from io import StringIO

import numpy as np
from numpy.testing import assert_almost_equal, assert_equal

import mdpow.config

@pytest.fixture
def cfg():
    # default bundled config
    return mdpow.config.get_configuration()

@pytest.fixture
def minicfg():
    s = StringIO("setup:\n  name: 'Alice'\n  gromacsoutput: False\n")
    yield s
    s.close()

class TestConfigurationParser:
    def test_get_NoSectionError(self):
        cfg = mdpow.config.POWConfigParser()
        with pytest.raises(mdpow.config.NoSectionError,
                           match="Config file has no section Jabberwocky"):
            cfg.get('Jabberwocky', 'elump')

    def test_get_NoOptionWarning_gives_None(self, cfg):
        with pytest.warns(mdpow.config.NoOptionWarning,
                           match="Config file section FEP contains no "
                                 "option Jabberwocky. Using 'None'."):
            item = cfg.get("FEP", "Jabberwocky")
        assert item is None

    def test_get(self, cfg):
        item = cfg.get("setup","solventmodel")
        assert isinstance(item, str)
        assert item == "tip4p"

    def test_get_None(self, cfg):
        item = cfg.get("setup", "prm")
        assert item is None

    def test_getstr(self, cfg):
        args = "setup","solventmodel"
        item = cfg.getstr(*args)
        assert isinstance(item, str)
        assert item == "tip4p"

        assert item == cfg.get(*args)

    def test_getboolean(self, cfg):
        args = "setup", "gromacsoutput"
        item = cfg.getboolean(*args)
        assert isinstance(item, bool)
        assert item is True

        assert item == cfg.get(*args)

    def test_getint(self, cfg):
        args = "setup", "maxwarn"
        item = cfg.getboolean(*args)
        assert isinstance(item, int)
        assert item is 0

        assert item == cfg.get(*args)

    def test_getfloat(self, cfg):
        args = "MD_relaxed", "runtime"
        item = cfg.getboolean(*args)
        assert isinstance(item, float)
        assert item == pytest.approx(5.0)

        assert item == cfg.get(*args)

    def test_getpath(self, cfg):
        pth = "~/mirrors/jbwck.itp"
        cfg.conf['setup']['itp'] = pth

        item = cfg.getpath("setup", "itp")

        assert item == os.path.expanduser(pth)

    def test_getpath_None(self, cfg):
        # None gives None
        assert cfg.getpath("setup", "prm") is None

    def test_findfile(self, cfg):
        pth = cfg.findfile("FEP", "mdp")

        assert pth.endswith(".mdp")
        assert os.path.exists(pth)

    def test_findfile_None(self, cfg):
        # None gives None
        assert cfg.findfile("setup", "itp") is None

    def test_getlist(self, cfg):
        item = cfg.getlist("FEP_schedule_Coulomb", "lambdas")

        assert isinstance(item, list)
        assert item == ['0', '0.25', '0.5', '0.75', '1.0']

    def test_getlist_empty(self, cfg):
        # get an option with None for this test
        assert cfg.getlist("setup", "name") == []

    def test_getarray(self, cfg):
        item = cfg.getarray("FEP_schedule_Coulomb", "lambdas")

        assert isinstance(item, np.ndarray)
        assert_almost_equal(item, [0, 0.25, 0.5, 0.75, 1.0])

    def test_getarray_empty(self, cfg):
        # get an option with None for this test
        item = cfg.getarray("setup", "name")
        assert isinstance(item, np.ndarray)
        assert_equal(item, np.array([]))

    def test_write(self, cfg, tmp_path):
        pth = tmp_path / "new.yaml"
        cfg.write(pth)

        new_cfg = mdpow.config.POWConfigParser().readfp(pth.open())
        # compare dicts
        assert new_cfg.conf == cfg.conf

    def test_readfp_stream(self, minicfg):
        cfg = mdpow.config.POWConfigParser().readfp(minicfg)
        assert cfg.get("setup", "name") == "Alice"
        assert cfg.get("setup", "gromacsoutput") == False

    def test_merge(self, cfg, minicfg):
        cfg.merge(minicfg)
        # minicfg changes
        assert cfg.get("setup", "name") == "Alice"
        assert cfg.get("setup", "gromacsoutput") == False
        # original
        assert cfg.get("setup", "forcefield") == "OPLS-AA"
        assert cfg.get("FEP", "method") == "BAR"
