import pytest

from numpy.testing import assert_equal

import mdpow.run
import mdpow.config

@pytest.fixture
def cfg():
    # default bundled config
    return mdpow.config.get_configuration()

@pytest.mark.parametrize("protocols", [["energy_minimize"],
                                       ["MD_relaxed"],
                                       ["MD_NPT", "FEP"],
                                      ])
def test_get_mdp_files(cfg, protocols):
    mdpfiles = mdpow.run.get_mdp_files(cfg, protocols)
    assert len(mdpfiles) == len(protocols)
    assert set(mdpfiles.keys()) == set(protocols)
    assert all([mdp.endswith(".mdp") for mdp in mdpfiles.values()])

@pytest.mark.parametrize("protocols", [["FEP"],
                                       ["Jabberwocky", "Mad Hatter"]
                                      ])
def test_get_mdp_files_None(cfg, protocols):
    # modify cfg
    del cfg.conf['FEP']['mdp']
    with pytest.warns(mdpow.config.NoOptionWarning):
        mdpfiles = mdpow.run.get_mdp_files(cfg, ["FEP"])
    assert mdpfiles == {}

def test_get_mdp_files_ValueError(cfg):
    # modify cfg with a non-existant file
    cfg.conf['FEP']['mdp'] = "smoke_and_mirror.mdp"
    with pytest.raises(ValueError):
        mdpow.run.get_mdp_files(cfg, ["MD_NPT", "FEP"])
