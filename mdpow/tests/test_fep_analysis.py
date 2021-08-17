import pytest

import os.path
import shutil
import copy

import mdpow.fep

from . import STATES

@pytest.fixture
def FEP_dir(tmpdir):
    name = STATES['FEP'].basename  # a py.path.local
    fepdir = tmpdir.join(name)     # a py.path.local
    shutil.copytree(STATES['FEP'].strpath, fepdir.strpath)
    assert os.path.isdir(fepdir.strpath)
    return fepdir

def setup_Ghyd(fepdir):
    basedir = fepdir.join("benzene")
    gsolv = basedir.join("FEP", "water", "Gsolv.fep")
    G = mdpow.fep.Ghyd(filename=gsolv.strpath)
    # patch paths
    G.basedir = basedir.strpath
    G.filename = gsolv.strpath
    return G


@pytest.fixture
def Ghyd(FEP_dir):
    return setup_Ghyd(FEP_dir)

@pytest.fixture
def Ghyd_other(FEP_dir):
    return setup_Ghyd(FEP_dir)

def test_load_Ghyd(Ghyd):
    assert isinstance(Ghyd, mdpow.fep.Ghyd)

@pytest.mark.parametrize("kwargs", (
    {},
    {'SI': True, 'estimator': 'alchemlyb', 'method': 'TI'},
    {'SI': False, 'estimator': 'alchemlyb', 'method': 'MBAR', 'force': False},
    {'SI': False, 'estimator': 'mdpow', 'method': 'TI', 'force': True},
    ),
    ids=["defaults",
         "SI=True, estimator='alchemlyb', method='TI'",
         "SI=False, estimator='alchemlyb', method='MBAR', force=False",
         "SI=False, estimator='mdpow', method='TI', force=True",
         ]
    )
def test_p_transfer(Ghyd, Ghyd_other, kwargs):
    """Test transfer water <-> water with same data."""
    G1 = Ghyd
    G2 = Ghyd_other

    transferFE, logPow = mdpow.fep.p_transfer(G1, G2, **kwargs)

    assert transferFE == pytest.approx(0.0)
    assert logPow == pytest.approx(0.0)

def test_p_transfer_wrong_method(Ghyd, Ghyd_other):
    """Test transfer water <-> water with same data."""
    G1 = Ghyd
    G2 = Ghyd_other

    with pytest.raises(ValueError,
                       match="Method MBAR is not implemented in MDPOW, use estimator='alchemlyb'"):
        mdpow.fep.p_transfer(G1, G2, estimator="mdpow", method="MBAR")


def test_pOW_error(Ghyd, Ghyd_other):
    with pytest.raises(ValueError):
        mdpow.fep.pOW(Ghyd, Ghyd_other)

def test_pCW_error(Ghyd, Ghyd_other):
    with pytest.raises(ValueError):
        mdpow.fep.pCW(Ghyd, Ghyd_other)
