import pytest
import pybol

import gromacs.run
import gromacs.exceptions

from numpy.testing import assert_equal

import mdpow.run
import mdpow.config


@pytest.fixture
def cfg():
    # default bundled config
    return mdpow.config.get_configuration()


@pytest.mark.parametrize(
    "protocols",
    [
        ["energy_minimize"],
        ["MD_relaxed"],
        ["MD_NPT", "FEP"],
    ],
)
def test_get_mdp_files(cfg, protocols):
    mdpfiles = mdpow.run.get_mdp_files(cfg, protocols)
    assert len(mdpfiles) == len(protocols)
    assert set(mdpfiles.keys()) == set(protocols)
    assert all([mdp.endswith(".mdp") for mdp in mdpfiles.values()])


@pytest.mark.parametrize("protocols", [["FEP"], ["Jabberwocky", "Mad Hatter"]])
def test_get_mdp_files_None(cfg, protocols):
    # modify cfg
    del cfg.conf["FEP"]["mdp"]
    with pytest.warns(mdpow.config.NoOptionWarning):
        mdpfiles = mdpow.run.get_mdp_files(cfg, ["FEP"])
    assert mdpfiles == {}


def test_get_mdp_files_ValueError(cfg):
    # modify cfg with a non-existant file
    cfg.conf["FEP"]["mdp"] = "smoke_and_mirror.mdp"
    with pytest.raises(ValueError):
        mdpow.run.get_mdp_files(cfg, ["MD_NPT", "FEP"])


# To test the failure modes of runMD_or_exit() we mock the functions
# and methods that would return failures, so that we don't have to
# actually run simulations.


@pytest.fixture
def MDrunner_failure(monkeypatch):
    # mock gromacs.run.MDrunner: pretend that the simulation failed
    def mock_run_check(*args, **kwargs):
        return False

    monkeypatch.setattr(gromacs.run.MDrunner, "run_check", mock_run_check)


# mock gromacs.run.check_mdrun_success(logfile)
@pytest.fixture
def check_mdrun_success_failure(monkeypatch):
    # pretend simulation has not completed as indicated by log file
    def mock_check_mdrun_success(arg):
        return False

    monkeypatch.setattr(gromacs.run, "check_mdrun_success", mock_check_mdrun_success)


@pytest.fixture
def check_mdrun_success_none(monkeypatch):
    # pretend no simulation has been run so there's no logfile to
    # check and check_mdrun_success() returns None
    def mock_check_mdrun_success(arg):
        return None

    monkeypatch.setattr(gromacs.run, "check_mdrun_success", mock_check_mdrun_success)


@pytest.mark.parametrize(
    "runlocal,exception",
    [
        (True, gromacs.exceptions.GromacsError),
        (False, gromacs.exceptions.MissingDataError),
    ],
)
def test_runMD_or_exit_exceptions(
    runlocal,
    exception,
    cfg,
    MDrunner_failure,
    check_mdrun_success_failure,
    monkeypatch,
    tmpdir,
):
    params = {"deffnm": "md"}
    S = {}

    def mock_getboolean(*args):
        return runlocal

    monkeypatch.setattr(cfg, "getboolean", mock_getboolean)

    with pytest.raises(exception):
        mdpow.run.runMD_or_exit(
            S, "FEP", params, cfg, dirname=str(tmpdir), exit_on_error=False
        )


def test_runMD_or_exit_None(cfg, check_mdrun_success_none, monkeypatch, tmpdir):
    # special case where runlocal=False and no simulation has been run
    # so there's no logfile to check and check_mdrun_success() returns
    # None
    params = {"deffnm": "md"}
    S = {}

    def mock_getboolean(*args):
        return False

    monkeypatch.setattr(cfg, "getboolean", mock_getboolean)

    return_value = mdpow.run.runMD_or_exit(
        S, "FEP", params, cfg, dirname=str(tmpdir), exit_on_error=False
    )
    assert return_value is None


@pytest.mark.parametrize("runlocal", [True, False])
def test_runMD_or_exit_SysExit(
    runlocal, cfg, MDrunner_failure, check_mdrun_success_failure, monkeypatch, tmpdir
):
    params = {"deffnm": "md"}
    S = {}

    def mock_getboolean(*args):
        return runlocal

    monkeypatch.setattr(cfg, "getboolean", mock_getboolean)

    with pytest.raises(SystemExit):
        mdpow.run.runMD_or_exit(
            S, "FEP", params, cfg, dirname=str(tmpdir), exit_on_error=True
        )
