# equil.py
# Copyright (c) 2010-2011 Oliver Beckstein

"""
:mod:`mdpow.run` --- Performing complete simulation protocols
=============================================================

The module provides building blocks for complete simulation protocols
(or pipelines). Each protocol is written as a function that takes a
run input file and the solvent type as input.

:ref:`mdpow-scripts-label` make use of the building blocks.

Protocols
---------

.. autofunction:: equilibrium_simulation
.. autofunction:: fep_simulation

Support
-------

.. autoclass:: MDrunnerSimple

"""

import sys
import os

import gromacs.run

import mdpow.equil
import mdpow.fep
from mdpow.config import get_configuration
from mdpow.restart import checkpoint

import logging
logger = logging.getLogger('mdpow.run')


class MDrunnerSimple(gromacs.run.MDrunner):
    """Gromacs mdrun.

    Use Gromacs 4.5.x with threaded mdrun to get max performance and
    just leave everything as it is. For older version you can change
    the name and add a mpiexec binary to launch a multiprocessor
    job. See :mod:`gromacs.run` for details.

    """
    mdrun = "mdrun"
    mpiexec = None

def setupMD(S, protocol, cfg):
    """setup MD simulation *protocol* using the runinput file *cfg*"""

    simulation_protocol = S.get_protocol(protocol)
    params = simulation_protocol(runtime=cfg.getfloat(protocol, "runtime"),
                                 qscript=cfg.getlist(protocol, "qscript"))
    return params

def runMD_or_exit(S, protocol, params, cfg):
    """run simulation

    Can launch mdrun itself or exit so that the user can run the
    simulation independently. Checks if the simulation has
    completed and sets the return value to ``True`` if this is the
    case.

    It never returns ``False`` but instead does a sys.exit.
    """
    simulation_done = False
    if cfg.getboolean(protocol, "runlocal"):
        logger.info("Running %s (%s.log) ... stand by.", protocol, params['deffnm'])
        logger.info("Run directory: %s", S.dirs[protocol])
        mdrun = MDrunnerSimple(dirname=S.dirs[protocol], deffnm=params['deffnm'],
                               v=cfg.getboolean('mdrun','verbose'),
                               stepout=cfg.getint('mdrun','stepout'),
                               nice=cfg.getint('mdrun','nice'),
                               nt=cfg.get('mdrun','maxthreads'),
                               cpi=True, append=True)
        simulation_done = mdrun.run_check()
        if not simulation_done:
            # should probably stop
            logger.critical("Failed %(protocol)s, investigate manually.", vars())
            sys.exit(1)
    else:
        # must check if the simulation was run externally
        logfile = os.path.join(S.dirs[protocol], params['deffnm']+os.extsep+"log")
        simulation_done = gromacs.run.check_mdrun_success(logfile)
        if simulation_done is None:
            logger.info("Now go and run %s in directory %r.", protocol, S.dirs[protocol])
            sys.exit(0)
        elif simulation_done is False:
            logger.warn("Simulation %s in directory %r is incomplete (log=%s).",
                        protocol, S.dirs[protocol], logfile)
            sys.exit(1)
        logger.info("Simulation %s seems complete (log=%s)", protocol, logfile)
    return simulation_done

def equilibrium_simulation(cfg, solvent, **kwargs):
    """Set up and run equilibrium simulation.

    See tutorial for the individual steps.

    .. Note:: Depending on the settings in the run input file
       (``runlocal``), :program:`mdrun` is executed at various stages,
       and hence this process can take a while.
    """
    Simulations = {
        'water': mdpow.equil.WaterSimulation,
        'octanol': mdpow.equil.OctanolSimulation,
        }
    try:
        Simulation = Simulations[solvent]
    except KeyError:
        raise ValueError("solvent must be 'water' or 'octanol'")

    # generate a canonical path under dirname
    topdir = kwargs.get("dirname", None)
    if topdir is None:
        topdir = cfg.get("setup", "name")
    dirname = os.path.join(topdir, Simulation.dirname_default)
    savefilename = os.path.join(topdir, "%(solvent)s.simulation" % vars())

    # start pipeline
    if os.path.exists(savefilename):
        S = Simulation(filename=savefilename)
    else:
        S = Simulation(molecule=cfg.get("setup", "molecule"), dirname=dirname)

    if S.journal.has_not_completed("energy_minimize"):
        maxwarn = cfg.getint("setup", "maxwarn") or None
        S.topology(itp=cfg.getpath("setup", "itp"))
        S.solvate(struct=cfg.getpath("setup", "structure"), maxwarn=maxwarn)
        S.energy_minimize(maxwarn=maxwarn)
        checkpoint('energy_minize', S, savefilename)
    else:
        logger.info("Fast-forwarding: setup + energy_minimize done")

    if S.journal.has_not_completed("MD_relaxed"):
        params = setupMD(S, "MD_relaxed", cfg)
        checkpoint("MD_relaxed", S, savefilename)
    else:
        # need to fudge params for the moment!!! (deffnm = md hard coded in mdpow...)
        params = {'deffnm': 'md'}
        logger.info("Fast-forwarding: MD_relaxed (setup) done")

    if S.journal.has_not_completed("MD_relaxed_run"):
        wrapper = S.get_protocol("MD_relaxed_run")
        success = wrapper(runMD_or_exit, S, "MD_relaxed", params, cfg)  # note: MD_relaxed!
        checkpoint("MD_relaxed_run", S, savefilename)
    else:
        logger.info("Fast-forwarding: MD_relaxed (run) done")

    if S.journal.has_not_completed("MD_NPT"):
        params = setupMD(S, "MD_NPT", cfg)
        checkpoint("MD_NPT", S, savefilename)
    else:
        logger.info("Fast-forwarding: MD_NPT (setup) done")

    if S.journal.has_not_completed("MD_NPT_run"):
        wrapper = S.get_protocol("MD_NPT_run")
        success = wrapper(runMD_or_exit, S, "MD_NPT", params, cfg)   # note: MD_NPT
        checkpoint("MD_NPT_run", S, savefilename)
    else:
        logger.info("Fast-forwarding: MD_NPT (run) done")

    logger.info("Equilibrium simulation phase complete, use %(savefilename)r to continue.",
                vars())
    return savefilename


def fep_simulation(cfg, solvent, **kwargs):
    """Set up and run FEP simulation.

    See tutorial for the individual steps.

    .. Note:: Depending on the settings in the run input file
       (``runlocal``), :program:`mdrun` is executed sequemtially for
       all windows and hence this can take a long time. It is
       recommended to use ``runocal = False`` in the run input file
       and submit all window simulations to a cluster.
    """
    Simulations = {
        'water': mdpow.fep.Ghyd,
        'octanol': mdpow.fep.Goct,
        }
    try:
        Simulation = Simulations[solvent]
    except KeyError:
        raise ValueError("solvent must be 'water' or 'octanol'")

    EquilSimulations = {
        'water': mdpow.equil.WaterSimulation,
        'octanol': mdpow.equil.OctanolSimulation,
        }
    try:
        EquilSimulation = Simulations[solvent]
    except KeyError:
        raise ValueError("solvent must be 'water' or 'octanol'")


    # generate a canonical path under dirname
    topdir = kwargs.get("dirname", None)
    if topdir is None:
        topdir = cfg.get("setup", "name")
    dirname = os.path.join(topdir, Simulation.dirname_default)
    # XXX nasty ... use same recipe to construct default save fle name as in
    # Gsolv ... should be a static method or something else I can use before
    # the class is instantiated. Note that the pickle files live under dirname
    # and NOT topdir (bit of an historic inconsistency)
    savefilename = os.path.join(dirname, Simulation.__name__ + os.extsep + 'fep')

    # need pickle files for the equilibrium simulation ... another nasty guess:
    equil_savefilename = os.path.join(topdir, "%(solvent)s.simulation" % vars())
    equil_S = EquilSimulation(filename=equil_savefilename)


    # start pipeline
    if os.path.exists(savefilename):
        S = Simulation(filename=savefilename)
    else:
        S = Simulation(simulation=equil_S, runtime=cfg.getfloat("FEP", "runtime"),
                       dirname=dirname)

#####
    if S.journal.has_not_completed("energy_minimize"):
        S.topology(itp=cfg.getpath("setup", "itp"))
        S.solvate(struct=cfg.getpath("setup", "structure"))
        S.energy_minimize()
        checkpoint('energy_minize', S, savefilename)
    else:
        logger.info("Fast-forwarding: setup + energy_minimze done")

    if S.journal.has_not_completed("MD_relaxed"):
        params = setupMD(S, "MD_relaxed", cfg)
        checkpoint("MD_relaxed", S, savefilename)
    else:
        # need to fudge params for the moment!!! (deffnm = md hard coded in mdpow...)
        params = {'deffnm': 'md'}
        logger.info("Fast-forwarding: MD_relaxed (setup) done")

    if S.journal.has_not_completed("MD_relaxed_run"):
        wrapper = S.get_protocol("MD_relaxed_run")
        success = wrapper(runMD_or_exit, S, "MD_relaxed", params, cfg)  # note: MD_relaxed!
        checkpoint("MD_relaxed_run", S, savefilename)
    else:
        logger.info("Fast-forwarding: MD_relaxed (run) done")

    if S.journal.has_not_completed("MD_NPT"):
        params = setupMD(S, "MD_NPT", cfg)
        checkpoint("MD_NPT", S, savefilename)
    else:
        logger.info("Fast-forwarding: MD_NPT (setup) done")

    if S.journal.has_not_completed("MD_NPT_run"):
        wrapper = S.get_protocol("MD_NPT_run")
        success = wrapper(runMD_or_exit, S, "MD_NPT", params, cfg)   # note: MD_NPT
        checkpoint("MD_NPT_run", S, savefilename)
    else:
        logger.info("Fast-forwarding: MD_NPT (run) done")

    logger.info("Equilibrium simulation phase complete, use %(savefilename)r to continue.",
                vars())
    return savefilename


