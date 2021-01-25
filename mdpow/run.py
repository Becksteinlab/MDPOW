# equil.py
# Copyright (c) 2010-2011 Oliver Beckstein

"""
:mod:`mdpow.run` --- Performing complete simulation protocols
=============================================================

The module provides building blocks for complete simulation protocols
(or pipelines). Each protocol is written as a function that takes a
run input file and the solvent type as input.

:ref:`mdpow-scripts-label` make use of the building blocks.

Typically, *journalling* is enabled, i.e. the tasks remember which
stages have already been completed and can be restarted directly from
the last completed stage. (Restarts are only implemeneted at the level
of individual steps in a MDPOW protocol, not at the level of
continuing interrupted simulations using the Gromacs restart files.)

Input is read from the run input *cfg* file.

.. SeeAlso:: :mod:`mdpow.restart` for the journalling required for restarts.


Protocols
---------

.. autofunction:: equilibrium_simulation
.. autofunction:: fep_simulation

Support
-------

.. autofunction:: setupMD
.. autofunction:: get_mdp_files
.. autofunction:: runMD_or_exit

"""

from __future__ import absolute_import

import sys
import os
import errno

import gromacs.run

from .config import get_configuration, set_gromacsoutput
from . import equil
from . import fep
from .restart import checkpoint

import logging
logger = logging.getLogger('mdpow.run')


def setupMD(S, protocol, cfg):
    """setup MD simulation *protocol* using the runinput file *cfg*"""

    try:
        maxwarn = cfg.getint(protocol, "maxwarn")
    except KeyError:
        maxwarn = 0

    simulation_protocol = S.get_protocol(protocol)
    params = simulation_protocol(runtime=cfg.getfloat(protocol, "runtime"),
                                 qscript=cfg.getlist(protocol, "qscript"),
                                 maxwarn=maxwarn,
                                 )
    return params

def get_mdp_files(cfg, protocols):
    """Get file names of MDP files from *cfg* for all *protocols*"""
    import ConfigParser

    mdpfiles = {}
    for protocol in protocols:
        try:
            mdp = cfg.getpath(protocol, 'mdp')
        except KeyError:
            # skip anything for which we do not define sections, such as
            # the dummy run protocols
            mdp = None
        except ConfigParser.NoOptionError:
            # Should not happen... let's continue and wait for hard-coded defaults
            logger.error("No 'mdp' config file entry for protocol [%s]---check input files!", protocol)
            mdp = None
        except ValueError:
            # But it is a problem if we can't find a file!
            logger.critical("Failed to find custom MDP file %r for protocol [%s]",
                            cfg.get(protocol, 'mdp'), protocol)
            raise
        if mdp:
            mdpfiles[protocol] = mdp
            logger.debug("%(protocol)s: Using MDP file %(mdp)r from config file", vars())
    return mdpfiles

def runMD_or_exit(S, protocol, params, cfg, **kwargs):
    """run simulation

    Can launch :program:`mdrun` itself (:class:`gromacs.run.MDrunner`) or exit so
    that the user can run the simulation independently. Checks if the
    simulation has completed and sets the return value to ``True`` if this is
    the case.

    - Configuration parameters are taken from the section *protocol* of the
      runinput configuration *cfg*.

    - The directory in which the simulation input files reside can be provided
      as keyword argument *dirname* or taken from `S.dirs[protocol]`.

    - Other *kwargs* are interpreted as options for
      :class:`~gromacs.tools.Mdrun`.

    It never returns ``False`` but instead does a :func:`sys.exit`.
    """
    dirname = kwargs.pop("dirname", None)
    if dirname is None:
        try:
            dirname = S.dirs[protocol]
        except KeyError:
            raise ValueError("S.dirs does not have protocol %r" % protocol)
        except AttributeError:
            raise ValueError("supply dirname as a keyword argument")
    simulation_done = False
    if cfg.getboolean(protocol, "runlocal"):
        logger.info("Running %s (%s.log) ... stand by.", protocol, params['deffnm'])
        logger.info("Run directory: %(dirname)s", vars())
        mdrun = gromacs.run.MDrunner(
            dirname=dirname, deffnm=params['deffnm'],
            v=cfg.getboolean('mdrun','verbose'),
            stepout=cfg.getint('mdrun','stepout'),
            nice=cfg.getint('mdrun','nice'),
            nt=cfg.get('mdrun','maxthreads'),
            cpi=True)
        simulation_done = mdrun.run_check()
        if not simulation_done:
            # should probably stop
            logger.critical("Failed %(protocol)s, investigate manually.", vars())
            sys.exit(1)
    else:
        # must check if the simulation was run externally
        logfile = os.path.join(dirname, params['deffnm']+os.extsep+"log")
        logger.debug("Checking logfile %r if simulation has been completed.", logfile)
        simulation_done = gromacs.run.check_mdrun_success(logfile) ### broken??
        if simulation_done is None:
            logger.info("Now go and run %(protocol)s in directory %(dirname)r.", vars())
            sys.exit(0)
        elif simulation_done is False:
            logger.warn("Simulation %(protocol)s in directory %(dirname)r is incomplete (log=%)logfile)s).", vars())
            sys.exit(1)
        logger.info("Simulation %(protocol)s seems complete (log=%(logfile)s)", vars())
    return simulation_done


def equilibrium_simulation(cfg, solvent, **kwargs):
    """Set up and run equilibrium simulation.

    See tutorial for the individual steps.

    .. Note:: Depending on the settings in the run input file
       (``runlocal``), :program:`mdrun` is executed at various stages,
       and hence this process can take a while.
    """
    deffnm = kwargs.pop('deffnm', "md")
    Simulations = {
        'water': equil.WaterSimulation,
        'octanol': equil.OctanolSimulation,
        'wetoctanol': equil.WetOctanolSimulation,
        'cyclohexane':equil.CyclohexaneSimulation,
        }
    try:
        Simulation = Simulations[solvent]
    except KeyError:
        raise ValueError("solvent must be one of {0}".format(
            ", ".join(Simulations.keys())))

    # generate a canonical path under dirname
    topdir = kwargs.get("dirname", None)
    if topdir is None:
        topdir = cfg.get("setup", "name")
    dirname = os.path.join(topdir, Simulation.dirname_default)
    savefilename = os.path.join(topdir, "%(solvent)s.simulation" % vars())
    # output to screen or hidden?
    set_gromacsoutput(cfg)

    # start pipeline
    if os.path.exists(savefilename):
        S = Simulation(filename=savefilename)
    else:
        # custom mdp files
        mdpfiles = get_mdp_files(cfg, Simulation.protocols)
        try:
            distance = cfg.get('setup', 'distance')
        except KeyError:
            distance = None # if no distance is specified, None = default
        try:
            boxtype = cfg.get('setup', 'boxtype')
        except KeyError:
            boxtype = None # if no distance is specified, None = default


        solventmodel = None
        try:
            solventmodel = cfg.get('setup', 'solventmodel')
            logger.info("Selected solvent model: {0}".format(solventmodel))
        except KeyError:
            solventmodel = None
            logger.info("Using default water model")
        # for other solvents we currently only have a single OPLS-AA
        # parameterization included and hence there is no mechanism to
        # choose between different models.

        S = Simulation(molecule=cfg.get("setup", "molecule"),
                       forcefield=cfg.get("setup", "forcefield"),
                       dirname=dirname, deffnm=deffnm, mdp=mdpfiles,
                       distance=distance,
                       solventmodel=solventmodel)

    if S.journal.has_not_completed("energy_minimize"):
        maxwarn = cfg.getint("setup", "maxwarn") or None
        prm = cfg.get("setup", "prm") or None
        maxthreads = cfg.get('mdrun', 'maxthreads') or None
        S.topology(itp=cfg.getpath("setup", "itp"), prm=prm)
        S.solvate(struct=cfg.getpath("setup", "structure"),
                  bt=cfg.get("setup", "boxtype"),
                  maxwarn=maxwarn)
        S.energy_minimize(maxwarn=maxwarn, mdrun_args={'nt': maxthreads})
        checkpoint('energy_minize', S, savefilename)
    else:
        logger.info("Fast-forwarding: setup + energy_minimize done")

    if S.journal.has_not_completed("MD_relaxed"):
        params = setupMD(S, "MD_relaxed", cfg)
        checkpoint("MD_relaxed", S, savefilename)
    else:
        params = {'deffnm': deffnm}
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
       (``runlocal``), :program:`mdrun` is executed sequentially for
       all windows and hence this can take a long time. It is
       recommended to use ``runlocal = False`` in the run input file
       and submit all window simulations to a cluster.
    """

    deffnm = kwargs.pop('deffnm', "md")
    EquilSimulations = {
        'water': equil.WaterSimulation,
        'octanol': equil.OctanolSimulation,
        'wetoctanol': equil.WetOctanolSimulation,
        'cyclohexane':equil.CyclohexaneSimulation
        }
    Simulations = {
        'water': fep.Ghyd,
        'octanol': fep.Goct,
        'wetoctanol': fep.Gwoct,
        'cyclohexane':fep.Gcyclo
        }
    try:
        EquilSimulation = EquilSimulations[solvent]
        Simulation = Simulations[solvent]
    except KeyError:
        raise ValueError("solvent must be 'water', 'octanol', 'wetoctanol' or 'cyclohexane'")
    # generate a canonical path under dirname
    topdir = kwargs.get("dirname", None)
    if topdir is None:
        topdir = cfg.get("setup", "name")
    dirname = os.path.join(topdir, Simulation.dirname_default)
    # XXX nasty ... use same recipe to construct default save file name as in
    # Gsolv ... should be a static method or something else I can use before
    # the class is instantiated. Note that the pickle files live under dirname
    # and NOT topdir (bit of an historic inconsistency)
    savefilename = os.path.join(dirname, Simulation.__name__ + os.extsep + 'fep')
    # need pickle files for the equilibrium simulation ... another nasty guess:
    equil_savefilename = os.path.join(topdir, "%(solvent)s.simulation" % vars())
    try:
        equil_S = EquilSimulation(filename=equil_savefilename)
    except IOError, err:
        if err.errno == errno.ENOENT:
            logger.critical("Missing the equilibrium simulation %(equil_savefilename)r.", vars())
            logger.critical("Run `mdpow-equilibrium -S %s %s'  first!", solvent, "RUNINPUT.cfg")
        raise

    # output to screen or hidden?
    set_gromacsoutput(cfg)
    # start pipeline
    if os.path.exists(savefilename):
        S = Simulation(filename=savefilename, basedir=topdir)
        logger.info("Running FEP with saved settings from %r", savefilename)
        logger.info("Note: all configuration file options are ignored!")
    else:
                # method to be used "TI"/"BAR"
        method = cfg.get("FEP", "method")
        logger.info("Running FEP with the %r implementation in Gromacs", method)
        # custom mdp file
        mdp = cfg.findfile("FEP", "mdp")
        logger.debug("Using [FEP] MDP file %r from config file", mdp)

        # lambda schedules can be read from [FEP_schedule_*] sections
        schedules = {'coulomb': fep.FEPschedule.load(cfg, "FEP_schedule_Coulomb"),
                     'vdw': fep.FEPschedule.load(cfg, "FEP_schedule_VDW"),
                     }
        logger.debug("Loaded FEP schedules %r from config file", schedules.keys())

        # Note that we set basedir=topdir (and *not* dirname=dirname!)...FEP is a bit convoluted
        S = Simulation(simulation=equil_S, runtime=cfg.getfloat("FEP", "runtime"),
                       basedir=topdir, deffnm=deffnm, mdp=mdp, schedules=schedules,
                       method=method)
    if S.journal.has_not_completed("setup"):
        params = S.setup(qscript=cfg.getlist("FEP", "qscript"),
                         maxwarn=cfg.getint("FEP", "maxwarn"))
        checkpoint("setup", S, savefilename)
    else:
        params = {'deffnm': deffnm}
        logger.info("Fast-forwarding: FEP setup done")

    if S.journal.has_not_completed("fep_run"):
        def run_all_FEPS():
            for wdir in S.fep_dirs():
                runMD_or_exit(S, "FEP", params, cfg, dirname=wdir, dgdl=True)
        wrapper = S.get_protocol("fep_run")
        wrapper(run_all_FEPS)
        checkpoint("fep_run", S, savefilename)
    else:
        logger.info("Fast-forwarding: fep (run) done")

    logger.info("FEP simulation phase complete, use %(savefilename)r to continue.",
                vars())
    return savefilename
