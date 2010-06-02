# POW package __init__.py
# Copyright (c) 2010 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""
=====================================================================
:mod:`mdpow` --- Computing the octanol/water partitioning coefficient
=====================================================================

The :mod:`mdpow` module helps in setting up and analyzing absolute
free energy calculations of small molecules by molecular dynamics (MD)
simulations. By computing the hydration free energy and the solvation
free energy in octanol one can compute the octanol/water partitioning
coefficient, an important quantity that is used to characterize
drug-like compounds.

The MD simulations are performed with Gromacs_ 4.x

.. _Gromacs: http://www.gromacs.org


How to use the module
=====================

Before you can start you will need

  -  a coordinate file for the small molecule
  -  a Gromacs OPLS/AA topology (itp) file
  -  an installation of Gromacs_ 4.0.x.

Basic work flow
---------------

You will typically calculate two solvation free energies (free
energy of transfer of the solute from the liquid into the vacuum
phase):

  1. solvent = *water*

     1. set up a short equilibrium simulation of the molecule in a *water*
        box (and run the MD simulation);

     2. set up a free energy perturbation calculation of the ligand in
        water , which will yield the hydration free energy;

  2. solvent = *octanol*

     1. set up a short equilibrium simulation of the molecule in a *octanol*
        box (and run the MD simulation);

     2. set up a free energy perturbation calculation of the ligand in
        octanol , which will yield the solvation free energy in octanol;

  3. run these simulations on a cluster;

  4. analyze the output and combine the free energies to arrive at an
     estimate of the octanol-water partition coefficient.


Customized submission scripts for queuing systems
-------------------------------------------------

One can also generate run scripts for various queuing systems; check
the documentation for :mod:`gromacs.qsub` and in particular the
section on `writing queuing system templates`_ . You will have to

- add a template script to your private GromacsWrapper template
  directory (``~/.gromacswrapper/qscripts``); in this example we call
  it ``my_script.sge``;

- add the keyword *qscript* to the :meth:`mdpow.equil.Simulation.MD`
  and :meth:`mdpow.fep.Gsolv.setup` invocations; e.g. as

    *qscript* = ['my_script.sge', 'local.sh']

- submit the generated queuing system script to your queuing system, e.g. ::

    cd Equilibrium/water
    qsub my_script.sh

.. _writing queuing system templates: 
   http://sbcb.bioch.ox.ac.uk/oliver/software/GromacsWrapper/html/gromacs/blocks/qsub.html#queuing-system-templates

    

Example session: 1-octanol as a solute
--------------------------------------

In the following interactive python session we use octanol as an
example for a solute; all files are present in the package so one can
work through the example immediately.

Before starting :program:`python` (preferrably ipython_) make sure that the
Gromacs_ 4.0.x tools can be found, e.g. ``which grompp`` should show you the
path to :program:`grompp`.

.. _ipython: http://ipython.scipy.org
.. _Gromacs: http://www.gromacs.org


Water
~~~~~

Equilibrium simulation
......................

Make a directory ``octanol`` and copy the ``octanol.itp`` and
``octanol.gro`` file into it. Launch :program:`ipython` from this
directory and type::

   import mdpow.equil
   S = mdpow.equil.WaterSimulation(molecule="OcOH")
   S.topology(itp="octanol.itp")
   S.solvate(struct="octanol.gro")
   S.energy_minimize()
   S.MD_relaxed()
   # run the simulation in the MD_relaxed/ directory
   S.MD(runtime=50, qscript=['my_script.sge', 'local.sh'])  # only run for 50 ps in this tutorial
   S.save("water.simulation")                               # save setup for later (analysis stage)

Background (:kbd:`Ctrl-Z`) or quit (:kbd:`Ctrl-D`) python and run the
simulations in the ``MD_relaxed`` and ``MD_NPT`` subdirectory. You can modify
the ``local.sh`` script to your ends or use *qscript* to generate queuing system
scripts.

.. note:: Here we only run 50 ps equilibrium MD for testing. For production this
          should be substantially longer, maybe even 50 ns if you want to
          extract thermodynamic data.


Hydration free energy
.....................

Reopen the python session and set up a :class:`~mdpow.fep.Ghyd` object::
   
   import mdpow.fep
   gwat = mdpow.fep.Ghyd(molecule="OcOH", top="Equilibrium/water/top/system.top", struct="Equilibrium/water/MD_NPT/md.pdb", ndx="Equilibrium/water/solvation/main.ndx", runtime=100)
   
Alternatively, one can save some typing if we continue the last session and use
the :class:`mdpow.equil.Simulation` object (which we can re-load from its saved
state file from disk)::

  import mdpow.equil
  S = mdpow.equil.WaterSimulation(filename="water.simulation")  # only needed when quit
  gwat = mdpow.fep.Ghyd(simulation=S, runtime=100)

This generates all the input files under ``FEP/water``.

.. note:: Here we only run 100 ps per window for testing. For production this
          should be rather something like 5-10 ns (the default is 5 ns).

Then set up all input files::

  gwat.setup(qscript=['my_script.sge', 'local.sh'])

(The details of the FEP runs can be customized by setting some keywords (such as
*lambda_vdw*, *lamda_coulomb*, see :class:`mdpow.fep.Gsolv` for details) or by
deriving a new class from the :class:`mdpow.fep.Ghyd` base class but this is not
covered in this tutorial.)


Octanol
~~~~~~~

Equilibrium simulation
......................

Almost identical to the water case::

   O = mdpow.equil.OctanolSimulation(molecule="OcOH")
   O.topology(itp="octanol.itp")
   O.solvate(struct="octanol.gro")
   O.energy_minimize()
   O.MD_relaxed()
   O.MD(runtime=50, qscript=['my_script.sge', 'local.sh'])    # only run for 50 ps in this tutorial
   O.save()

.. note:: Here we only run 50 ps equilibrium MD for testing. For production this
          should be substantially longer, maybe even 50 ns if you want to
          extract thermodynamic data.


Octanol solvation free energy
.............................

Almost identical setup as in the water case::

   goct = mdpow.fep.Goct(simulation=O, runtime=100)
   goct.setup(qscript=['my_script.sge', 'local.sh'])

This generates all the input files under ``FEP/octanol``.

.. note:: Here we only run 100 ps per window for testing. For
          production this should be rather something like 5-10 ns (the
          default is 5 ns)


Running the FEP simulations
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The files are under the ``FEP/water`` and ``FEP/octanol`` directories
in separate sub directories. 

Either run job arrays that should have been generated from the
``my_script.sge`` template ::

   qsub Coul_my_script.sge
   qsub VDW_my_script.sge

Or run each job in its own directory. Note that :program:`mdrun`
should be called with at least the following options ::

   mdrun -deffnm $DEFFNM  -dgdl

where DEFFNM is typically "md"; see the run ``local.sh`` script in
each direcory for hints on what needs to be done.


Analyze output
~~~~~~~~~~~~~~

For the water and octanol FEPs do ::

 gwat.collect()
 gwat.analyze()

 goct.collect()
 goct.analyze()

The analyze step reports the estimate for the free energy
difference. 

Calculate the free energy for transferring the solute from water to
octanol and octanol-water partition coefficient log P_OW ::

 mdpow.fep.pOW(gwat, goct)

(see :func:`mdpow.fep.pOW` for details and definitions).

All individual results can also accessed as a dictionary ::

 gwat.results.DeltaA

Free energy of transfer from water to octanol::

 goct.results.DeltaA.total - gwat.results.DeltaA.total

The individual components are

 total
   total standard free energy difference in kJ/mol; 
   DeltaA0 = (A_solv - A_vac) + DA_v
 standardstate
   correction DA_v for the standard state concentration *c* (depends
   on the volume of the simulation cell); in kJ/mol. Different *c* can
   be supplied as an argument to analyze().
 coulomb
   contribution of the de-charging process to DeltaA
 vdw
   contribution of the de-coupling process to DeltaA

To plot the data (de-charging and de-coupling)::

 import pylab
 gwat.plot()
 pylab.figure()
 goct.plot()

The error data points are the average of dV/dl over each windows. The
error bars are the standard deviations of the data points from the
average.

If the graphs do not look smooth then a longer *runtime* is definitely
required. It might also be necessary to add additional lambda values
in regions where the function changes rapidly.


The mdpow scripts
=================

Some tasks are simplified by using scripts, which are installed in a
bin directory (or the directory pointed to by
``--install-scripts``). See :doc:`scripts` for details.

 
"""
VERSION = 0,2,0

__all__ = ['fep', 'equil']

import log
logger = log.create('mdpow', 'mdpow.log')

import config

def get_version():
    """Return current package version as a string."""
    return ".".join(map(str,VERSION))

def get_version_tuple():
    """Return current package version as a (MAJOR,MINOR,PATCHLEVEL)."""
    return tuple(VERSION)

# commented out so that one can get at version without importing
# GromacsWrapper etc
#import fep, equil
