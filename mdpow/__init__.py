# POW package __init__.py
# Copyright (c) 2010 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""
:mod:`mdpow` --- Computing the octanol/water partitioning coefficient
=====================================================================

The :mod:`mdpow` module helps in setting up and analyzing absolute
free energy calculations of small molecules by molecular dynamics (MD)
simulations. By computing the hydration free energy and the solvation
free energy in octanol one can compute the octanol/water partitioning
coefficient, and important quantity that is used to characterize
drug-like compounds.

The MD simulations are performed with Gromacs_ 4.x

.. _Gromacs: http://www.gromacs.org


How to use the module
=====================

Before you can start you will need

  1. a coordinate file for the small molecule
  2. a Gromacs OPLS/AA topology (itp) file

Basic work flow
---------------

Then you will typically calculate two solvation free energies (free
energy of transfer of the olute from the liquid into the vacuum
phase):

  1. solvent = water

  1.1 set up a short equilibrium simulation of the molecule in a *water*
      box (and run the MD simulation);

  1.2 set up a free energy perturbation calculation of the ligand in
      water , which will yield the hydration free energy;

  2. solvent = octanol

  2.1 set up a short equilibrium simulation of the molecule in a *octanol*
      box (and run the MD simulation);

  2.2 set up a free energy perturbation calculation of the ligand in
      octanol , which will yield the solvation free energy in octanol;

  3. run these simulations on a cluster;

  4. analyze the output and combine the free energies to arrive at an
     estimate of the octanol-water partition coefficient.
    

Example session: Octanol
------------------------

In the following interactive python session we use octanol as an
example; all files are present in the package so one can work through
the example immediately.

Before starting python (preferrably ipython) make sure that the
gromacs tools can be found, e.g. ``which grompp`` should show you the
path to :program:`grompp`.


Water
~~~~~

Equilibrium simulation
......................

Make a directory ``octanol`` and copy the ``octanol.itp`` and
``octanol.gro`` file into it. Launch :program:`python` from this
directory and type::

   import mdpow
   S = mdpow.equil.WaterSimulation(molecule="OcOH")
   S.topology(itp="octanol.itp")
   S.solvate(struct="octanol.gro")
   S.energy_minimize()
   S.MD(runtime=50)            # only run for 50 ps in this tutorial
   S.save("water.simulation")  # save setup for later (analysis stage)

Background (Ctrl-Z) or quit (Ctrl-D) python and run the simulation in
the ``MD_NPT`` subdirectory. You can modify the ``local.sh`` script to
your ends. (One can also generate run scripts for various queuing
systems, check the documentation for :func:`gromacs.setup.MD` for
details).


Hydration free energy
.....................

Reopen the python session (if you quit it you will have to ``import
mdpow`` again) and set up a Ghyd object::

   G = mdpow.ghyd.Ghyd(molecule="OcOH", top="Equilibrium/water/top/system.top", struct="Equilibrium/water/MD_NPT/md.pdb", ndx="Equilibrium/water/solvation/main.ndx")
   
Alternatively, one can save some typing if we continue the last
session and use the Simulation object::

  S = mdpow.equil.WaterSimulation(filename="water.simulation")  # only needed when quit
  G = mdpow.ghyd.Ghyd(simulation=S)

Then set up all input files::

  G.setup()

(The details of the FEP runs can be customized by setting some
keywords (lambda_vdw, lamda_coulomb) or by deriving a new class from
the :class:`mdpow.ghyd.Ghyd` base class but this is not covered in
this tutorial.)


Octanol
~~~~~~~

Equilibrium simulation
......................

XXX: need to have a different molname than OcOH for the octanol that
     is used as a solvent!

Almost identical to the water case

   H = mdpow.equil.OctanolSimulation(molecule="OcOH")
   H.topology(itp="octanol.itp")
   H.solvate(struct="octanol.gro")
   H.energy_minimize()
   H.MD(runtime=50)          # only run for 50 ps in this tutorial


Octanol solvation free energy
.............................

TODO

 H = mdpow.gsolv.Goct(simulation=S)


Running the FEP simulations (4)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The files are under the ``FEP/water`` and ``FEP/octanol`` directories
in separate sub directories. Run each job in its own directory. Note
that :program:`mdrun` should be called with at least the following
options ::

     mdrun -deffnm $DEFFNM  -dgdl

where DEFFNM is typically "md"; see the run script in each direcory
for hints,


Analyze output (5)
~~~~~~~~~~~~~~~~~~




"""

__all__ = ['ghyd', 'equil']

import log
logger = log.create('mdpow', 'mdpow.log')

import config, ghyd, equil

