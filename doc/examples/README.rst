.. -*- coding: utf-8 -*-

========================
 MDPOW Example: Benzene
========================

The ``examples`` directory contains input file for the solvation free
energy calculation of benzene in solvent (water and octanol), using
the OPLS-AA force field and the TIP4P water model.

There are two ways to use MDPOW:

1. commandline interface (CLI) through scripts
2. Python API


The online documentation goes into greater details. This README just
gives the briefest of summaries.



Included files
==============

The files in the examples directory are ::

    benzene.yml
    benzene/
    ├── README.txt
    ├── benzene.itp
    ├── benzene.pdb
    ├── benzene.pdb.png
    ├── dVdl.pdf
    └── session.py  

Run input file ``benzene.yml``
------------------------------

For running the simulations with the CLI, the run input file
``benzene.yml`` and the GROMACS input files ``benzene/benzene.itp``
and ``benzene/benzene.pdb`` are needed. The path to these files is
already in the runinput file so leave all files in place.


Benzene directory
-----------------

The ``benzene`` directory contains input files for GROMACS

* ``benzene.pdb``: coordinates in PDB format
* ``benzene.itp``: topology and parameters

an image of the chemical structure ``benzene.pdb.png``, output from
plotting dH/dlambda (for the MDPOW TI estimator) ``dVdl.pdf`` for
benzene in water.
  
The  ``README.txt`` file explains how to try out the Python API, with
the commands recorded in ``session.py``.



Using the CLI
=============

Make sure that GROMACS commands can be run (e.g., ``source GMXRC`` or
``module load gromacs``). All commands are run from the directory that
contains the run input file ``benzene.yml``.

To calculate the water solvation free energy::

  mdpow-equilibrium --solvent water benzene.yml
  mdpow-fep --solvent water benzene.yml
  mdpow-solvationenergy --solvent water benzene

See output in file ``energies.txt`` where all energies are in kJ/mol::

  ITP mol    solvent              dG   err_dG   Coulomb  err_Cou       VDW  err_VDW directory 
  ---------- ---------------- ------ -------- --------- -------- --------- -------- ---------
  BNZ        water             -8.32     1.49     +7.53     0.51     +0.78     1.40  benzene


To calculate the octanol solvation free energy::

  mdpow-equilibrium --solvent octanol benzene.yml
  mdpow-fep --solvent octanol benzene.yml
  mdpow-solvationenergy --solvent octanol benzene
 
The octanol energies are appended to the ``energies.txt`` file::

  BNZ        water             -8.32     1.49     +7.53     0.51     +0.78     1.40  benzene
  BNZ        octanol          -17.03     1.24     +1.93     0.25    +15.10     1.22  benzene  

Based on these *very* short simulations, the solvation free energy of benzene
in water of -8.32±1.49 kJ/mol is less favorable than the solvation free energy
in octanol of –17.03±1.24 kJ/mol.

Using these numbers, one can directly calculate the water-octanol partition
coefficient as

.. math::

   \log P_{OW} = (\Delta G_W - \Delta G_O)/kT \log e

or use a script ::

  mdpow-pow benzene/

which outputs ::

  mdpow.fep   : INFO     [BNZ] Free energy of transfer water --> octanol: -8.712 (1.941) kJ/mol
  mdpow.fep   : INFO     [BNZ] log P_ow: 1.517 (0.338)

and appends a line to the file ``pow.txt`` with content ::


  ITP mol    transferFE error          logPow err_logPow directory 
  ---------- --------- ------ ------ -------- ---------- --------- 
  BNZ         -8.71     1.94  logPow    +1.52     0.34   benzene

(It also writes output to ``energies.txt``. See online docs.)

Instead of running ``mdpow-solvationenergy`` for different solvents,
one can also just run ``mdpow-pow DIRECTORY`` and it will run the
separate free energy calculations automatically and then compute the
octanol-water partition coefficient.

  
The example CLI runs execute GROMACS immediately and run the FEP calculations
serially. For production calculations you will likely want to run these in
parallel on a computer cluster or HPC. The online documentation explains how to
change the runinput file to stop at various steps instead of running MD
simulations immediately.
  

Using the Python API
====================

See the ``benzene/README.txt`` and ``benzene/session.py`` files as well as the
online documentation.
