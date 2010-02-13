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

"""

__all__ = ['ghyd', 'equil']

import log
logger = log.create('mdpow', 'mdpow.log')

import config, ghyd, equil

