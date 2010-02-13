# equil.py
# Copyright (c) 2010 Oliver Beckstein

"""
:mod:`pow.equil` --- Setting up and running equilibrium MD
==========================================================

The :mod:`pow.equil` module facilitates the setup of equilibrium
molecular dynamics simulations of compound molecule in a simulation
box of water.

It requires as input
- the itp file for the compound
- a coordinate (structure) file (in pdb or gro format)

By default it uses the *OPLS/AA* forcefield and the *TIP4P* water
model.
"""

from __future__ import with_statement

import os
import shutil
import logging

import gromacs.setup
import gromacs.cbook
from gromacs.utilities import in_dir, realpath, asiterable

import config

logger = logging.getLogger('pow.equil')


class Simulation(object):
    """Simple MD simulation of a compound in water."""
    
    def __init__(self, name='DRUG', **kwargs):
        self.name = name
        self.dirs = {'basedir': realpath(os.path.curdir),
                     'includes': list(asiterable(kwargs.pop('includes',[]))) + [config.includedir],
                     }
        self.topology = kwargs.pop('top', None)

    def create_topology(self, itp='drug.itp', **kwargs):
        """Generate a topology for compound *name*.

        :Keywords:
           *name*
               identifier, should be the same as used in the itp
            *struct*
               coordinate file (gro or pdb) (is copied to topology dir)
            *itp*
               Gromacs itp file; will be copied to topology dir and
               included in topology
            *dirname*
               name of the topology directory ["top"]
            *kwargs*
               see source for *top_template*, *topol*
        """
        dirname = kwargs.pop('dirname', 'top')
        top_template = config.get_template(kwargs.pop('top_template', 'system.top'))
        topol = kwargs.pop('topol', os.path.basename(top_template))
        itp = os.path.realpath(itp)
        _itp = os.path.basename(itp)

        with in_dir(dirname):
            shutil.copy(itp, _itp)
            gromacs.cbook.edit_txt(top_template,
                                   [('#include', 'compound.itp', _itp),
                                    ('Compound', 'DRUG', self.name),
                                    ('DRUG\s*1', 'DRUG', self.name),
                                    ],
                                   newname=topol)
        logger.info('[%(dirname)s] Created topology %(topol)r that includes %(_itp)r', vars())

        self.topology = realpath(dirname, topol)
        self.dirs['topology'] = realpath(dirname)
        
        return {'dirname': dirname, 'topol': topol}


    def solvate(self, struct='drug.pdb', **kwargs):
        """Solvate structure *struct* in a dodecahedral box of water.

        *kwargs* are passed on to :func:`gromacs.setup.solvate`.
        """
        self.structure = realpath(struct)
        kwargs['struct'] = struct
        kwargs['top'] = self.topology
        kwargs.setdefault('water', 'tip4p')
        kwargs.setdefault('mainselection', self.name)
        kwargs.setdefault('dirname', 'solvation')
        kwargs['includes'] = self.dirs['includes']  # XXX: add user ones
        self.dirs['solvation'] = realpath(kwargs['dirname'])

        return gromacs.setup.solvate(**kwargs)
