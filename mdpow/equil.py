# equil.py
# Copyright (c) 2010 Oliver Beckstein

"""
:mod:`mdpow.equil` --- Setting up and running equilibrium MD
============================================================

The :mod:`mdpow.equil` module facilitates the setup of equilibrium
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
from gromacs.utilities import in_dir, realpath, asiterable, AttributeDict

import config

logger = logging.getLogger('mdpow.equil')


class Simulation(object):
    """Simple MD simulation of a compound in water."""
    
    def __init__(self, name='DRUG', **kwargs):
        self.name = name
        self.dirs = AttributeDict(
            basedir=realpath(os.path.curdir),
            includes=list(asiterable(kwargs.pop('includes',[]))) + [config.includedir],
            )
        self.files = AttributeDict(
            topology=kwargs.pop('top', None),     # realpath?
            structure=kwargs.pop('struct', None), # realpath?
            solvated=kwargs.pop('solvated', None), # realpath?
            )

        if self.files.topology:
            # assume that a user-supplied topology lives in a 'standard' top dir
            # that includes the necessary itp file(s)
            self.dirs.topology = realpath(os.path.dirname(self.files.topology))
            self.dirs.includes.append(self.dirs.topology)

    def topology(self, itp='drug.itp', **kwargs):
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
        self.dirs.topology = realpath(dirname)
        
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

        # update known files and dirs
        self.files.topology = realpath(dirname, topol)
        if not self.dirs.topology in self.dirs.includes:
            self.dirs.includes.append(self.dirs.topology)
        
        return {'dirname': dirname, 'topol': topol}


    def solvate(self, **kwargs):
        """Solvate structure *struct* in a box of water.

        :Keywords:
          *struct*
              pdb or gro coordinate file (if not supplied, the value is used
              that was supplied to the constructor of :class:`~mdpow.equil.Simulation`)
          *kwargs*
              All other arguments are passed on to :func:`gromacs.setup.solvate`, but
              set to sensible default values. *top* is always fixed.
        """
        self.dirs.solvation = realpath(kwargs.setdefault('dirname', 'solvation'))        
        kwargs.setdefault('struct', self.files.structure)
        kwargs['top'] = self.files.topology
        kwargs.setdefault('water', 'tip4p')
        kwargs.setdefault('mainselection', '"%s"' % self.name)  # quotes are needed for make_ndx
        kwargs['includes'] = asiterable(kwargs.pop('includes',[])) + self.dirs.includes

        params = gromacs.setup.solvate(**kwargs)
        
        self.files.solvated = params['struct']
        return params

    def energy_minimize(self, **kwargs):
        """Energy minimize the solvatet structure on the local machine.

        *kwargs* are passed to :func:`gromacs.setup.energ_minimize` but if
        :meth:`~mdpow.equil.Simulation.solvate` step has been carried out
        previously all the defaults should just work.
        """

        self.dirs.energy_minimization = realpath(kwargs.setdefault('dirname', 'em'))
        kwargs['top'] = self.files.topology
        kwargs.setdefault('struct', self.files.solvated)
        # construct the complete include path
        kwargs['includes'] = asiterable(kwargs.pop('includes',[])) + self.dirs.includes

        params = gromacs.setup.energy_minimize(**kwargs)

        return params
