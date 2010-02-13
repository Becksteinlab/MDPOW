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
import cPickle
import logging

import gromacs.setup
import gromacs.cbook
from gromacs.utilities import in_dir, realpath, asiterable, AttributeDict

import config

logger = logging.getLogger('mdpow.equil')


class Simulation(object):
    """Simple MD simulation of a single compound molecule in water.

    Typical use ::
       S = Simulation(name='DRUG')
       S.topology(itp='drug.itp')
       S.solvate(struct='DRUG-H.pdb')
       S.energy_minimize()

    .. Note:: The OPLS/AA force field and the TIP4P water molecule is the
              default; changing this is possible but will require provision of
              customized itp and mdp files at various stages.
    """
    
    def __init__(self, name=None, filename=None, **kwargs):
        """Set up Simulation instance.

        The *name* of the compound molecule should be supplied. Existing files
        (which have been generated in previous runs) can also be supplied.

        :Keywords:
          *name*
              Identifier for the compound molecule. This is the same as the
              entry in the ``[ molecule ]`` section of the itp file. ["DRUG"]
          *filename*
              If provided and *name* is ``None`` then load the instance from
              the pickle file *filename*, which was generated with
              :meth:`~mdpow.equil.Simulation.save`.
          *kwargs*
              advanced keywords for short-circuiting; see the source
        """
        if name is None and not filename is None:
            # load from pickle file
            self.load(filename)
            self.filename = filename            
            kwargs = {}    # for super
        else:
            self.name = name or 'DRUG'
            self.dirs = AttributeDict(
                basedir=realpath(os.path.curdir),
                includes=list(asiterable(kwargs.pop('includes',[]))) + [config.includedir],
                )
            self.files = AttributeDict(
                topology=kwargs.pop('top', None),     # realpath?
                processed_topology = None,
                structure=kwargs.pop('struct', None), # realpath?
                solvated=kwargs.pop('solvated', None), # realpath?
                ndx=kwargs.pop('ndx', None),
                energy_minimized=kwargs.pop('energy_minimized', None),
                MD_restrained=kwargs.pop('MD_restrained', None),
                )

            if self.files.topology:
                # assume that a user-supplied topology lives in a 'standard' top dir
                # that includes the necessary itp file(s)
                self.dirs.topology = realpath(os.path.dirname(self.files.topology))
                self.dirs.includes.append(self.dirs.topology)

        super(Simulation, self).__init__(**kwargs)

    def save(self, filename=None):
        """Save instance to a pickle file.

        The default filename is the name of the file that was last loaded from
        or saved to.
        """
        if filename is None:
            filename = self.filename
        else:
            self.filename = filename
        with open(filename, 'wb') as f:
            cPickle.dump(self, f, protocol=cPickle.HIGHEST_PROTOCOL)
        logger.debug("Instance pickled to %(filename)r" % vars())
        
    def load(self, filename=None):
        """Re-instantiate class from pickled file."""
        if filename is None:
            filename = self.filename        
        with open(filename, 'rb') as f:
            instance = cPickle.load(f)
        self.__dict__.update(instance.__dict__)
        logger.debug("Instance loaded from %(filename)r" % vars())


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

        self.files.structure = kwargs['struct']
        self.files.solvated = params['struct']
        self.files.ndx = params['ndx']

        # we can also make a processed topology right now
        self.processed_topology()
        
        return params

    def processed_topology(self):
        """Create a portable topology file from the topology and the solvated system."""
        if not os.path.exists(self.files.solvated):
            self.solvate()
        self.files.processed_topology = gromacs.cbook.create_portable_topology(
            self.files.topology, self.files.solvated, includes=self.dirs.includes)
        return self.files.processed_topology

    def energy_minimize(self, **kwargs):
        """Energy minimize the solvated structure on the local machine.

        *kwargs* are passed to :func:`gromacs.setup.energ_minimize` but if
        :meth:`~mdpow.equil.Simulation.solvate` step has been carried out
        previously all the defaults should just work.
        """

        self.dirs.energy_minimization = realpath(kwargs.setdefault('dirname', 'em'))
        kwargs['top'] = self.files.topology
        kwargs.setdefault('struct', self.files.solvated)
        kwargs['includes'] = asiterable(kwargs.pop('includes',[])) + self.dirs.includes

        params = gromacs.setup.energy_minimize(**kwargs)

        self.files.energy_minimized = params['struct']
        return params

    def MD_restrained(self, **kwargs):
        """Short MD simulation with position restraints on compound.

        See documentation of :func:`gromacs.setup.MD_restrained` for
        details. The following keywords can not be changed: top, mdp, ndx,
        mainselection

        .. Note:: Position restraints are activated with ``-DPOSRES`` directives
                  for :func:`gromacs.grompp`. Hence this will only work if the
                  compound itp file does indeed contain a ``[ posres ]``
                  section that is protected by a ``#ifdef POSRES`` clause.
        """
        self.dirs.MD_restrained = realpath(kwargs.setdefault('dirname', 'MD_restrained'))
        kwargs['top'] = self.files.topology
        kwargs.setdefault('struct', self.files.energy_minimized)
        kwargs['includes'] = asiterable(kwargs.pop('includes',[])) + self.dirs.includes
        kwargs['mdp'] = config.get_template('NPT_opls.mdp')
        kwargs['ndx'] = self.files.ndx
        kwargs['mainselection'] = None    # important for SD (use custom mdp and ndx!)

        params = gromacs.setup.MD_restrained(**kwargs)

        self.files.MD_restrained = params['struct']
        return params

    def MD(self, **kwargs):
        """Short NPT MD simulation.

        See documentation of :func:`gromacs.setup.MD` for
        details. The following keywords can not be changed: top, mdp, ndx,
        mainselection
        """
        self.dirs.MD_NPT = realpath(kwargs.setdefault('dirname', 'MD_NPT'))
        kwargs['top'] = self.files.topology
        kwargs.setdefault('struct', self.files.MD_restrained)
        kwargs['includes'] = asiterable(kwargs.pop('includes',[])) + self.dirs.includes
        kwargs['mdp'] = config.get_template('NPT_opls.mdp')
        kwargs['ndx'] = self.files.ndx
        kwargs['mainselection'] = None    # important for SD (use custom mdp and ndx!)

        params = gromacs.setup.MD(**kwargs)

        self.files.MD_NPT = params['struct']
        return params
        
