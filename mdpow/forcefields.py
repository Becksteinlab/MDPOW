# POW package __init__.py
# Copyright (c) 2012 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""
Configuration settings related to force fields
==============================================

At the moment, only the OPLS-AA force field is directly supported
(although in the principle it is possible to switch to a different
force field by supplying alternative template files). However, in the
future we want to support a simple configuration based switch.

.. autodata:: DEFAULT_FORCEFIELD

Different **water models** are already supported.

.. autodata:: DEFAULT_WATER_MODEL
.. autodata:: GROMACS_WATER_MODELS

"""

from __future__ import absolute_import

from collections import defaultdict

import logging
logger = logging.getLogger("mdpow.forecefields")


#: Default force field. At the moment, only OPLS-AA is directly
#: supported.
DEFAULT_FORCEFIELD = "OPLS-AA"

#------------------------------------------------------------
# Gromacs water models
#------------------------------------------------------------

#: See the file ``top/oplsaa.ff/watermodels.dat`` for a description of
#: available water models that are bundled with MDPOW.
GMX_WATERMODELS_DAT="""
tip4p   TIP4P  TIP 4-point, recommended
tip3p   TIP3P  TIP 3-point
tip5p   TIP5P  TIP 5-point
spc     SPC    simple point charge
spce    SPC/E  extended simple point charge
"""

class GromacsWaterModel(object):
    """Data for a water model."""
    def __init__(self, identifier, name=None, itp=None, coordinates=None,
                 description=None,
                 forcefield="OPLS-AA"):
        self.identifier = identifier
        self.name = name if name is not None else str(identifier).upper()
        self.itp = itp if itp is not None else self.guess_filename('itp')
        self.coordinates = coordinates if coordinates is not None else self.guess_filename('gro')
        self.description = description
        self.forcefield = forcefield

    def guess_filename(self, extension):
        return self.identifier.lower() + '.' + str(extension)

    def __repr__(self):
        return "<{0[name]} water: identifier={0[identifier]}, ff={0[forcefield]}>".format(vars(self))

#: For some water models we cannot derive the filename for the equilibrated box
#: so we supply them explicitly.
SPECIAL_WATER_COORDINATE_FILES = defaultdict(
    lambda: None,
    spc='spc216.gro',
    spce='spc216.gro',
    tip3p='spc216.gro',
)

def _create_water_models(watermodelsdat):
    models = {}
    for line in GMX_WATERMODELS_DAT.split('\n'):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        fields = line.split()
        identifier, name = fields[:2]
        description = " ".join(fields[2:])
        models[identifier] = GromacsWaterModel(identifier, name=name,
                                               coordinates=SPECIAL_WATER_COORDINATE_FILES[identifier],
                                               description=description)
    return models

#: Use the default water model unless another water model is chosen in the runinput
#: file (``setup.watermodel``).
DEFAULT_WATER_MODEL = "tip4p"

#: Dictionary of :class:`GromacsWaterModel` instances, one for each Gromacs water
#: model  available under the force field directory. The keys are the water model
#: identifiers.
#: For OPLS-AA the following ones are available.
GROMACS_WATER_MODELS = _create_water_models(GMX_WATERMODELS_DAT)

def get_water_model(watermodel=DEFAULT_WATER_MODEL):
    """Return a :class:`GromacsWaterModel` corresponding to identifier *watermodel*"""

    try:
        return GROMACS_WATER_MODELS[watermodel]
    except KeyError:
        msg = "{0} is not a valid water model: choose one from {1}".format(
            watermodel, ", ".join(GROMACS_WATER_MODELS.keys()))
        logger.error(msg)
        raise ValueError(msg)
