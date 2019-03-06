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
.. autodata:: DEFAULT_WATER_MODEL


Solvent models
--------------

Different **water models** are already supported

.. autodata:: GROMACS_WATER_MODELS

as well as different general **solvent models**

.. autodata:: GROMACS_SOLVENT_MODELS


Internal data
-------------

.. autodata:: SPECIAL_WATER_COORDINATE_FILES
.. autodata:: GROMACS_WATER_MODELS
.. autodata:: GROMACS_SOLVENT_MODELS

Internal classes and functions
------------------------------

.. autoclass:: GromacsSolventModel
   :members:

.. autofunction:: get_water_model

.. autofunction:: get_solvent_identifier

.. autofunction:: get_solvent_model
"""

from __future__ import absolute_import

import os
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
tip4p       TIP4P      TIP 4-point, recommended
tip3p       TIP3P      TIP 3-point
tip5p       TIP5P      TIP 5-point
spc         SPC        simple point charge
spce        SPC/E      extended simple point charge
m24         M24        TIP 3-point with modified LJ (M24)
"""

class GromacsSolventModel(object):
    """Data for a solvent model in Gromacs."""
    def __init__(self, identifier, name=None, itp=None, coordinates=None,
                 description=None, forcefield="OPLS-AA"):
        self.identifier = identifier
        self.name = name if name is not None else str(identifier).upper()
        self.itp = itp if itp is not None else self.guess_filename('itp')
        self.coordinates = coordinates if coordinates is not None else self.guess_filename('gro')
        self.description = description
        self.forcefield = forcefield

    def guess_filename(self, extension):
        """Guess the filename for the model and add *extension*"""
        return self.identifier.lower() + os.extsep + str(extension)

    def __repr__(self):
        return "<{0[name]} water: identifier={0[identifier]}, ff={0[forcefield]}>".format(vars(self))

#: For some water models we cannot derive the filename for the equilibrated box
#: so we supply them explicitly.
SPECIAL_WATER_COORDINATE_FILES = defaultdict(
    lambda: None,
    spc='spc216.gro',
    spce='spc216.gro',
    tip3p='spc216.gro',
    m24='spc216.gro',
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
        models[identifier] = GromacsSolventModel(
            identifier, name=name,
            coordinates=SPECIAL_WATER_COORDINATE_FILES[identifier],
            description=description)
    return models

#: Use the default water model unless another water model is chosen in the runinput
#: file (``setup.watermodel``).
DEFAULT_WATER_MODEL = "tip4p"

#: Dictionary of :class:`GromacsSolventModel` instances, one for each Gromacs water
#: model  available under the force field directory. The keys are the water model
#: identifiers.
#: For OPLS-AA the following ones are available.
GROMACS_WATER_MODELS = _create_water_models(GMX_WATERMODELS_DAT)

def get_water_model(watermodel=DEFAULT_WATER_MODEL):
    """Return a :class:`GromacsSolventModel` corresponding to identifier *watermodel*"""

    try:
        return GROMACS_WATER_MODELS[watermodel]
    except KeyError:
        msg = "{0} is not a valid water model: choose one from {1}".format(
            watermodel, ", ".join(GROMACS_WATER_MODELS.keys()))
        logger.error(msg)
        raise ValueError(msg)

#: Other solvents (not water, see :data:`GROMACS_WATER_MODELS` for those).
OPLS_SOLVENT_MODELS = {
    'octanol': GromacsSolventModel(
        identifier="octanol", itp="1oct.itp", coordinates="1oct.gro"),
    'cyclohexane': GromacsSolventModel(
        identifier="cyclohexane", itp="1cyclo.itp", coordinates="1cyclo.gro"),
    'wetoctanol': GromacsSolventModel(
        identifier="wetoctanol", itp="1octwet.itp", coordinates="1octwet.gro"),
    }

CHARMM_SOLVENT_MODELS = {
    'octanol': GromacsSolventModel(
        identifier="octanol", itp="1oct.itp", coordinates="1oct_charmm.gro"),
    'wetoctanol': GromacsSolventModel(
        identifier="wetoctanol", itp="1octwet.itp", coordinates="1octwet_charmm.gro"),
    }

AMBER_SOLVENT_MODELS = {
    'octanol': GromacsSolventModel(
        identifier="octanol", itp="1oct.itp", coordinates="1oct_charmm.gro"),
    'wetoctanol': GromacsSolventModel(
        identifier="wetoctanol", itp="1octwet.itp", coordinates="1octwet_amber.gro"),
    }

GROMACS_SOLVENT_MODELS = {
    'OPLS-AA': OPLS_SOLVENT_MODELS,
    'CHARMM': CHARMM_SOLVENT_MODELS,
    'AMBER': AMBER_SOLVENT_MODELS,
    }

def get_solvent_identifier(solvent_type, model=None, forcefield='OPLS-AA'):
    """Get the identifier for a solvent model.

    The identifier is needed to access a water model (i.e., a
    :class:`GromacsSolventModel`) through
    :func:`get_solvent_model`. Because we have multiple water models
    but only limited other solvents, the organization of these models
    is a bit convoluted and it is best to obtain the desired water
    model in these two steps::

      identifier = get_solvent_identifier("water", model="tip3p")
      model = get_solvent_model(identifier)


    For ``solvent_type`` *water*: either provide ``None`` or "water"
    for the specific ``model`` (and the default
    :data:`DEFAULT_WATER_MODEL` will be selected, or a specific water
    model such as "tip3p" or "spce" (see
    :data:`GROMACS_WATER_MODELS`). For other solvent types, the
    ``model`` is ignored.

    :Returns: Either an identifier or ``None``

    """
    if solvent_type is "water":
        identifier = model if not model in (None, 'water') else DEFAULT_WATER_MODEL
        return identifier if identifier in GROMACS_WATER_MODELS else None
    return solvent_type if solvent_type in GROMACS_SOLVENT_MODELS[forcefield] else None


def get_solvent_model(identifier, forcefield='OPLS-AA'):
    """Return a :class:`GromacsSolventModel` corresponding to identifier *identifier*.

    If identifier is "water" then the :data:`DEFAULT_WATER_MODEL` is assumed.
    """

    if identifier is "water":
        identifier = DEFAULT_WATER_MODEL
    try:
        return GROMACS_WATER_MODELS[identifier]
    except KeyError:
        try:
            return GROMACS_SOLVENT_MODELS[forcefield][identifier]
        except KeyError:
            msg = "No solvent model with name {0} is available.".format(identifier)
            logger.critical(msg)
            raise ValueError(msg)

