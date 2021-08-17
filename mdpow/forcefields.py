# -*- coding: utf-8 -*-
# POW package __init__.py
# Copyright (c) 2012 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""
Force field selection
=====================

The :mod:`mdpow.forcefields` module contains settings for selecting
different force fields and the corresponding solvent topologies.

The OPLS-AA, CHARMM/CGENFF and the AMBER/GAFF force field are directly
supported. In the principle it is possible to switch to a
different force field by supplying alternative template
files.

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
   :noindex:
.. autodata:: GROMACS_SOLVENT_MODELS
   :noindex:

Internal classes and functions
------------------------------

.. autoclass:: GromacsSolventModel
   :members:

.. autofunction:: get_water_model

.. autofunction:: get_solvent_identifier

.. autofunction:: get_solvent_model

"""
import os
from collections import defaultdict

import logging
logger = logging.getLogger("mdpow.forecefields")

#: Default force field. At the moment, OPLS-AA, CHARMM/CGENFF, and AMBER/GAFF
#: are directly supported. However, it is not recommended to change the
#: default here as this behavior is not tested.
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
tip4pd      TIP4P-D    TIP 4-point with modified dispersion (TIP4P-D)
tip4pew     TIP4PEW    TIP 4-point modified for use with Ewald techniques (TIP4PEW)
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
    tip4pd='tip4p.gro',
    tip4pew='tip4p.gro',
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

#: Use the default water model unless another water model is chosen in the
#: :ref:`run input file <runinput-file>` file by setting the
#: ``setup.watermodel`` parameter.
#:
#: .. warning::
#:    Select the native water model **manually** and do not rely on the default
#:    set here. For CHARMM/GAFF the CHARMM TIP3P model is recommended.
#:    For AMBER/GAFF the TIP3P mdeol is often used. Choosing the correct water model
#:    is a *scientific* decision that *you* have to make conscientiously.
#:
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
new_octanol = '''Zangi R (2018) Refinement of the OPLSAA force-field
                for liquid alcohols.; ACS Omega 3(12):18089-18099.
                doi: 10.1021/acsomega.8b03132'''

OPLS_SOLVENT_MODELS = {
    'octanol': GromacsSolventModel(
        identifier="octanol", itp="1oct.itp", coordinates="1oct.gro"),
    'octanolnew': GromacsSolventModel(
        identifier="octanol", itp="1octnew.itp", coordinates="1oct.gro",
        description=new_octanol),
    'cyclohexane': GromacsSolventModel(
        identifier="cyclohexane", itp="1cyclo.itp", coordinates="1cyclo.gro"),
    'wetoctanol': GromacsSolventModel(
        identifier="wetoctanol", itp="1octwet.itp", coordinates="1octwet.gro"),
    'wetoctanolnew': GromacsSolventModel(
        identifier="wetoctanol", itp="1octwetnew.itp", coordinates="1octwet.gro",
        description=new_octanol)
    }

CHARMM_SOLVENT_MODELS = {
    'octanol': GromacsSolventModel(
        identifier="octanol", itp="1oct.itp", coordinates="1oct_charmm.gro"),
    'wetoctanol': GromacsSolventModel(
        identifier="wetoctanol", itp="1octwet.itp", coordinates="1octwet_charmm.gro"),
    'cyclohexane': GromacsSolventModel(
        identifier="cyclohexane", itp="1cyclo.itp", coordinates="1cyclo_charmm.gro"),
    }

AMBER_SOLVENT_MODELS = {
    'octanol': GromacsSolventModel(
        identifier="octanol", itp="1oct.itp", coordinates="1oct_amber.gro"),
    'wetoctanol': GromacsSolventModel(
        identifier="wetoctanol", itp="1octwet.itp", coordinates="1octwet_amber.gro"),
    'cyclohexane': GromacsSolventModel(
        identifier="cyclohexane", itp="1cyclo.itp", coordinates="1cyclo_amber.gro"),
    }

#: Solvents available in GROMACS; the keys of the dictionary
#: are the forcefields.
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
    :data:`GROMACS_WATER_MODELS`). For other "octanol" or "wetoctanol"
    of OPLS-AA forcefield, the ``model`` is used to select a specific
    model. For other solvents and forcefields, "model" is not required.

    :Returns: Either an identifier or ``None``
    """
    
    if solvent_type == "water":
        identifier = model if not model in (None, 'water') else DEFAULT_WATER_MODEL
        return identifier if identifier in GROMACS_WATER_MODELS else None
    if model not in GROMACS_SOLVENT_MODELS[forcefield]:
        if solvent_type in GROMACS_SOLVENT_MODELS[forcefield]:
            model = solvent_type
        else:
            model = None
    return model


def get_solvent_model(identifier, forcefield='OPLS-AA'):
    """Return a :class:`GromacsSolventModel` corresponding to identifier *identifier*.

    If identifier is "water" then the :data:`DEFAULT_WATER_MODEL` is assumed.
    """

    if identifier == "water":
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


def get_ff_paths(forcefield='OPLS-AA'):
    """Return a :list: containing the forcefield directory, paths of ions
    and default watermodel itp files.
    """
    
    settings = {
                'OPLS-AA': ['oplsaa.ff/', 'oplsaa.ff/ions_opls.itp',
                            'oplsaa.ff/tip4p.itp'],
                'AMBER': ['amber99sb.ff/', 'amber99sb.ff/ions.itp',
                           'amber99sb.ff/tip3p.itp'],
                'CHARMM': ['charmm36-mar2019.ff/', 'charmm36-mar2019.ff/ions.itp',
                          'charmm36-mar2019.ff/tip3p.itp'],
                }
    try:
        return settings[forcefield]
    except KeyError:
        msg = "No forcefield with name {0} is available".format(forcefield)
        logger.critical(msg)
        raise ValueError(msg)


def get_top_template(identifier):
    """Return the topology file template suitable for the solvent model."""
      
    templates = {'water': 'system.top', 'octanol': 'system.top',
                 'cyclohexane': 'system.top', 'wetoctanol': 'system_octwet.top'}
    try:
        return templates[identifier]
    except KeyError:
        msg = "No template for solvent {0} is available".format(identifier)
        logger.critical(msg)
        raise ValueError(msg)
