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

.. autoclass:: Forcefield
   :members:

.. autofunction:: get_water_model

.. autofunction:: get_solvent_identifier

.. autofunction:: get_solvent_model

"""
import os
from dataclasses import dataclass
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Union

import logging

logger = logging.getLogger("mdpow.forecefields")

HERE = Path(__file__).parent
TOP_DIR = HERE / "top"

#: Default force field. At the moment, OPLS-AA, CHARMM/CGENFF, and AMBER/GAFF
#: are directly supported. However, it is not recommended to change the
#: default here as this behavior is not tested.
DEFAULT_FORCEFIELD = "OPLS-AA"


@dataclass
class GromacsSolventModel(object):
    """Data for a solvent model in Gromacs."""

    def __init__(
        self,
        identifier,
        name=None,
        itp=None,
        coordinates=None,
        description=None,
        forcefield="OPLS-AA",
    ):
        self.identifier = identifier
        self.name = name if name is not None else str(identifier).upper()
        self.itp = itp if itp is not None else self.guess_filename("itp")
        self.coordinates = (
            coordinates if coordinates is not None else self.guess_filename("gro")
        )
        self.description = description
        self.forcefield = forcefield

    def guess_filename(self, extension):
        """Guess the filename for the model and add *extension*"""
        return self.identifier.lower() + os.extsep + str(extension)

    def __repr__(self):
        return (
            "<{0[name]} water: identifier={0[identifier]}, ff={0[forcefield]}>".format(
                vars(self)
            )
        )


# ------------------------------------------------------------
# Gromacs water models
# ------------------------------------------------------------

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
#: For some water models we cannot derive the filename for the equilibrated box
#: so we supply them explicitly.
SPECIAL_WATER_COORDINATE_FILES = defaultdict(
    lambda: None,
    spc="spc216.gro",
    spce="spc216.gro",
    tip3p="spc216.gro",
    m24="spc216.gro",
    tip4pd="tip4p.gro",
    tip4pew="tip4p.gro",
)

def _create_water_models(watermodelsdat: Path) -> Dict[str, GromacsSolventModel]:
    models = {}
    data = watermodelsdat.read_text()
    for line in data.split('\n'):
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

#: See the file ``top/oplsaa.ff/watermodels.dat`` for a description of
#: available water models that are bundled with MDPOW.

#: Dictionary of :class:`GromacsSolventModel` instances, one for each Gromacs water
#: model  available under the force field directory. The keys are the water model
#: identifiers.
#: For OPLS-AA the following ones are available.
GROMACS_WATER_MODELS: Dict[str, GromacsSolventModel] = {
    "tip4p": GromacsSolventModel(
        "tip4p",
        name="TIP4P",
        coordinates=SPECIAL_WATER_COORDINATE_FILES["tip4p"],
        description="TIP 4-point, recommended",
    ),
    "tip3p": GromacsSolventModel(
        "tip3p",
        name="TIP3P",
        coordinates=SPECIAL_WATER_COORDINATE_FILES["tip3p"],
        description="TIP 3-point",
    ),
    "tip5p": GromacsSolventModel(
        "tip5p",
        name="TIP5P",
        coordinates=SPECIAL_WATER_COORDINATE_FILES["tip5p"],
        description="TIP 5-point",
    ),
    "spc": GromacsSolventModel(
        "spc",
        name="SPC",
        coordinates=SPECIAL_WATER_COORDINATE_FILES["spc"],
        description="simple point charge",
    ),
    "spce": GromacsSolventModel(
        "spce",
        name="SPC/E",
        coordinates=SPECIAL_WATER_COORDINATE_FILES["spce"],
        description="extended simple point charge",
    ),
    "m24": GromacsSolventModel(
        "m24",
        name="M24",
        coordinates=SPECIAL_WATER_COORDINATE_FILES["m24"],
        description="TIP 3-point with modified LJ (M24)",
    ),
    "tip4pd": GromacsSolventModel(
        "tip4pd",
        name="TIP4P-D",
        coordinates=SPECIAL_WATER_COORDINATE_FILES["tip4pd"],
        description="TIP 4-point with modified dispersion (TIP4P-D)",
    ),
    "tip4pew": GromacsSolventModel(
        "tip4pew",
        name="TIP4PEW",
        coordinates=SPECIAL_WATER_COORDINATE_FILES["tip4pew"],
        description="TIP 4-point modified for use with Ewald techniques (TIP4PEW)",
    ),
}

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


def get_water_model(watermodel=DEFAULT_WATER_MODEL):
    """Return a :class:`GromacsSolventModel` corresponding to identifier *watermodel*"""

    try:
        return GROMACS_WATER_MODELS[watermodel]
    except KeyError:
        msg = "{0} is not a valid water model: choose one from {1}".format(
            watermodel, ", ".join(GROMACS_WATER_MODELS.keys())
        )
        logger.error(msg)
        raise ValueError(msg)


@dataclass
class Forcefield:
    """Contains information about files corresponding to a forcefield."""

    name: str
    solvent_models: Dict[str, GromacsSolventModel]
    forcefield_dir: Path
    ions_itp: Path
    default_water_itp: Path

    def __post_init__(self):
        """Check that the provided paths exist and populate the :attr:`water_models`."""
        for path_ in self.ff_paths:
            if not path_.exists():
                raise ValueError(f"Could not find required forcefield files: {path_}")
        
        self.water_models: Dict[str, GromacsSolventModel] = _create_water_models(self.forcefield_dir / "watermodels.dat")

    @property
    def ff_paths(self) -> List[Path]:
        """Get the path to the forcefield directory and the ITP files for the ions and default water model."""
        return [self.forcefield_dir, self.ions_itp, self.default_water_itp]
    
    def __repr__(self) -> str:
        """Get the name of the force field."""
        return self.name


#: Other solvents (not water, see :data:`GROMACS_WATER_MODELS` for those).
NEW_OCTANOL_DESC = "Zangi R (2018) Refinement of the OPLSAA force-field for liquid alcohols.; ACS Omega 3(12):18089-18099.  doi: 10.1021/acsomega.8b03132"

OPLS_AA_FF_DIR = TOP_DIR / "oplsaa.ff"
OPLS_AA = Forcefield(
    "OPLS-AA",
    {
        "octanol": GromacsSolventModel(
            identifier="octanol", itp="1oct.itp", coordinates="1oct.gro"
        ),
        "octanolnew": GromacsSolventModel(
            identifier="octanol",
            itp="1octnew.itp",
            coordinates="1oct.gro",
            description=NEW_OCTANOL_DESC,
        ),
        "cyclohexane": GromacsSolventModel(
            identifier="cyclohexane", itp="1cyclo.itp", coordinates="1cyclo.gro"
        ),
        "wetoctanol": GromacsSolventModel(
            identifier="wetoctanol", itp="1octwet.itp", coordinates="1octwet.gro"
        ),
        "wetoctanolnew": GromacsSolventModel(
            identifier="wetoctanol",
            itp="1octwetnew.itp",
            coordinates="1octwet.gro",
            description=NEW_OCTANOL_DESC,
        ),
        "toluene": GromacsSolventModel(
            identifier="toluene", itp="1tol.itp", coordinates="1tol_oplsaa.gro"
        ),
    },
    OPLS_AA_FF_DIR,
    OPLS_AA_FF_DIR / "ions_opls.itp",
    OPLS_AA_FF_DIR / "tip4p.itp",
)

CHARMM_FF_DIR = TOP_DIR / "charmm36-mar2019.ff"
CHARMM = Forcefield(
    "CHARMM",
    {
        "octanol": GromacsSolventModel(
            identifier="octanol", itp="1oct.itp", coordinates="1oct_charmm.gro"
        ),
        "wetoctanol": GromacsSolventModel(
            identifier="wetoctanol", itp="1octwet.itp", coordinates="1octwet_charmm.gro"
        ),
        "cyclohexane": GromacsSolventModel(
            identifier="cyclohexane", itp="1cyclo.itp", coordinates="1cyclo_charmm.gro"
        ),
        "toluene": GromacsSolventModel(
            identifier="toluene", itp="1tol.itp", coordinates="1tol_charmm.gro"
        ),
    },
    CHARMM_FF_DIR,
    CHARMM_FF_DIR / "ions.itp",
    CHARMM_FF_DIR / "tip3p.itp",
)

AMBER_FF_DIR = TOP_DIR / "amber99sb.ff"
AMBER = Forcefield(
    "AMBER",
    {
        "octanol": GromacsSolventModel(
            identifier="octanol", itp="1oct.itp", coordinates="1oct_amber.gro"
        ),
        "wetoctanol": GromacsSolventModel(
            identifier="wetoctanol", itp="1octwet.itp", coordinates="1octwet_amber.gro"
        ),
        "cyclohexane": GromacsSolventModel(
            identifier="cyclohexane", itp="1cyclo.itp", coordinates="1cyclo_amber.gro"
        ),
        "toluene": GromacsSolventModel(
            identifier="toluene", itp="1tol.itp", coordinates="1tol_amber.gro"
        ),
    },
    AMBER_FF_DIR,
    AMBER_FF_DIR / "ions.itp",
    AMBER_FF_DIR / "tip3p.itp",
)

ALL_FORCEFIELDS: Dict[str, Forcefield] = {
    "OPLS-AA": OPLS_AA,
    "CHARMM": CHARMM,
    "AMBER": AMBER,
}

#: Solvents available in GROMACS; the keys of the dictionary
#: are the forcefields.
GROMACS_SOLVENT_MODELS = {
    "OPLS-AA": OPLS_AA.solvent_models,
    "CHARMM": CHARMM.solvent_models,
    "AMBER": AMBER.solvent_models,
}

def _get_forcefield(ff_name: str) -> Forcefield:
    """Get the :class:`Forcefield` instance corresponding to a given name."""
    try:
        return ALL_FORCEFIELDS[ff_name]
    except KeyError:
        raise ValueError(f"Forcefield `{ff_name}` not found.")

def get_solvent_identifier(solvent_type, model=None, forcefield: Union[Forcefield, str] = OPLS_AA):
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

    :Raises ValueError: If there is no identifier found for the combination.

    :Returns: An identifier
    """
    if isinstance(forcefield, str):
        forcefield = _get_forcefield(forcefield)

    if solvent_type == "water":
        identifier = model if not model in (None, "water") else DEFAULT_WATER_MODEL
        return identifier if identifier in forcefield.water_models else None
    if model not in forcefield.solvent_models:
        if solvent_type in forcefield.solvent_models:
            model = solvent_type
        else:
            model = None
    return model


def get_solvent_model(identifier, forcefield: Union[Forcefield, str] = OPLS_AA):
    """Return a :class:`GromacsSolventModel` corresponding to identifier *identifier*.

    If identifier is "water" then the :data:`DEFAULT_WATER_MODEL` is assumed.
    """
    if isinstance(forcefield, str):
        forcefield = _get_forcefield(forcefield)

    if identifier == "water":
        identifier = DEFAULT_WATER_MODEL
    try:
        return forcefield.water_models[identifier]
    except KeyError:
        try:
            return forcefield.solvent_models[identifier]
        except KeyError:
            msg = "No solvent model with name {0} is available for forcefield {1}.".format(identifier, forcefield.name)
            logger.critical(msg)
            raise ValueError(msg)

def get_top_template(identifier):
    """Return the topology file template suitable for the solvent model."""

    templates = {
        "water": "system.top",
        "octanol": "system.top",
        "cyclohexane": "system.top",
        "wetoctanol": "system_octwet.top",
        "toluene": "system.top",
    }
    try:
        return templates[identifier]
    except KeyError:
        msg = "No template for solvent {0} is available".format(identifier)
        logger.critical(msg)
        raise ValueError(msg)
