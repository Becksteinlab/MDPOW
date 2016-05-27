# POW package __init__.py
# Copyright (c) 2012 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""
Tables of hard-coded values used in mdpow
=========================================

TODO: Move these data into files in the top directory.
"""

from __future__ import absolute_import

from numkit.observables import QuantityWithError

import logging
logger = logging.getLogger("mdpow.tables")

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

def _create_water_models(watermodelsdat):
    models = {}
    for line in GMX_WATERMODELS_DAT.split('\n'):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        fields = line.split()
        identifier, name = fields[:2]
        description = " ".join(fields[2:])
        models[identifier] = GromacsWaterModel(identifier, name=name, description=description)
    return models

#: Use the default water model unless another water model is chosen in the runinput
#" file (``setup.watermodel``).
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

#------------------------------------------------------------
# constants
#------------------------------------------------------------

#: Avogadro's constant from http://physics.nist.gov/cgi-bin/cuu/Value?na
#: in mol\ :sup:`-1`.
N_Avogadro = 6.02214129e23 # mol^-1

#: Boltzmann's constant *k* from http://physics.nist.gov/cgi-bin/cuu/Value?k
#: in J/K.
k_Boltzmann = 1.3806488e-23

#: Molecular weight in g/mol from PubChem for `water`_ and `1-octanol`_.
#: .. _water:       http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid=962
#: .. _`1-octanol`: http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid=957
molecular_weight = {
    "water":     18.015280,
    "octanol":   130.227920,
    "1-octanol": 130.227920,
    }

#: :program:`make_ndx` commands to select either the OW in TIP4P
#: water (key *water*) or OH in 1-octanol (*octanol*). This is rather
#: specific for mdpow because we assume that the residue name is
#: "SOL" for water or "1OCT" for octanol.
#: These commands are used in :func:`count_solvent_molecules`.
solvent_selections = {
    "water": "r SOL & a OW",
    "octanol": "r 1OCT & a OH",
    }


#: preliminary table for solvent density 1t T=300 K and P=1 bar; needs
#: to be made more general, possibly by using a spline interpolation
#: for the equation of state. Unit kg/m^3.
solvent_density = {
    "water": QuantityWithError(992.342, 0.18),     # TIP4P
    # xfer3-336/Equilibrium/water/MD_NPT
    # 992.342 (0.18), 652 w, 15ns (disc 5ns from 20ns), T=296.009 (0.019), P=0.993607 (0.024) bar
    "octanol": QuantityWithError(819.491, 0.24),
    # octanol_box_compresibility/Equilibrium/octanol/MD_NPT
    # 819.491 (0.24), 496 octanols, 303.5010 (0.0038) K, 1.00669 (0.0026) bar, 95ns (disc.5ns of 100ns)

    # experimental values
    # http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid=957
    "octanol_exp": 827.0,  #  at 20 deg C/4 deg C
    # Matsuo 1989 (see kappaT); tabulated for T=298.15 K ... 348.15 K
    "octanol_Matsuo1989_298K": 822.8, # T=298.15 K and P=1 atm
    "octanol_Matsuo1989_303K": 819.3, # T=303.15 K and P=1 atm

    # other results (not used at the moment)
    # --------------------------------------
    "water_296": QuantityWithError(993.379, 0.13),     # TIP4P,
    # xfer3-302/Equilibrium/water/MD_NPT
    # 993.379 (0.13), 296 w, 15ns (discarded 5ns from 20ns)
    "octanol_123": QuantityWithError(818.654, 0.39),   # OPLS-AA
    # xfer3-336/Equilibrium/octanol/MD_NPT
    # 818.654 (0.39), 123 octanols
    }

# In the future, use the table to compute a proper interpolating function
# density...
_density = """
Table[water_density]: TIP4P water density from simulations
======== ========== ============= ============ =============
solvent  forcefield temperature   pressure     density
======== ========== ============= ============ =============
water    TIP4P            300.0          1.0      993.379
octanol  OPLS-AA          300.0          1.0      200
======== ========== ============= ============ =============
"""

#: volume of a solvent molecule at standard conditions, in nm^3
solvent_molecular_volume = {
    "water": QuantityWithError(0.0301475, 5.52147e-06),    # TIP4P sim with 652 waters
    "octanol": QuantityWithError(0.263893, 7.66129e-05),   # 1-oct with 496 oct
    #--------------------------------------------------
    "water_296": QuantityWithError(0.0301174, 4.05405e-06),# TIP4 sim with 296 w???
}


#: 1 bar in kJ nm\ :sup:`-3`
bar_kJnm = 1e5 * 1e-3 * (1e9)**(-3)
#: 1 bar in kJ mol\ :sup:`-1` nm\ :sup:`-3`
bar = bar_kJnm * N_Avogadro
#: 1 atm in kJ mol\ :sup:`-1` nm\ :sup:`-3`
atm = 1.01325 * bar
#: 1 MPa in kJ mol\ :sup:`-1` nm\ :sup:`-3`
MPa = 0.1 * bar
#: isothermal compressibility coefficient (in 1/bar)
kappaT = {
    'water':
        {
        #: isothermal compressibility of water at 25C and 1 atm [Fine et al, J
        #: Chem Phys 59 (1973)], in 1/bar
        'exp': 4.5248e-5,
        #: TIP4P compressibility at 25C and 1atm [Mahoney & Jorgensen, JChemPhys 112 2000]
        'Mahoney2000': 6.0e-5 / (atm/bar),
        #: TIP4P compressibility at 10C and 1atm [Matubayasi & Levy, J Phys Chem 100 (1996)]
        #: in 1/bar
        'Matubayasi1996': 5.00e-5 / (atm/bar),
        },
    'octanol':
        {
        # ISOTHERMAL COMPRESSIBILITY: 6.82 @ 1 ATM & 0 DEG C
        # http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid=957
        # [Weast, R.C. (ed.) Handbook of Chemistry and Physics. 69th
        # ed. Boca Raton, FL: CRC Press Inc., 1988-1989., p. F-14]
        # **PEER REVIEWED**
        'CRC': 6.82e-4 / MPa, # UNITS?? -- multiplied by 1e-4 from comparison to Matsuo1989 (OB)
        # S. Matsuo and T. Makita. Volumetric properties of 1-alkanols
        # at temperatures in the range 298--348 K and pressures up to
        # 40 mpa. International Journal of Thermophysics, 10:885--897,
        # 1989.  10.1007/BF00514483. URL
        # http://dx.doi.org/10.1007/BF00514483.
        'Matsuo1989': 7.61e-4 / MPa,  # 7.61 in 1e4 MPa^-1 @ 298 K and 1 atm

        # 1/MPa = 10 * 1/bar --> 7.61e-4 / MPa = 7.61e-3 / bar;
        # i.e. 1-octanol is 2 orders of magnitude more compressible
        # than water and hence the NVT/NPT correction will be about
        # 100 times smaller than for water
        },
}
kappaT['water']['DEFAULT'] = kappaT['water']['Mahoney2000']
kappaT['octanol']['DEFAULT'] = kappaT['octanol']['Matsuo1989']


