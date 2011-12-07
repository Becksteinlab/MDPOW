"""
Tables of hard-coded values used in mdpow.

TODO: Move these data into files in the top directory.
"""
from numkit.observables import QuantityWithError

#: Avogadro's constant from http://physics.nist.gov/cgi-bin/cuu/Value?na
#: in mol<sup>-1</sup>.
N_Avogadro = 6.02214129e23 # mol^-1

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
    "water": QuantityWithError(992.342, 0.18)     # TIP4P
    # xfer3-336/Equilibrium/water/MD_NPT
    # 992.342 (0.18), 652 w, 15ns (disc 5ns from 20ns), T=296.009 (0.019), P=0.993607 (0.024) bar
    "octanol": QuantityWithError(819.491, 0.24),
    # octanol_box_compresibility/Equilibrium/octanol/MD_NPT
    # 819.491 (0.24), 496 octanols, 303.5010 (0.0038) K, 1.00669 (0.0026) bar

    # other results (not used at the moment)
    # --------------------------------------
    "water_296": QuantityWithError(993.379, 0.13)     # TIP4P,
    # xfer3-302/Equilibrium/water/MD_NPT
    # 993.379 (0.13), 296 w, 15ns (discarded 5ns from 20ns)
    "octanol_123": QuantityWithError(818.654, 0.39)   # OPLS-AA
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

solvent_molecular_volume = {
    "water": 0.0301174,  #  0.0301174 (4.05405e-06) same TIP4P sim as above
    "octanol": None,
}


bar_kJnm = 1e5 * 1e-3 * (1e9)**(-3)  # in kJ*nm**-3
bar = bar_kJnm * N_Avogadro          # in kJ mol**-1 nm**-3
atm = 1.01325 * bar
kappaT = {
    #: isothermal compressibility of water at 25C and 1 atm [Fine et al, J
    #: Chem Phys 59 (1973)], in 1/bar
    'exp': 4.5248e-5,
    #: TIP4P compressibility at 25C and 1atm [Mahoney & Jorgensen, JChemPhys 112 2000]
    'Mahoney2000': 6.0e-5 / (atm/bar),
    #: TIP4P compressibility at 10C and 1atm [Matubayasi & Levy, J Phys Chem 100 (1996)]
    #: in 1/bar
    'Matubayasi1996': 5.00e-5 / (atm/bar),
    }
kappaT['DEFAULT'] = kappaT['Mahoney2000']



