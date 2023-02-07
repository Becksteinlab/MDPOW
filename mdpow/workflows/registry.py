# MDPOW: registry.py
# 2023 Cade Duckworth

# import analysis
from mdpow.workflows import dihedrals_217 as dihedrals
 
registry = {
    'DihedralAnalysis' : dihedrals.automated_dihedral_analysis,
    'SolvationAnalysis' : None,
    'HBondAnalysis' : None
}



    
    
    
    
    