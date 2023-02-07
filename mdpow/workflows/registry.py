# MDPOW: registry.py
# 2023 Cade Duckworth

# new workflows registry module
# holds dictionary of available automated analysis workflows
# for use with workflows base module, directory iteration function

# holds kwarg dictionaries for each automated analyses

import .dihedrals
# import analysis type and respective kwarg dict when used in workflows base module
# use tuple assignment, analysis, _ = registry.DihedralAnalysis, etc
registry = {
    'DihedralAnalysis' : (dihedrals.automated_dihedral_analysis,
                          dihedrals_kwargs = {
                              'dirname' : None,
                              'df_save_dir' : None,
                              'figdir' : None,
                              'resname' : None,
                              'molname' : None,
                              'SMARTS' : SMARTS_DEFAULT,
                              'solvents' : SOLVENTS_DEFAULT,
                              'interactions' : INTERACTIONS_DEFAULT,
                              'dataframe' : None,
                              'padding' : 45,
                              'width' : 0.9,
                              'start' : None,
                              'stop' : None,
                              'step' : None
                          }
                         ),
    'SolvationAnalysis' : None,
    'HBondAnalysis' : None
}


    