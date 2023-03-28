# MDPOW: registry.py
# 2023 Cade Duckworth

"""
:mod:`mdpow.workflows.registry` --- Registry of currently supported automated workflows
=======================================================================================

Each entry in :mod:`mdpow.workflows.registry` corresponds to an :class:`~mdpow.analysis.ensemble.EnsembleAnalysis`
for which exists a corresponding automated workflow.

Currently supported workflows :class:`~mdpow.analysis.dihedral.DihedralAnalysis`

.. property:: `~mdpow.workflows.registry.registry`

Intended for use with :mod:`mdpow.workflows.base` to specify which
:class:`~mdpow.analysis.ensemble.EnsembleAnalysis` should run iteratively over
the provided project data.

.. seealso:: :mod:`~mdpow.workflows.base`, :mod:`~mdpow.workflows.dihedrals`

"""

# import analysis
from mdpow.workflows import dihedrals
 
registry = {
    'DihedralAnalysis' : dihedrals.automated_dihedral_analysis
}



    
    
    
    
    