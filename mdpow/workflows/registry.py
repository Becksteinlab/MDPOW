# MDPOW: registry.py
# 2023 Cade Duckworth

"""
:mod:`mdpow.workflows.registry` --- Registry of currently supported automated workflows
=======================================================================================

.. autodata:: registry

Each entry in :mod:`mdpow.workflows.registry` corresponds to an :class:`~mdpow.analysis.ensemble.EnsembleAnalysis`
for which exists a corresponding automated workflow.

Currently supported workflows - :mod:`~mdpow.workflows.dihedrals`

Intended for use with :mod:`mdpow.workflows.base` to specify which
:class:`~mdpow.analysis.ensemble.EnsembleAnalysis` should run iteratively over
the provided project data directory. To include a new automated workflow for use with
:mod:`mdpow.workflows.base`, create a key that is the name of the corresponding
:class:`~mdpow.analysis.ensemble.EnsembleAnalysis`, with the value defined as
`<module name>.<top-level function>`.

.. seealso:: :mod:`~mdpow.workflows.base`, :mod:`~mdpow.workflows.dihedrals`

"""

# import analysis
from mdpow.workflows import dihedrals
 
registry = {
    'DihedralAnalysis' : dihedrals.automated_dihedral_analysis
}
