# MDPOW: registry.py
# 2023 Cade Duckworth

"""
:mod:`mdpow.workflows.registry` --- Registry of currently supported automated workflows
=======================================================================================

:mod:`mdpow.workflows.registry` hosts a dictionary with keys that correspond to an
:class:`~mdpow.analysis.ensemble.EnsembleAnalysis` for which exists a corresponding automated workflow.

.. table:: Currently supported automated workflows.
   :widths: auto
   :name: workflows_registry

   +-------------------------------+------------------------------------------------------------------------------------------------------+
   | key/keyword: EnsembleAnalysis | value: <workflow module>.<top-level automated analysis function>                                     |
   +===============================+======================================================================================================+
   | DihedralAnalysis              | :any:`dihedrals.automated_dihedral_analysis <mdpow.workflows.dihedrals.automated_dihedral_analysis>` |
   +-------------------------------+------------------------------------------------------------------------------------------------------+

.. autodata:: registry

.. seealso:: :mod:`~mdpow.workflows.base`

"""

# import analysis
from mdpow.workflows import dihedrals

registry = {

    'DihedralAnalysis' : dihedrals.automated_dihedral_analysis

}

"""
Each entry corresponds to an :class:`~mdpow.analysis.ensemble.EnsembleAnalysis`
for which exists a corresponding automated workflow.

Intended for use with :mod:`mdpow.workflows.base` to specify which
:class:`~mdpow.analysis.ensemble.EnsembleAnalysis` should run iteratively over
the provided project data directory. To include a new automated workflow for use with
:mod:`mdpow.workflows.base`, create a key that is the name of the corresponding
:class:`~mdpow.analysis.ensemble.EnsembleAnalysis`, with the value defined as
`<workflow module>.<top-level automated analysis function>`.

The currently available workflows are listed in the
:any:`Currently supported automated workflows. <workflows_registry>`

"""