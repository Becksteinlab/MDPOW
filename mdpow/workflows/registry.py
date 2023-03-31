# MDPOW: registry.py
# 2023 Cade Duckworth

"""
:mod:`mdpow.workflows.registry` --- Registry of currently supported automated workflows
=======================================================================================

Each entry in :mod:`mdpow.workflows.registry` corresponds to an :class:`~mdpow.analysis.ensemble.EnsembleAnalysis`
for which exists a corresponding automated workflow.

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

"""
.. data:: mdpow.workflows.registry.registry

   :type: dictionary

   Intended for use with :mod:`mdpow.workflows.base` to specify which
   :class:`~mdpow.analysis.ensemble.EnsembleAnalysis` should run iteratively over
   the provided project data directory. To include a new automated workflow for use with
   :mod:`mdpow.workflows.base`, create a key that is the name of the corresponding
   :class:`~mdpow.analysis.ensemble.EnsembleAnalysis`, with the value defined as
   `<workflow module>.<top-level automated analysis function>`.

   The currently available workflows are listed in the
   :ref:`Currently supported automated workflows. <mdpow.workflows.registry.workflows_registry>`

"""
registry = {
    'DihedralAnalysis' : dihedrals.automated_dihedral_analysis
}
