# MDPOW: solvation.py
# 2021 Alia Lescoulie

from typing import List

import pandas as pd

import MDAnalysis as mda

from .ensemble import EnsembleAnalysis, Ensemble, EnsembleAtomGroup

import logging

logger = logging.getLogger('mdpow.dihedral')


class SolvationAnalysis(EnsembleAnalysis):
    """Measures the number of solvent molecules withing the given distances
    in an :class:`~mdpow.analysis.ensemble.Ensemble` .

    :keyword:

    *ensemble*
        The :class:`~mdpow.analysis.ensemble.Ensemble` used in the analysis.

    *distances*
        The cutoff distances around the solute measured in Angstroms.

    The data is returned in a :class:`pandas.DataFrame` with observations sorted by
    distance, solvent, interaction, lambda, time.

    .. ruberic:: Example

    Typical Workflow::

        ens = Ensemble(dirname='Mol')

        solv = SolvationAnalysis(ens, [1.2, 2.4]).run(start=0, stop=10, step=1)

    """
    def __init__(self, solute: EnsembleAtomGroup, solvent: EnsembleAtomGroup, distances: List[float]):
        self.check_groups_from_common_ensemble([solute, solvent])
        super(SolvationAnalysis, self).__init__(solute.ensemble)
        self._solu = solute
        self._solv = solvent
        self._dists = distances

    def _prepare_ensemble(self):
        self._sel = 'name '
        keys = [k for k in self._solu.keys()]
        for n in self._solu[keys[0]].names:
            self._sel += f' {n}'
        self._col = ['distance', 'solvent', 'interaction',
                     'lambda', 'time', 'quantity']
        self.results = pd.DataFrame(columns=self._col)
        self._res_dict = {key: [] for key in self._col}

    def _single_universe(self):
        self._temp_sys = self._solv[self._key] + self._solu[self._key]

    def _single_frame(self):
        for d in self._dists:
            solvs = len(self._temp_sys.select_atoms(f'around {d} ({self._sel})').residues)
            result = [d, self._key[0], self._key[1], self._key[2], self._ts.time, solvs]
            for i in range(len(result)):
                self._res_dict[self._col[i]].append(result[i])

    def _conclude_ensemble(self):
        for k in self._col:
            self.results[k] = self._res_dict[k]
