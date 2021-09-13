# MDPOW: solvation.py
# 2021 Alia Lescoulie

from typing import List

import pandas as pd

from .ensemble import EnsembleAnalysis, Ensemble

import logging

logger = logging.getLogger('mdpow.dihedral')


class SolvationAnalysis(EnsembleAnalysis):
    """Measures the number of solvent molecules withing the given distances
    in an :class:`~mdpow.analysis.ensemble.Ensemble` .

    :keyword:

    *ensemble*
        The :class:`~mdpow.analysis.ensemble.Ensemble` used in the analysis.

    *distances*
        The cutoff distances around the solute measured in Angstroms

    The data is returned in a :class:`pandas.DataFrame` with observations sorted by
    distance, solvent, interaction, lambda, time.

    .. ruberic:: Example

    Typical Workflow::

        ens = Ensemble(dirname='Mol')

        solv = SolvationAnalysis(ens, [1.2, 2.4]).run(start=0, stop=10, step=1)

    """
    def __init__(self, ensemble: Ensemble, distances: List[float]):
        super(SolvationAnalysis, self).__init__(ensemble)
        self._ens = ensemble
        self._dists = distances

    def _prepare_ensemble(self):
        self._col = ['distance', 'solvent', 'interaction',
                     'lambda', 'time', 'quantity']
        self.results = pd.DataFrame(columns=self._col)
        self._res_dict = {key: [] for key in self._col}

    def _single_frame(self):
        for d in self._dists:
            solvs = len(self._system.select_atoms(f'around {d} not resname SOL'))
            result = [d, self._key[0], self._key[1], self._key[2], self._ts.time, solvs]
            for i in range(len(result)):
                self._res_dict[self._col[i]].append(result[i])

    def _conclude_ensemble(self):
        for k in self._col:
            self.results[k] = self._res_dict[k]
