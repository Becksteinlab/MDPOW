# MDPOW: dihedral.py
# 2021 Alia Lescoulie

import numpy as np
import pandas as pd

from MDAnalysis.analysis.dihedrals import calc_dihedrals


from .ensemble import Ensemble, EnsembleAtomGroup, EnsembleAnalysis

import logging
logger = logging.getLogger('mdpow.analysis.dihedral')


class DihedralAnalysis(EnsembleAnalysis):
    """Analyzes dihedral angles of a simulation.

    Accepts an :class:`~mdpow.analysis.ensemble.EnsembleAtomGroup`
    with four atoms selected on each.
    """

    def __init__(self, DihedralGroups: EnsembleAtomGroup):
        super(DihedralAnalysis, self).__init__(DihedralGroups.ensemble())
        self._sel = DihedralGroups

    def _prepare_ensemble(self):
        self._col = ['solvent', 'interaction',
                     'lambda', 'time', 'dihedral']
        self.results = pd.DataFrame(columns=self._col)
        self._res_dict = {key: [] for key in self._col}

    def _single_frame(self):
        angle = calc_dihedrals(self._sel[self._key].positions[0], self._sel[self._key].positions[1],
                               self._sel[self._key].positions[2], self._sel[self._key].positions[3])
        angle_data = [self._key[0], self._key[1], self._key[2], self._ts.time, angle]
        for i in range(len(self._col)):
            self._res_dict[self._col[i]].append(angle_data[i])

    def _conclude_ensemble(self):
        for k in self._col:
            self.results[k] = self._res_dict[k]
