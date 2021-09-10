# MDPOW: dihedral.py
# 2021 Alia Lescoulie

from typing import List

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

    def __init__(self, dihedralgroups: List[EnsembleAtomGroup]):
        super(DihedralAnalysis, self).__init__(dihedralgroups[0].ensemble)
        self._sel = dihedralgroups

    def _prepare_ensemble(self):
        self._col = ['selection', 'solvent', 'interaction',
                     'lambda', 'time', 'dihedral']
        self.results = pd.DataFrame(columns=self._col)
        self._res_dict = {key: [] for key in self._col}

    def _single_frame(self):
        for group in self._sel:
            angle = calc_dihedrals(group[self._key].positions[0], group[self._key].positions[1],
                                   group[self._key].positions[2], group[self._key].positions[3])
            name = ''.join([x.name for x in group[self._key].atoms])
            angle_data = [name, self._key[0], self._key[1],
                          self._key[2], self._ts.time, angle]
            for i in range(len(self._col)):
                self._res_dict[self._col[i]].append(angle_data[i])

    def _conclude_ensemble(self):
        for k in self._col:
            self.results[k] = self._res_dict[k]
