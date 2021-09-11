# MDPOW: dihedral.py
# 2021 Alia Lescoulie

from typing import List

import pandas as pd

from MDAnalysis.exceptions import SelectionError
from MDAnalysis.analysis.dihedrals import calc_dihedrals

from .ensemble import EnsembleAtomGroup, EnsembleAnalysis

import logging
logger = logging.getLogger('mdpow.analysis.dihedral')


class DihedralAnalysis(EnsembleAnalysis):
    """Analyzes dihedral angles of selections from a single
    :class:`~mdpow.analysis.ensemble.Ensemble` .

    :keywords:

    *dihedral_groups
        list of :class:`~mdpow.analysis.ensemble.EnsembleAtomGroup`
        with four atoms selected on each. All selections must be from the same
        :class:`~mdpow.analysis.ensemble.Ensemble` .

    Data is returned in a :class:`pandas.DataFrame` with observations sorted by
    selection, solvent, interaction, lambda, time.
    """

    def __init__(self, dihedral_groups: List[EnsembleAtomGroup]):
        self.check_groups_from_common_ensemble(dihedral_groups)
        self.check_dihedral_inputs(dihedral_groups)
        super(DihedralAnalysis, self).__init__(dihedral_groups[0].ensemble)
        self._sel = dihedral_groups

    @staticmethod
    def check_dihedral_inputs(selections):
        for group in selections:
            for k in group.keys():
                if len(group[k]) != 4:
                    msg = ''''Dihedral calculations require AtomGroups with
                              only 4 atoms, %s selected''' % len(group)
                    logger.error(msg)
                    raise SelectionError(msg)

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
