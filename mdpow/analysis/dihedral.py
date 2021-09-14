# MDPOW: dihedral.py
# 2021 Alia Lescoulie

from typing import List

import pandas as pd
import numpy as np

import MDAnalysis as mda
from MDAnalysis.exceptions import SelectionError
from MDAnalysis.analysis.dihedrals import calc_dihedrals

from .ensemble import EnsembleAtomGroup, EnsembleAnalysis

import logging

logger = logging.getLogger('mdpow.analysis.dihedral')


class DihedralAnalysis(EnsembleAnalysis):
    """Analyzes dihedral angles of selections from a single
    :class:`~mdpow.analysis.ensemble.Ensemble` .

    :keywords:

    *dihedral_groups*
        list of :class:`~mdpow.analysis.ensemble.EnsembleAtomGroup`
        with four atoms selected on each. All selections must be from the same
        :class:`~mdpow.analysis.ensemble.Ensemble` .

    Data is returned in a :class:`pandas.DataFrame` with observations sorted by
    selection, solvent, interaction, lambda, time.

    .. ruberic:: Example

    Typical Workflow::

        ens = Ensemble(dirname='Mol')

        dihedral1 = Ens.select_atoms('name C1 or name C2 or name C3 or name C4')
        dihedral2 = Ens.select_atoms('name C5 or name C8 or name C10 or name C12')

        dih_run = DihedralAnalysis([dihedral1, dihedral2]).run(start=0, stop=10, step=1)

    """

    def __init__(self, dihedral_groups: List[EnsembleAtomGroup]):
        self.check_groups_from_common_ensemble(dihedral_groups)
        self.check_dihedral_inputs(dihedral_groups)
        super(DihedralAnalysis, self).__init__(dihedral_groups[0].ensemble)
        self.g1, self.g2, self.g3, self.g4, self.names = self._reorg_groups(dihedral_groups)

    @staticmethod
    def _reorg_groups(groups: List[EnsembleAtomGroup]):
        ag1 = []
        ag2 = []
        ag3 = []
        ag4 = []
        ag_keys = []
        names = []
        for group in groups:
            ag1 += [mda.AtomGroup([ag[0]]) for ag in [group[k] for k in group.keys()]]
            ag2 += [mda.AtomGroup([ag[1]]) for ag in [group[k] for k in group.keys()]]
            ag3 += [mda.AtomGroup([ag[2]]) for ag in [group[k] for k in group.keys()]]
            ag4 += [mda.AtomGroup([ag[3]]) for ag in [group[k] for k in group.keys()]]
            names.append('-'.join([ag1[-1].atoms[0].name, ag2[-1].atoms[0].name,
                                   ag3[-1].atoms[0].name, ag4[-1].atoms[0].name]))
            for k in group.keys():
                ag_keys.append((names[-1], k[0], k[1], k[2]))
        eag1 = EnsembleAtomGroup({ag_keys[i]: ag1[i] for i in range(len(ag_keys))}, groups[0].ensemble)
        eag2 = EnsembleAtomGroup({ag_keys[i]: ag2[i] for i in range(len(ag_keys))}, groups[0].ensemble)
        eag3 = EnsembleAtomGroup({ag_keys[i]: ag3[i] for i in range(len(ag_keys))}, groups[0].ensemble)
        eag4 = EnsembleAtomGroup({ag_keys[i]: ag4[i] for i in range(len(ag_keys))}, groups[0].ensemble)
        return eag1, eag2, eag3, eag4, names

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
        key_list = [(n, self._key[0], self._key[1], self._key[2]) for n in self.names]
        cord_dict1 = self.g1.positions(keys=key_list)
        cord_dict2 = self.g2.positions(keys=key_list)
        cord_dict3 = self.g3.positions(keys=key_list)
        cord_dict4 = self.g4.positions(keys=key_list)
        cord1 = np.concatenate(tuple([cord_dict1[k] for k in key_list]))
        cord2 = np.concatenate(tuple([cord_dict2[k] for k in key_list]))
        cord3 = np.concatenate(tuple([cord_dict3[k] for k in key_list]))
        cord4 = np.concatenate(tuple([cord_dict4[k] for k in key_list]))
        angle = calc_dihedrals(cord1, cord2, cord3, cord4,
                               box=self.g1[key_list[0]].dimensions)
        angle = np.rad2deg(angle)
        for i in range(len(self.names)):
            result = list(key_list[i]) + [self._ts.time, angle[i]]
            for j in range(len(self._col)):
                self._res_dict[self._col[j]].append(result[j])

    def _conclude_ensemble(self):
        for k in self._col:
            self.results[k] = self._res_dict[k]
