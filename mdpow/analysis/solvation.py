# MDPOW: solvation.py
# 2021 Alia Lescoulie

from typing import List

import numpy as np
import pandas as pd

import MDAnalysis as mda
from MDAnalysis.lib.distances import capped_distance

from .ensemble import EnsembleAnalysis, Ensemble, EnsembleAtomGroup

import logging

logger = logging.getLogger('mdpow.dihedral')


class SolvationAnalysis(EnsembleAnalysis):
    """Measures the number of solvent molecules withing the given distances
    in an :class:`~mdpow.analysis.ensemble.Ensemble` .

    :Parameters:

    *solute*
        An :class:`~mdpow.analysis.ensemble.EnsembleAtom` containing the solute
        used to measure distance.

    *solvent*
        An :class:`~mdpow.analysis.ensemble.EnsembleAtom` containing the solvents
        counted in by the distance measurement. Each solvent atom is counted by the
        distance calculation.

    *distances*
        Array like of the cutoff distances around the solute measured in Angstroms.

    The data is returned in a :class:`pandas.DataFrame` with observations sorted by
    distance, solvent, interaction, lambda, time.

    .. rubric:: Example

    Typical Workflow::

        ens = Ensemble(dirname='Mol')
        solvent = ens.select_atoms('resname SOL and name OW')
        solute = ens.select_atoms('resname UNK')

        solv_dist = SolvationAnalysis(solute, solvent, [1.2, 2.4]).run(stop=10)

    """
    def __init__(self, solute: EnsembleAtomGroup, solvent: EnsembleAtomGroup, distances: List[float]):
        self.check_groups_from_common_ensemble([solute, solvent])
        super(SolvationAnalysis, self).__init__(solute.ensemble)
        self._solute = solute
        self._solvent = solvent
        self._dists = distances

    def _prepare_ensemble(self):
        self._col = ['distance', 'solvent', 'interaction',
                     'lambda', 'time', 'N_solvent']
        self.results = pd.DataFrame(columns=self._col)
        self._res_dict = {key: [] for key in self._col}

    def _single_frame(self):
        solute = self._solute[self._key]
        solvent = self._solvent[self._key]
        pairs, distances = capped_distance(solute.positions, solvent.positions,
                                          max(self._dists), box=self._ts.dimensions)
        solute_i, solvent_j = np.transpose(pairs)
        for d in self._dists:
            close_solv_atoms = solvent[solvent_j[distances < d]]
            result = [d, self._key[0], self._key[1],self._key[2],
                      self._ts.time, close_solv_atoms.n_atoms]
            for i in range(len(self._col)):
                self._res_dict[self._col[i]].append(result[i])

    def _conclude_ensemble(self):
        for k in self._col:
            self.results[k] = self._res_dict[k]
