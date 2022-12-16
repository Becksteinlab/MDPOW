# MDPOW: hbond.py
# 2022 Cade Duckworth

#NEED TO FIX THE WAY IT DETERMINES TIMES FOR RESULTS
#Review and update old docs

import pandas as pd
import numpy as np

import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA

from .ensemble import Ensemble, EnsembleAtomGroup, EnsembleAnalysis

import logging

logger = logging.getLogger('mdpow.analysis.hbond')

class HBondAnalysis(EnsembleAnalysis):
    """Analyzes potential hydrogen bonds of solute from a single
    :class:`~mdpow.analysis.ensemble.Ensemble` .
    :keywords:
    *solute*
        single :class:`~mdpow.analysis.ensemble.EnsembleAtomGroup`
        selected by resname provided for Universe. Selections must be from the same
        :class:`~mdpow.analysis.ensemble.Ensemble` .
    Data is returned in a :class:`pandas.DataFrame` with observations sorted by
        solvent, interaction, lambda, frame, donor index, 
        hydrogen index, acceptor index, distance, angle.
    .. ruberic:: Example
    Typical Workflow::
        ens = Ensemble(dirname='Mol', solvents=('water', 'octanol'), 
                                      interactions=('Coulomb', 'VDW'))
        solute = ens.select_atoms('resname UNK')
        hb = HBondAnalysis(solute, acceptors_sel='resname UNK', d_a_cutoff=2.5)
        hb.run(start=0, stop=1000, step=10)
        
        HydrogenBondAnalysis
        numpy array results format:
        <frame>, <donor index (0-based)>,
        <hydrogen index (0-based)>,
        <acceptor index (0-based)>,
        <distance>, <angle>
        
        HydrogenBondAnalysis
        :keywords:
        acceptors_sel=None
        donors_sel=None
        hydrogens_sel=None
        d_h_cutoff=1.2
        d_a_cutoff=3.0
        d_h_a_angle_cutoff=150
        update_selections=True
        
        https://docs.mdanalysis.org/stable/documentation_pages/analysis/hydrogenbonds.html

        ***NEED TO REFINE THIS AND ADD MORE DETAIL***
    """

    def __init__(self, solute: EnsembleAtomGroup, **kwargs):
        super(HBondAnalysis, self).__init__(solute.ensemble)
        self._solute = solute
        self.hb_kwargs = kwargs
        

    def _prepare_ensemble(self):
        self._col = ['solvent', 'interaction',
                     'lambda', 'time', 'frame', 
                     'donor index', 'hydrogen index', 
                     'acceptor index', 'distance', 'angle']
        self.results = pd.DataFrame(columns=self._col)
        self._res_dict = {key: [] for key in self._col}

    def _single_universe(self):
        solute = self._solute[self._key]
        hb = HBA(universe=self._solute.ensemble[self._key], **self.hb_kwargs)
        hb.run(verbose=True, start=self.start, stop=self.stop, step=self.step)
        results = hb.results.hbonds
        #self._times = np.array([self._trajectory[ts].time for ts in self._trajectory[self.start:self.stop:self.step]])

        
        for h in results:
            result = [self._key[0], self._key[1],self._key[2],
                      (self._trajectory.dt*h[0]), int(h[0]), int(h[1]),
                      int(h[2]), int(h[3]), h[4], h[5]]
            for i in range(len(self._col)):
                self._res_dict[self._col[i]].append(result[i])

    def _conclude_ensemble(self):
        for k in self._col:
            self.results[k] = self._res_dict[k]         
