# MDPOW: solvations.py
# 2022 Cade Duckworth

"""
:mod:`mdpow.workflows.solvations` --- Automation for :class:`SolvationAnalysis`
==============================================================================
:mod:`~mdpow.workflows.solvations` module with functions
useful for automated use of
:class:`~mdpow.analysis.solvation.SolvationAnalysis`.
See each function for usage, output, and examples. 
Most functions can be used as standalone or in combination
depending on the desired results. Complete automation encompassed in
:func:`~mdpow.workflows.solvations.automated_solvation_analysis`.

.. autofunction:: solvation_ensemble
.. autofunction:: solvation_analysis
.. autofunction:: asa_save_df
.. autofunction:: automated_solvation_analysis
"""

import os
import numpy as np
import pandas as pd

import mdpow
from mdpow.analysis.solvation import SolvationAnalysis

import logging

logger = logging.getLogger('mdpow.workflows.solvations')

def solvation_ensemble(dirname, resname, solvents=('water', 'octanol'),
                       interactions=('Coulomb', 'VDW')):

    ens = mdpow.analysis.ensemble.Ensemble(dirname=dirname,
                                           interactions=interactions,
                                           solvents=solvents)
    solute = ens.select_atoms(f'resname {resname}')
    solvent = ens.select_atoms(f'not resname {resname}')
    return solute, solvent
                              
def solvation_analysis(solute=None, solvent=None, distances=None,
                       start=None, stop=None, step=None):
    
    solv = SolvationAnalysis(solute=solute, solvent=solvent, distances=distances)
    ds = solv.run(start=start, stop=stop, step=step)
    df = solv.results
    return df

def asa_save_df(df, df_save_dir=None, resname=None, molname=None):
    '''Takes a :class:`pandas.DataFrame` of results from 
       :class:`~mdpow.analysis.solvation.SolvationAnalysis`
       as input to optionaly save the data.
       Given a parent directory, creates subdirectory
       for molecule, saves fully sampled csv.
       :keywords:
       *df*
           results :class:`pandas.DataFrame` from
           :class:`~mdpow.analysis.solvation.SolvationAnalysis`
       *df_save_dir*
           path to parent directory to create
           subdirectory for saving the .csv files
    '''

    if molname is None:
        molname = resname

    if df_save_dir is not None:
        subdir = molname
        newdir = os.path.join(df_save_dir, subdir)
        os.mkdir(newdir)

    df = df.sort_values(by=["solvent",
                            "interaction",
                            "lambda"]).reset_index(drop=True)

    if df_save_dir is not None:
        df.to_csv(f'{newdir}/{molname}_full_df.csv', index=False, compression='bz2')
        # this part might need some work
    return
    
def automated_solvation_analysis(dirname, df_save_dir=None, resname=None, molname=None,
                                 solvents=('water', 'octanol'), interactions=('Coulomb', 'VDW'),
                                 distances=[1.2, 2.4], figdir=None,
                                 start=None, stop=None, step=None,
                                 SMARTS=None, padding=None, width=None):
    # figure out kwargs for each analysis type
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
    """
                              
    components = solvation_ensemble(dirname=dirname, resname=resname, solvents=solvents)
    df = solvation_analysis(solute=components[0], solvent=components[1],
                            distances=distances, start=start, stop=stop, step=step)
    if df_save_dir is not None:
         asa_save_df(df, df_save_dir=df_save_dir, resname=resname, molname=molname)
                              
    return df