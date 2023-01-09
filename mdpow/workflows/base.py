# MDPOW: base.py
# 2022 Cade Duckworth

"""
:mod:`mdpow.workflows.base` --- Automated workflow base functions
=================================================================

.. autofunction:: directory_paths
.. autofunction:: directory_iteration
"""

import os
import re
import pandas as pd

#import dihedrals as ada
import solvations as asa

import logging

logger = logging.getLogger('mdpow.workflows.base')

def directory_paths(parent_directory=None, csv=None):
    """Takes a parent directory containing MDPOW project subdirectories,
       or .csv file containing :code:`molname`, :code:`resname`, and
       simulation data directory paths as argument and returns a
       :class:`pandas.DataFrame` for use with
       :func:`~mdpow.workflows.base.directory_iteration`, which iterates
       :func:`~mdpow.workflows.dihedrals.automated_dihedral_analysis`
       over the project directories included in the
       :func:`~mdpow.workflows.base.directory_paths` :class:`pandas.DataFrame`.
       
       :keywords:
       
       *parent_directory*
           the path for the location of the top directory 
           under which the subdirectories of MDPOW simulation
           data exist, additionally creates a 'dir.csv' file
           for user manipulation of metadata and future reference
           
       *csv*
           .csv file containing the molecule names, resnames,
           and paths, in that order, for the MDPOW simulation
           data to be iterated over
           must contain header of the form:
           format: molecule,resname,path

       .. rubric:: Examples

       directory_paths = directory_paths(parent_directory='/foo/bar/MDPOW_projects')
       directory_iteration(directory_paths)

       or

       directory_paths = directory_paths(csv='/foo/bar/MDPOW.csv')
       directory_iteration(directory_paths)
    """

    if parent_directory is not None:

        locations = []

        reg_compile = re.compile('FEP')              
        for dirpath, dirnames, filenames in os.walk(parent_directory):
            result = [dirpath.strip() for dirname in dirnames if  reg_compile.match(dirname)]
            if result:
                locations.append(result[0])

        resnames = []

        for loc in locations:
            res_temp = loc.strip().split('/')
            resnames.append(res_temp[-1])

        directory_paths = pd.DataFrame(
        {
            "molecule": resnames,
            "resname": resnames,
            "path": locations
        }
    )
        directory_paths.to_csv('dir.csv', index=False)

    elif csv is not None:
        locations = pd.read_csv(csv)
        directory_paths = locations.sort_values(by=['molecule', 'resname', 'path']).reset_index(drop=True)

    return directory_paths

def directory_iteration(directory_paths, df_save_dir=None, figdir=None,
                        solvents=('water','octanol'), interactions=('Coulomb','VDW'),
                        ensemble_analysis=None, SMARTS=None, padding=None, width=None,
                        start=None, stop=None, step=None, distances=None):
    """Takes a :class:`pandas.DataFrame` created by
       :func:`~mdpow.workflows.base.directory_paths`
       as input and iterates over the provided projects to implement
       :func:`~mdpow.workflows.dihedrals.automated_dihedral_analysis`
       for each project directory. Optionally accepts a figure directory for
       saving plots. Extracts :code:`molname`, :code:`resname`, and :code:`dirname`
       from :func:`~mdpow.workflows.base.directory_paths` :class:`pandas.DataFrame`
       for use in obtaining dihedral groups and plotting dihedral angle frequency KDEs.

       :keywords:
       
       *ensemble_analysis*
           name of the :class:`~mdpow.analysis.ensemble.EnsembleAnalysis`
           that corresponds to the desired automation script

       *figdir*
           optional, path to an existing directory where plots
           can be saved for each dihedral atom group, for each project

       *df_save_dir*
           path to the location to save results :class:`pandas.DataFrame`

       *padding*
           must be in degrees, values for
           :func:`~mdpow.workflows.dihedrals.periodic_angle`
           used for KDE violin plots of dihedral angle frequencies
           recommended default = 45

       *width*
           used for violin plots
           width of violins, (>1 overlaps)
           see :func:`~mdpow.workflows.dihedrals.dihedral_violins`
           recommended default = 0.9

       *solvents*
           Solvents from directory given to the new instance
           default :code:`solvents=('water', 'octanol')`

       *interactions*
           Interactions from directory given to the instance
           default :code:`interactions=('Coulomb', 'VDW')`

       *SMARTS*
           optional user input of different SMARTS string selection 
           recommended default = '[!#1]~[!$(*#*)&!D1]-!@[!$(*#*)&!D1]~[!#1]'
           
       *distances*
           Array like of the cutoff distances around the solute measured in Angstroms.

       .. rubric:: Examples

       directory_paths = directory_paths(parent_directory='/foo/bar/MDPOW_projects')
       directory_iteration(directory_paths, figdir='/foo/bar/figure_directory,
       padding=45, width=0.9,
       solvents=('water','octanol'), interactions=('Coulomb','VDW'),
       start=0, stop=100, step=10)
    """

    analyses = {
        'DihedralAnalysis' : 0,
        'SolvationAnalysis' : asa.automated_solvation_analysis,
        'HBondAnalysis' : 0
    }
    
    if ensemble_analysis is not None:
        try:
            for row in directory_paths.itertuples():
                molname = row.molecule
                resname = row.resname
                dirname = row.path
                #analyses[ensemble_analysis](dirname, solvents=solvents,
                #                            interactions=interactions, **kwargs)
                
                analyses[ensemble_analysis](dirname=dirname, df_save_dir=df_save_dir, figdir=figdir,
                                            molname=molname, resname=resname, SMARTS=SMARTS,
                                            padding=padding, width=width, distances=distances,
                                            solvents=solvents, interactions=interactions,
                                            start=start, stop=stop, step=step)
                
        except NotImplementedError as err:
            logger.error("invalid EnsembleAnalysis selection", err, msg)
            msg = ('An EnsembleAnalysis type that corresponds to an existing '
                   'automated workflow module must be input as a kwarg. '
                   'ex. ensemble_analysis=DihedralAnalysis')
            raise NotImplementedError(msg)
    else:
        msg = ('Expected keyword argument for ensemble_analysis is '
               'missing. An EnsembleAnalysis type that corresponds '
               'to an existing automated workflow module must be '
               'input as a kwarg. ex. ensemble_analysis=DihedralAnalysis')
        raise ValueError(msg)

    return
