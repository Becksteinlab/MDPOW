# MDPOW: base.py
# 2022 Cade Duckworth

"""
:mod:`mdpow.workflows.base` --- Automated workflow base functions
=================================================================
   
Iterative Functions
-------------------

.. autofunction:: directory_paths
.. autofunction:: directory_iteration
"""

import os
import re
import pandas as pd

#import will need to change based on new naming convention for ada
#and once the function changes for more general use
import mdpow.workflows.dihedrals as ada

import logging

logger = logging.getLogger('mdpow.workflows.base')

def directory_paths(parent_directory=None, csv=None):
    """Takes a parent directory containing MDPOW simulation project subdirectories,
       or .csv file containing :code:`molname`, :code:`resname`, and
       simulation data directory paths as argument and returns a
       :class:`pandas.DataFrame` for use with
       :func:`~mdpow.analysis.workflows.base.directory_iteration`, which iterates
       :func:`~mdpow.analysis.workflows.dihedrals.automated_dihedral_analysis`
       over the project directories included in the
       :func:`~mdpow.analysis.workflows.base.directory_paths` :class:`pandas.DataFrame`.
       
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

#needs to be changed for use with other automated workflows
def directory_iteration(directory_paths, df_save_dir=None, figdir=None, padding=45, width=0.9,
                        solvents=('water','octanol'), interactions=('Coulomb','VDW'),
                        start=None, stop=None, step=None,
                        SMARTS='[!#1]~[!$(*#*)&!D1]-!@[!$(*#*)&!D1]~[!#1]'):
    """Takes a :class:`pandas.DataFrame` created by
       :func:`~mdpow.analysis.workflows.base.directory_paths`
       as input and iterates over the provided projects to implement
       :func:`~mdpow.analysis.workflows.dihedrals.automated_dihedral_analysis`
       for each project directory. Optionally accepts a figure directory for
       saving plots. Extracts :code:`molname`, :code:`resname`, and :code:`dirname`
       from :func:`~mdpow.analysis.workflows.base.directory_paths` :class:`pandas.DataFrame`
       for use in obtaining dihedral groups and plotting dihedral angle frequency KDEs.

       :keywords:

       *figdir*
           optional, path to an existing directory where plots
           can be saved for each dihedral atom group, for each project

       *df_save_dir*
           path to the location to save results :class:`pandas.DataFrame`

       *padding*
           must be in degrees, values for
           :func:`~mdpow.analysis.workflows.dihedrals.periodic_angle`
           used for KDE violin plots of dihedral angle frequencies

       *width*
           used for violin plots
           width of violins, (>1 overlaps)
           see :func:`~mdpow.analysis.workflows.dihedrals.dihedral_violins`

       *solvents*
           Solvents from directory given to the new instance. Default
           :code:`solvents=('water', 'octanol')`

       *interactions*
           Interactions from directory given to the instance. Default
           :code:`interactions=('Coulomb', 'VDW')`

       *SMARTS*
           optional user input of different SMARTS string selection, for
           default see :func:`~mdpow.analysis.workflows.dihedrals.dihedral_indices`

       .. rubric:: Examples

       directory_paths = directory_paths(parent_directory='/foo/bar/MDPOW_projects')
       directory_iteration(directory_paths, figdir='/foo/bar/figure_directory,
       padding=45, width=0.9,
       solvents=('water','octanol'), interactions=('Coulomb','VDW'),
       start=0, stop=100, step=10)
    """

    for row in directory_paths.itertuples():
            molname = row.molecule
            resname = row.resname
            dirname = row.path

            ada.automated_dihedral_analysis(dirname=dirname, df_save_dir=df_save_dir, figdir=figdir,
                                            molname=molname, resname=resname, SMARTS=SMARTS,
                                            padding=padding, width=width,
                                            solvents=solvents, interactions=interactions,
                                            start=start, stop=stop, step=step)

    return
