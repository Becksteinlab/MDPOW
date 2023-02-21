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

from mdpow.workflows import registry

import logging

logger = logging.getLogger('mdpow.workflows.base')

# DONT FORGET UPDATES FROM SAMPL9 REPO

def directory_paths(parent_directory=None, csv=None, csv_save_dir=None):
    """Takes a top directory containing MDPOW projects and determines
       the molname, resname, and path, of each MDPOW project within.

       Optionally takes a .csv file containing `molname`, `resname`, and
       `paths`, in that order. 

       :keywords:

       *parent_directory*
           the path for the location of the top directory 
           under which the subdirectories of MDPOW simulation
           data exist, additionally creates a 'dir_paths.csv' file
           for user manipulation of metadata and for future reference

       *csv*
           .csv file containing the molecule names, resnames,
           and paths, in that order, for the MDPOW simulation
           data to be iterated over
           must contain header of the form:
           format: molecule,resname,path

       *csv_save_dir*
           optionally provided directory to save .csv file, otherwise,
           data will be saved in current working directory

       :returns:

       *directory_paths*
           :class:`pandas.DataFrame` containing MDPOW project metadata

       .. rubric:: Example
       
       Typical Workflow::

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
        if csv_save_dir is not None:
            
            directory_paths.to_csv(f'{csv_save_dir}/dir_paths.csv', index=False)
            logger.info(f'dir_paths saved under {csv_save_dir}')

        else:
            current_directory = os.system('pwd')
            directory_paths.to_csv('dir_paths.csv', index=False)
            logger.info(f'dir_paths saved under {current_directory}')

    elif csv is not None:
        locations = pd.read_csv(csv)
        directory_paths = locations.sort_values(by=['molecule', 'resname', 'path']).reset_index(drop=True)

    return directory_paths

def directory_iteration(directory_paths, ensemble_analysis, **kwargs):
    """Takes a :class:`pandas.DataFrame` created by
       :func:`~mdpow.workflows.base.directory_paths`
       and iterates over the project paths to implement
       :func:`~mdpow.workflows.dihedrals.automated_dihedral_analysis`
       for each project directory.

       Compatibility with more automated analyses in development.

       :keywords:

       *directory_paths*
           :class:`pandas.DataFrame` that provides paths to MDPOW projects

       *ensemble_analysis*
           name of the :class:`~mdpow.analysis.ensemble.EnsembleAnalysis`
           that corresponds to the desired automated workflow module

       *kwargs*
           keyword arguments from :func:`~mdpow.workflows.dihedrals.automated_dihedral_analysis`

           .. autodata:: mdpow.workflows.dihedrals.automated_dihedral_analysis

       .. rubric:: Example

       Typical Workflow::

           directory_paths = directory_paths(parent_directory='/foo/bar/MDPOW_projects')
           directory_iteration(directory_paths, ensemble_analysis='DihedralAnalysis', **kwargs)

    """

    try:
        for row in directory_paths.itertuples():
            molname = row.molecule
            resname = row.resname
            dirname = row.path

            logger.info(f'starting {molname}')

            registry.registry[ensemble_analysis](dirname=dirname, resname=resname, molname=molname, **kwargs)

            logger.info(f'{molname} completed')

    except KeyError as err:
        msg = (f'Invalid ensemble_analysis {err}. An EnsembleAnalysis type that corresponds to an existing '
                'automated workflow module must be input as a kwarg. '
                'ex: ensemble_analysis=\'DihedralAnalysis\'')
        logger.error(f'{err} is an invalid selection')

        raise KeyError(msg)

    except TypeError as err:
        msg = (f'Invalid ensemble_analysis {ensemble_analysis}. An EnsembleAnalysis type that corresponds to an existing '
                'automated workflow module must be input as a kwarg. '
                'ex: ensemble_analysis=\'DihedralAnalysis\'')
        logger.error(f'workflow module for {ensemble_analysis} does not exist yet')

        raise TypeError(msg)

    return logger.info('all analyses completed')
