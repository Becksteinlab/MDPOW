# MDPOW: base.py
# 2022 Cade Duckworth

import os
import re
import pandas as pd

#import will need to change based on new naming convention for ada
#and once the function changes for more general use
import automated_dihedral_analysis as ada

def directory_paths(parent_directory=None, csv=None):
    '''Takes a parent directory containing MDPOW simulation project subdirectories,
       or .csv file containing molecule names, resnames, and
       simulation data directory paths as argument and returns
       a Pandas DataFrame for use with directory_iteration,
       which iterates automated analysis over the project directories
       included in the directory_paths DataFrame.
       
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
           
       Examples:
       
           directory_paths = directory_paths(parent_directory='/foo/bar/MDPOW_projects')
           directory_iteration(directory_paths)
           
           or
           
           directory_paths = directory_paths(csv='/foo/bar/MDPOW.csv')
           directory_iteration(directory_paths)
           
    '''

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
                        start=None, stop=None, step=None):
    '''Takes a Pandas DataFrame created by directory_paths as input
    and iterates over the provided projects to implement automated_dihedral_analysis
    for each project directory. Optionally accepts a figure directory for
    saving plots, and all automated_dihedral_analysis and DihedralAnalysis
    keyword arguments. Extracts molname, resname, and datadir from 
    directory_paths Pandas DataFrame for use in obtaining dihedral groups
    and plotting dihedral angle frequency KDEs.
    
    :keywords:
    
    *figdir : string*
        optional, path to existing directory where plots
        can be saved for each MDPOW project analyzed
        
    *df_save_dir : string*
           path to the location to save results DataFrame
        
    *padding : int*
        must be in degrees, values for angle padding function
        used for KDE violin plots of dihedral angle frequencies
        
    *width : float*
        used for violin plots
        width of violins, (>1 overlaps)
        see automated_dihedral_analysis.dihedral_violins
        
    *solvents : tuple*
        solvent specifications for use with DihedralAnalysis
        see automated_dihedral_analysis
        
    *interactions : tuple*
        interaction specifications for use with DihedralAnalysis
        see automated_dihedral_analysis
        
    *start,stop,step : int*
        run frame analysis parameters for use with DihedralAnalysis
        see automated_dihedral_analysis
        
    Example:
       
           directory_paths = directory_paths(parent_directory='/foo/bar/MDPOW_projects')
           directory_iteration(directory_paths, figdir='/foo/bar/figure_directory,
                               padding=45, width=0.9,
                               solvents=('water','octanol'), interactions=('Coulomb','VDW'),
                               start=0, stop=100, step=10)
           
    '''

    for row in dirpaths.itertuples():
            molname = row.molecule
            resname = row.resname
            datadir = row.path

            ada.automated_dihedral_analysis(datadir=datadir, df_save_dir=df_save_dir, figdir=figdir,
                                            molname=molname, resname=resname,
                                            padding=padding, width=width,
                                            solvents=solvents, interactions=interactions,
                                            start=start, stop=stop, step=step)

    return