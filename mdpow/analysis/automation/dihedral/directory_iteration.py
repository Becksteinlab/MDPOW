# MDPOW: directory_iteration.py
# 2022 Cade Duckworth

import pandas

import automated_dihedral_analysis as ada



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
