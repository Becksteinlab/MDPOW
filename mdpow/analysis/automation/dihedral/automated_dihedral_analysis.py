# MDPOW: automated_dihedral_analysis.py
# 2022 Cade Duckworth

import numpy as np
import pandas as pd
import os
import pathlib

from rdkit import Chem

import seaborn as sns
import matplotlib.pyplot as plt

import mdpow
from mdpow.analysis.dihedral import DihedralAnalysis

import MDAnalysis as mda
from MDAnalysis.topology.guessers import guess_atom_element

'''Contains an assortment of functions useful for automation
   of DihedralAnalysis workflow. See each function for
   requirements, usage, output, or examples. 
   
   Most can be used as standalone or in combination for desired purposes.
   Complete automation encompassed in automated_dihedral_analysis().
'''

def dihedral_indices(datadir, resname):
    '''Requires an MDPOW project directory and resname 
       as input. Uses the topology and trajectory from
       a location that should be present in all MDPOW
       projects and creates an MDAnalysis Universe object.
       Uses resname to select the solute for conversion 
       to RDKit, and will add Hydrogen if not listed
       in the topology. Uses a SMARTS string to obtain all 
       relevant dihedral atom groups and returns their indices
       for use by dihedral_groups and dihedral_groups_ensemble.
       
       :arguments:
       
       *datadir : string*
           path to the location of MDPOW project data
           
       *resname : string*
           resname for the molecule as defined in 
           the topology and trajectory
    '''

    path = pathlib.Path(datadir)
    topology = path / 'FEP/water/Coulomb/0000' / 'md.tpr'
    trajectory = path / 'FEP/water/Coulomb/0000' / 'md.xtc'
    u = mda.Universe(str(topology), str(trajectory))

    try:
        solute = u.select_atoms(f'resname {resname}')
        mol = solute.convert_to('RDKIT')
    except AttributeError:
        u_aug = add_hydrogens(u)
        solute = u_aug.select_atoms(f'resname {resname}')
        mol = solute.convert_to('RDKIT')

    pattern = Chem.MolFromSmarts('[!#1]~[!$(*#*)&!D1]-!@[!$(*#*)&!D1]~[!#1]')
    
    '''SMARTS string that identifies relevant dihedral atom groups
    
    [!#1] : any atom, not Hydrogen
    ~  : any bond
    [!$(*#*)&!D1] : any atom that is not part of linear triple bond
    and not atom with 1 explicit bond
    -!@ : single bond that is not ring bond

    [!$(*#*)&!D1]-!@[!$(*#*)&!D1] : the central portion selects two atoms
    that are not involved in a triple bond and are not terminal,
    that are connected by a single, non-ring bond
    
    [!#1] : the first and last portion specify any bond,
    to any atom that is not hydrogen
    '''
    
    bonds = mol.GetSubstructMatches(pattern)
    return bonds

def dihedral_groups(datadir, resname):
    '''Requires an MDPOW project directory and resname 
       as input. Expands upon usage of dihedral_indices()
       to return an array of the names of each atom within
       its respective dihedral atom group as identified by
       the SMARTS selection string. Not necessary for automation,
       useful for verifying all atom groups of interest are
       properly identified before analysis.
       
       :arguments:
       
       *datadir : string*
           path to the location of MDPOW project data
           
       *resname : string*
           resname for the molecule as defined in 
           the topology and trajectory
    '''

    path = pathlib.Path(datadir)
    topology = path / 'FEP/water/Coulomb/0000' / 'md.tpr'
    trajectory = path / 'FEP/water/Coulomb/0000' / 'md.xtc'

    u = mda.Universe(str(topology), str(trajectory))
    solute = u.select_atoms(f'resname {resname}')

    bonds = dihedral_indices(datadir, resname)
    dihedral_groups = [solute.atoms[list(dihedral_indices)].names for dihedral_indices in bonds]

    return dihedral_groups

def add_hydrogens(u): 
    '''Used by dihedral_indices() for proper conversion
       to RDKit. Adds Hydrogen if not listed in the topology.
    '''

    elements = [guess_atom_element(name) for name in u.atoms.names]
    u.add_TopologyAttr("elements", elements)

    return u

def dihedral_groups_ensemble(bonds, datadir,
                             solvents=('water','octanol'),
                             interactions=('Coulomb','VDW'),
                             start=None, stop=None, step=None): 
    '''Creates one ensemble for the MDPOW project and runs
       DihedralAnalysis for each dihedral atom group identified
       by the SMARTS selection string. For keyword arguments
       see automated_dihedral_analysis() or DihedralAnalysis.
       
       :arguments:
       
       *bonds*
           bond indices for dihedral atom groups
           see dihedral_indices()
           
       *datadir : string*
           path to the location of MDPOW project data
    '''
    
    dih_ens = mdpow.analysis.ensemble.Ensemble(dirname=datadir,
                                               solvents=solvents,
                                               interactions=interactions)

    all_dihedrals = [dih_ens.select_atoms(f'index {b[0]}',
                                          f'index {b[1]}',
                                          f'index {b[2]}',
                                          f'index {b[3]}' ) for b in bonds]

    da = DihedralAnalysis(all_dihedrals)
    da.run(start=start, stop=stop, step=step)
    df = da.results

    return df

def save_df(df, df_save_dir=None, resname=None, molname=None):
    '''Takes a Pandas DataFrame of results from DihedralAnalysis
       as input before padding the angles to optionaly save the raw
       data before plotting dihedral angle frequencies as KDE violins
       with dihedral_violins(). Given a parent directory, creates subdirectory
       for molecule, saves fully sampled csv and per dihedral atom group selection.
       
       :arguments:
       
       *df : Pandas DataFrame*
           results DataFrame from DihedralAnalysis
           
       :keywords:
       
       *df_save_dir : string*
           path to parent directory to create
           subdirectory for saving the csv files
           
    '''
    
    if molname is None:
        molname = resname

    if df_save_dir is not None:
        subdir = molname
        newdir = os.path.join(df_save_dir, subdir)
        os.mkdir(newdir)
        
    df = df.sort_values(by=["selection",
                            "solvent",
                            "interaction",
                            "lambda"]).reset_index(drop=True)
    
    section = df.groupby(by='selection')

    if df_save_dir is not None:
        df.to_csv(f'{newdir}/{molname}_full_df.csv', index=False)
        for name in section:
            sliced_df = df.loc[df['selection'] == name[0]]
            sliced_df.to_csv(f'{newdir}/{molname}_{name[0]}_df.csv', index=False)

    return

def periodic_angle(df, padding=45):
    '''Takes a Pandas DataFrame of results from DihedralAnalysis
       as input and pads the angles to maintain periodicity
       for properly plotting dihedral angle frequencies as KDE violins
       with dihedral_violins(). Creates two new DataFrames based on the 
       cutoff value specified, adds to the angle values, concatenates
       all three DataFrames, maintaining original data and adding padding,
       and returns new augmented DataFrame.
       
       :arguments:
       
       *df : Pandas DataFrame*
           results DataFrame from DihedralAnalysis
           
       :keywords:
       
       *padding : int*
           value must be in degrees
           default set to 45
           
       Example:
       
           da = DihedralAnalysis(all_dihedrals)
           da.run(start=0, stop=100, step=10)
           df = da.results
           df_aug = periodic_angle(df, padding=45)
           plot = dihedral_violins(df_aug, width=0.9)
    '''

    df1 = df[df.dihedral > 180 - padding]
    df1.dihedral -= 360
    df2 = df[df.dihedral < -180 + padding]
    df2.dihedral += 360
    df_aug = pd.concat([df1, df, df2]).reset_index()

    return df_aug

def dihedral_violins(df, width=0.9):
    '''Plots distributions of dihedral angles for 
       one dihedral atom group as violin plots,
       using as input the augmented DataFrame from periodic_angle().
       Returns a figure in the form of seaborn.FacetGrid.
       See periodic_angle() for example.

       :arguments:
       
       *df : Pandas DataFrame*
           augmented DataFrame from periodic_angle()
           
       :keywords:
       
       *width : float*
           width of the violin element (>1 overlaps)
    '''
    
    df['lambda'] = df['lambda'].astype('float') / 1000
    df = df.sort_values(by=["selection",
                            "solvent",
                            "interaction",
                            "lambda"]).reset_index(drop=True)

    # number of Coul windows + 1 / number of VDW windows
    # (+1 for additional space with axes)
    width_ratios = [len(df[df['interaction'] == "Coulomb"]["lambda"].unique()) + 1,
                    len(df[df['interaction'] == "VDW"]["lambda"].unique())]

    g = sns.catplot(data=df, x="lambda", y="dihedral", hue="solvent", col="interaction",
                    kind="violin", split=True, width=width, inner=None, cut=0,
                    linewidth=0.5,
                    hue_order=[df.solvent.iloc[0], df.solvent.iloc[-1]], col_order=["Coulomb", "VDW"],
                    sharex=False, sharey=True,
                    height=2, aspect=2.5,
                    facet_kws={'ylim': (-180, 180),
                               'gridspec_kws': {'width_ratios': width_ratios,
                                                # 'wspace': 0.03
                                                }})
    g.set_xlabels(r"$\lambda$")
    g.set_ylabels(r"dihedral angle $\phi$")
    g.despine(offset=5)

    axC = g.axes_dict['Coulomb']
    axC.yaxis.set_major_locator(plt.matplotlib.ticker.MultipleLocator(60))
    axC.yaxis.set_minor_locator(plt.matplotlib.ticker.MultipleLocator(30))
    axC.yaxis.set_major_formatter(plt.matplotlib.ticker.FormatStrFormatter(r"$%g^\circ$"))

    axV = g.axes_dict['VDW']
    axV.yaxis.set_visible(False)
    axV.spines["left"].set_visible(False)

    return g

def plot_violins(df, resname, figdir=None, molname=None, width=0.9):
    '''Coordinates plotting and optionally saving figures
       for all dihedral atom groups. Makes a subdirectory within
       the specified figure directory using rename or molname provided.
       See automated_dihedral_analysis(), dihedral_violins()
       Returns figures in the form of seaborn.FacetGrid.
    '''

    if molname is None:
        molname = resname

    if figdir is not None:
        subdir = molname
        newdir = os.path.join(figdir, subdir)
        os.mkdir(newdir)

    section = df.groupby(by='selection')

    for name in section:
        plot = dihedral_violins(name[1], width=width)
        plot.set_titles(f'{molname},{name[0]}, ''{col_name}')

        if figdir is not None:
            figfile = pathlib.Path(newdir) / f"{molname}_{name[0]}_violins.pdf"
            plot.savefig(figfile)

    return plot

def automated_dihedral_analysis(datadir=None, df_save_dir=None, figdir=None,
                                resname=None, molname=None,
                                dataframe=None, padding=45, width=0.9,
                                solvents=('water','octanol'),
                                interactions=('Coulomb','VDW'),
                                start=None, stop=None, step=None):
    '''For one MDPOW project, automatically determines all relevant
       dihedral atom groups, conducts DihedralAnalysis for each group,
       pads the dihedral angles from analysis results for all groups,
       creates KDE violin plots of dihedral angle frequencies for each group,
       and optionally saves figure for each group as pdf in provided figure directory.
       Returns KDE violin plot figure for each atom group as seaborn.FacetGrid.
       
       :keywords:
       
       *datadir : string*
           path to the location of MDPOW project data
           
       *df_save_dir : string*
           path to the location to save results DataFrame
           
       *figdir : string*
           optional, path to an existing directory where plots
           can be saved for each dihedral atom group
        
       *resname : string*
           resname for the molecule as defined in 
           the topology and trajectory
           
       *molname : string*
           molecule name to be used for labelling
           plots, if different from resname
           
       *dataframe : Pandas DataFrame*
           optional, if DihedralAnalysis was done prior,
           then results DataFrame can be input to utilize
           angle padding and violin plotting functionality
        
       *padding : int*
           must be in degrees, values for angle padding function
           used for KDE violin plots of dihedral angle frequencies
        
       *width : float*
           used for violin plots
           width of violins, (>1 overlaps)
           see automated_dihedral_analysis.dihedral_violins()
        
       *solvents : tuple*
           solvent specifications for use by DihedralAnalysis
           see DihedralAnalysis
        
       *interactions : tuple*
           interaction specifications for use by DihedralAnalysis
           see DihedralAnalysis
        
       *start,stop,step : int*
           run _single_frame() analysis parameters for use by DihedralAnalysis
           see DihedralAnalysis
       
       Example:
       
           import automated_dihedral_analysis as ada
           
           ada.automated_dihedral_analysis(datadir='/foo/bar/MDPOW_project_data',
                                           figdir='/foo/bar/MDPOW_figure_directory',
                                           resname='UNK', molname='benzene',
                                           padding=45, width=0.9,
                                           solvents=('water','octanol'),
                                           interactions=('Coulomb','VDW'),
                                           start=0, stop=100, step=10)
    '''

    bonds = dihedral_indices(datadir=datadir, resname=resname)

    if dataframe is None:
        df = dihedral_groups_ensemble(bonds, datadir=datadir,
                                      solvents=solvents, interactions=interactions,
                                      start=start, stop=stop, step=step)
    else:
        df = dataframe

    if df_save_dir is not None:
        save_df(df, df_save_dir=df_save_dir, resname=resname, molname=molname)
        
    df_aug = periodic_angle(df, padding=padding)

    return plot_violins(df_aug, resname=resname, figdir=figdir, molname=molname, width=width)          
