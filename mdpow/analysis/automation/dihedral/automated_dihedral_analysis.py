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



def dihedral_indices(datadir, resname): 

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
    '''
    [!#1] : any atom, not Hydrogen
    ~  : any bond
    !$(*#*)&!D1] : any atom that is not part of linear triple bond
    and not atom with 1 explicit bond
    -!@ : single bond that is not ring bond (or is it not aromatic?)

    [!$(*#*)&!D1]-!@[!$(*#*)&!D1] : the central portion selects two atoms
    that are not involved in a triple bond and are not terminal,
    that are connected by a single, non-ring bond
    
    [!#1] : the first and last portion specify any bond,
    to any atom that is not hydrogen
    '''
    #^verify and update this docstring
    
    bonds = mol.GetSubstructMatches(pattern)
    return bonds

def dihedral_groups(datadir, resname):

    #uses a standard, expected file location to
    #determine structure and obtain dihedral atom groups
    path = pathlib.Path(datadir)
    topology = path / 'FEP/water/Coulomb/0000' / 'md.tpr'
    trajectory = path / 'FEP/water/Coulomb/0000' / 'md.xtc'

    u = mda.Universe(str(topology), str(trajectory))
    solute = u.select_atoms(f'resname {resname}')

    bonds = dihedral_indices(datadir, resname)
    dihedral_groups = [solute.atoms[list(dihedral_indices)].names for dihedral_indices in bonds]

    return dihedral_groups

def add_hydrogens(u): 

    elements = [guess_atom_element(name) for name in u.atoms.names]
    u.add_TopologyAttr("elements", elements)

    return u

def dihedral_groups_ensemble(bonds, datadir,
                             solvents=('water','octanol'),
                             interactions=('Coulomb','VDW'),
                             start=None, stop=None, step=None): 
    #need to bounce out ensemble kwargs
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

def periodic_angle(df, padding=45): 

    df1 = df[df.dihedral > 180 - padding]
    df1.dihedral -= 360
    df2 = df[df.dihedral < -180 + padding]
    df2.dihedral += 360
    df_aug = pd.concat([df1, df, df2]).reset_index()

    return df_aug

def dihedral_violins(df, width=0.9):
    #need to fix plotting code to match solvents, interactions
    """Plot distributions of all dihedrals as violin plots.

    Parameters
    ----------
    df : DataFrame
    width : float
          width of the violin element (>1 overlaps)

    Returns
    -------
    seaborn.FacetGrid
    """
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

def automated_dihedral_analysis(datadir=None, figdir=None,
                                resname=None, molname=None,
                                dataframe=None, padding=45, width=0.9,
                                solvents=('water','octanol'),
                                interactions=('Coulomb','VDW'),
                                start=None, stop=None, step=None): 

    bonds = dihedral_indices(datadir=datadir, resname=resname)

    if dataframe is None:
        df = dihedral_groups_ensemble(bonds, datadir=datadir,
                                      solvents=solvents, interactions=interactions,
                                      start=start, stop=stop, step=step)
    else:
        df = dataframe

    df_aug = periodic_angle(df, padding=padding)

    return plot_violins(df_aug, resname=resname, figdir=figdir, molname=molname, width=width)          
