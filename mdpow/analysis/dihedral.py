# MDPOW: dihedral.py
# 2021 Alia Lescoulie

from typing import List

import pandas as pd
import numpy as np

import MDAnalysis as mda
from MDAnalysis.exceptions import SelectionError
from MDAnalysis.analysis.dihedrals import calc_dihedrals

from .ensemble import EnsembleAtomGroup, EnsembleAnalysis

import logging

logger = logging.getLogger('mdpow.analysis.dihedral')

# for directory_paths
import os
import re

# for automated_dihedral_analysis
import pathlib

from rdkit import Chem

import seaborn as sns
import matplotlib.pyplot as plt

#import mdpow
#from mdpow.analysis.dihedral import DihedralAnalysis
# fix the way DihedralAnalysis is called

from MDAnalysis.topology.guessers import guess_atom_element

class DihedralAnalysis(EnsembleAnalysis):
    """Analyzes dihedral angles of selections from a single
    :class:`~mdpow.analysis.ensemble.Ensemble` .

    :keywords:

    *dihedral_groups*
        list of :class:`~mdpow.analysis.ensemble.EnsembleAtomGroup`
        with four atoms selected on each. All selections must be from the same
        :class:`~mdpow.analysis.ensemble.Ensemble` .

    Data is returned in a :class:`pandas.DataFrame` with observations sorted by
    selection, solvent, interaction, lambda, time.

    .. ruberic:: Example

    Typical Workflow::

        ens = Ensemble(dirname='Mol')

        dihedral1 = Ens.select_atoms('name C1 or name C2 or name C3 or name C4')
        dihedral2 = Ens.select_atoms('name C5 or name C8 or name C10 or name C12')

        dih_run = DihedralAnalysis([dihedral1, dihedral2]).run(start=0, stop=10, step=1)

    """

    def __init__(self, dihedral_groups: List[EnsembleAtomGroup]):
        self.check_groups_from_common_ensemble(dihedral_groups)
        self.check_dihedral_inputs(dihedral_groups)
        super(DihedralAnalysis, self).__init__(dihedral_groups[0].ensemble)
        self.g1, self.g2, self.g3, self.g4, self.names = self._reorg_groups(dihedral_groups)

    @staticmethod
    def _reorg_groups(groups: List[EnsembleAtomGroup]):
        ag1 = []
        ag2 = []
        ag3 = []
        ag4 = []
        ag_keys = []
        names = []
        for group in groups:
            ag1 += [mda.AtomGroup([ag[0]]) for ag in [group[k] for k in group.keys()]]
            ag2 += [mda.AtomGroup([ag[1]]) for ag in [group[k] for k in group.keys()]]
            ag3 += [mda.AtomGroup([ag[2]]) for ag in [group[k] for k in group.keys()]]
            ag4 += [mda.AtomGroup([ag[3]]) for ag in [group[k] for k in group.keys()]]
            names.append('-'.join([ag1[-1].atoms[0].name, ag2[-1].atoms[0].name,
                                   ag3[-1].atoms[0].name, ag4[-1].atoms[0].name]))
            for k in group.keys():
                ag_keys.append((names[-1], k[0], k[1], k[2]))
        eag1 = EnsembleAtomGroup({ag_keys[i]: ag1[i] for i in range(len(ag_keys))}, groups[0].ensemble)
        eag2 = EnsembleAtomGroup({ag_keys[i]: ag2[i] for i in range(len(ag_keys))}, groups[0].ensemble)
        eag3 = EnsembleAtomGroup({ag_keys[i]: ag3[i] for i in range(len(ag_keys))}, groups[0].ensemble)
        eag4 = EnsembleAtomGroup({ag_keys[i]: ag4[i] for i in range(len(ag_keys))}, groups[0].ensemble)
        return eag1, eag2, eag3, eag4, names

    @staticmethod
    def check_dihedral_inputs(selections):
        for group in selections:
            for k in group.keys():
                if len(group[k]) != 4:
                    msg = ''''Dihedral calculations require AtomGroups with
                              only 4 atoms, %s selected''' % len(group)
                    logger.error(msg)
                    raise SelectionError(msg)

    def _prepare_ensemble(self):
        self._col = ['selection', 'solvent', 'interaction',
                     'lambda', 'time', 'dihedral']
        self.results = pd.DataFrame(columns=self._col)
        self._res_dict = {key: [] for key in self._col}

    def _single_frame(self):
        key_list = [(n, self._key[0], self._key[1], self._key[2]) for n in self.names]
        cord_dict1 = self.g1.positions(keys=key_list)
        cord_dict2 = self.g2.positions(keys=key_list)
        cord_dict3 = self.g3.positions(keys=key_list)
        cord_dict4 = self.g4.positions(keys=key_list)
        cord1 = np.concatenate(tuple([cord_dict1[k] for k in key_list]))
        cord2 = np.concatenate(tuple([cord_dict2[k] for k in key_list]))
        cord3 = np.concatenate(tuple([cord_dict3[k] for k in key_list]))
        cord4 = np.concatenate(tuple([cord_dict4[k] for k in key_list]))
        angle = calc_dihedrals(cord1, cord2, cord3, cord4,
                               box=self.g1[key_list[0]].dimensions)
        angle = np.rad2deg(angle)
        for i in range(len(self.names)):
            result = list(key_list[i]) + [self._ts.time, angle[i]]
            for j in range(len(self._col)):
                self._res_dict[self._col[j]].append(result[j])

    def _conclude_ensemble(self):
        for k in self._col:
            self.results[k] = self._res_dict[k]

# start automation section
def directory_paths(parent_directory=None, csv=None):
    '''Takes a parent directory containing MDPOW simulation project subdirectories,
       or .csv file containing molecule names, resnames, and
       simulation data directory paths as argument and returns
       a Pandas DataFrame for use with DirectoryIteration.
       
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
           must contain header
           format: molecule,resname,path
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

def directory_iteration(dirpaths, figdir=None, padding=45, width=0.9,
                        solvents=('water','octanol'), interactions=('Coulomb','VDW'),
                        start=None, stop=None, step=None):
    #fix to update how I end up handling arguments

    for row in dirpaths.itertuples():
            molname = row.molecule
            resname = row.resname
            datadir = row.path

            automated_dihedral_analysis(datadir=datadir, figdir=figdir,
                                            molname=molname, resname=resname,
                                            padding=padding, width=width,
                                            solvents=solvents, interactions=interactions,
                                            start=start, stop=stop, step=step)

    return

# automated_dihedral_analysis
# assess redundancies and optimize 
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
                    hue_order=["water", "octanol"], col_order=["Coulomb", "VDW"],
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
