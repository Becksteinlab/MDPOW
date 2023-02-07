# MDPOW: dihedrals.py
# 2022 Cade Duckworth

"""
:mod:`mdpow.workflows.dihedrals` --- Automation for :class:`DihedralAnalysis`
==============================================================================

:mod:`~mdpow.workflows.dihedrals` module provides functions for automated
workflows that encompass :class:`~mdpow.analysis.dihedral.DihedralAnalysis`.
See each function for requirements and examples. 

Most functions can be used as standalone, individually, or in combination
depending on the desired results. Details of the completely automated workflow
are discussed under :func:`~mdpow.workflows.dihedrals.automated_dihedral_analysis`.

.. autodata:: solvents_default
    :annotation: = ('water', 'octanol')
.. autodata:: interactions_default
    :annotation: = ('Coulomb', 'VDW')
.. autodata:: SMARTS_DEFAULT
    :annotation: = [!#1]~[!$(*#*)&!D1]-!@[!$(*#*)&!D1]~[!#1]
.. autofunction:: automated_dihedral_analysis
.. autofunction:: build_universe
.. autofunction:: rdkit_conversion
.. autofunction:: dihedral_indices
.. autofunction:: dihedral_groups
.. autofunction:: dihedral_groups_ensemble
.. autofunction:: save_df
.. autofunction:: periodic_angle
.. autofunction:: dihedral_violins
.. autofunction:: plot_violins

"""

import os
import pathlib
import numpy as np
import pandas as pd

from rdkit import Chem

import seaborn as sns
import matplotlib.pyplot as plt

import mdpow
from mdpow.analysis.dihedral import DihedralAnalysis

import MDAnalysis as mda
from MDAnalysis.topology.guessers import guess_atom_element

import logging

logger = logging.getLogger('mdpow.workflows.dihedrals')



SOLVENTS_DEFAULT = ('water', 'octanol')
"""Default solvents are water and octanol.

    * Must match solvents used in project directory.
    * One or two solvents can be specified.
    * Current solvents supported,
      .. seealso:: :mod:`mdpow.forcefields`

"""

INTERACTIONS_DEFAULT = ('Coulomb', 'VDW')
"""Default interactions set to Coulomb and VDW.

    * Default values should not be changed
    * Order should not be changed

"""

SMARTS_DEFAULT = '[!#1]~[!$(*#*)&!D1]-!@[!$(*#*)&!D1]~[!#1]'
"""Default SMARTS string to identify relevant dihedral atom groups.

     * ``[!#1]`` : any atom, not Hydrogen
     * ``~``  : any bond
     * ``[!$(*#*)&!D1]`` : any atom that is not part of linear triple
       bond and not atom with 1 explicit bond
     * ``-!@`` : single bond that is not ring bond
     * ``[!$(*#*)&!D1]-!@[!$(*#*)&!D1]`` : the central portion selects
       two atoms that are not involved in a triple bond and are not terminal,
       that are connected by a single, non-ring bond
     * ``[!#1]~` or `~[!#1]`` : the first and last portion specify any bond,
       to any atom that is not hydrogen

"""

def build_universe(dirname):
    """Builds :class:`~MDAnalysis.core.universe.Universe` from
       ``water/Coulomb/0000`` topology and trajectory for the specified project.
       
       Used by :func:`~mdpow.workflows.dihedrals.rdkit_conversion`
       and :func:`~mdpow.workflows.dihedrals.dihedral_indices` to obtain atom indices
       for each dihedral atom group.
       
       :keywords:

       *dirname*
           Molecule Simulation directory. Loads simulation files present in
           lambda directories into the new instance. With this method for
           generating an :class:`~mdpow.analysis.ensemble.Ensemble` the
           lambda directories are explored and
           :meth:`~mdpow.analysis.ensemble.Ensemble._load_universe_from_dir`
           searches for .gro, .gro.bz2, .gro.gz, and .tpr files for topology,
           and .xtc files for trajectory. It will default to using the tpr file
           available.
           
       :returns:

       *u*
           :class:`~MDAnalysis.core.universe.Universe` object
           
    """

    path = pathlib.Path(dirname)
    topology = path / 'FEP/water/Coulomb/0000' / 'md.tpr'
    trajectory = path / 'FEP/water/Coulomb/0000' / 'md.xtc'
    u = mda.Universe(str(topology), str(trajectory))

    return u

def rdkit_conversion(u, resname):
    """Converts the solute, ``resname``, of the
       :class:`~MDAnalysis.core.universe.Universe` to :class:`RDKit.Mol` object
       for use with a ``SMARTS`` selection string to identify dihedral atom groups.
       
       Accepts :class:`~MDAnalysis.core.universe.Universe` object made with
       :func:`~mdpow.workflows.dihedrals.build_universe` and a ``resname`` as input.
       Uses ``resname`` to select the solute for conversion to :class:`RDKit.Mol`,
       and will add element attributes for Hydrogen if not listed in the topology.
       
       :keywords:
       
       *u*
           :class:`~MDAnalysis.core.universe.Universe` object
           
       *resname*
           ``resname`` for the molecule as defined in 
           the topology and trajectory
           
       :returns:
       
       *tuple(mol, solute)*
           function call returns tuple, see below
       
       *mol*
           :class:`RDKit.Mol` object converted from :code:`solute`
           
       *solute*
           molecule specified by :func:`~MDAnalysis.core.universe.Universe.select_atoms`
           for :class:`~MDAnalysis.core.universe.Universe` object

    """

    try:
        solute = u.select_atoms(f'resname {resname}')
        mol = solute.convert_to('RDKIT')
    except AttributeError:
        elements = [guess_atom_element(name) for name in u.atoms.names]
        u.add_TopologyAttr("elements", elements)
        solute = u.select_atoms(f'resname {resname}')
        mol = solute.convert_to('RDKIT')

    return mol, solute

def dihedral_indices(dirname, resname, SMARTS=SMARTS_DEFAULT):
    '''Uses a `SMARTS` selection string to identify indices for
       relevant dihedral atom groups.
       
       Requires an MDPOW project directory and :code:`resname` 
       as input. With :func:`~mdpow.workflows.dihedrals.build_universe` and
       :func:`~mdpow.workflows.dihedrals.rdkit_conversion`, uses the topology
       and trajectory from :code:`water/Coulomb/0000` and creates a
       :class:`~MDAnalysis.core.universe.Universe` object.
       Uses a :code:`SMARTS` string to obtain all relevant dihedral atom groups.

       :keywords:

       *dirname*
           Molecule Simulation directory. Loads simulation files present in
           lambda directories into the new instance. With this method for
           generating an :class:`~mdpow.analysis.ensemble.Ensemble` the
           lambda directories are explored and
           :meth:`~mdpow.analysis.ensemble.Ensemble._load_universe_from_dir`
           searches for .gro, .gro.bz2, .gro.gz, and .tpr files for topology,
           and .xtc files for trajectory. It will default to using the tpr file
           available.
           
       *resname*
           :code:`resname` for the molecule as defined in 
           the topology and trajectory

       *SMARTS*
           The default SMARTS string is described in detail under :data:`SMARTS_DEFAULT`.
           
       :returns:
       
       *atom_group_indices*
           Tuple of tuples of indices for each dihedral atom group.

    '''

    u = build_universe(dirname=dirname)
    mol = rdkit_conversion(u=u, resname=resname)[0]
    pattern = Chem.MolFromSmarts(SMARTS)
    atom_group_indices = mol.GetSubstructMatches(pattern)
    
    return atom_group_indices

def dihedral_groups(dirname, resname, SMARTS=SMARTS_DEFAULT):
    '''Uses the indices of the relevant dihedral atom groups determined
       by :func:`~mdpow.workflows.dihedral.dihedral_indices`
       and returns the names for each atom in each group.
       
       Requires an MDPOW project directory and :code:`resname` 
       as input. Expands upon usage of
       :func:`~mdpow.workflows.dihedral.dihedral_indices`
       to return an array of the names of each atom within
       its respective dihedral atom group as identified by
       the :code:`SMARTS` selection string. Not necessary
       for automation, useful for verifying all atom groups of interest
       are properly identified before analysis.

       :keywords:

       *dirname*
           Molecule Simulation directory. Loads simulation files present in
           lambda directories into the new instance. With this method for
           generating an :class:`~mdpow.analysis.ensemble.Ensemble` the
           lambda directories are explored and
           :meth:`~mdpow.analysis.ensemble.Ensemble._load_universe_from_dir`
           searches for .gro, .gro.bz2, .gro.gz, and .tpr files for topology,
           and .xtc files for trajectory. It will default to using the tpr file
           available.

       *resname*
           :code:`resname` for the molecule as defined in 
           the topology and trajectory

       *SMARTS*
           The default SMARTS string is described in detail under :data:`SMARTS_DEFAULT`.
           
       :returns:
       
       *dihedral_groups*
           List of :func:`numpy.array` for atom names in each dihedral atom group.

    '''

    # revisit this (reminder 1/30/23)
    # temporary fix for current indexing method (solute.atoms...)
    u = build_universe(dirname=dirname)
    solute = rdkit_conversion(u=u, resname=resname)[1]
    atom_group_indices = dihedral_indices(dirname=dirname, resname=resname, SMARTS=SMARTS)
    dihedral_groups = [solute.atoms[list(dihedral_indices)].names for dihedral_indices
                       in atom_group_indices]

    return dihedral_groups

def dihedral_groups_ensemble(dirname, atom_group_indices,
                             solvents=SOLVENTS_DEFAULT,
                             interactions=INTERACTIONS_DEFAULT,
                             start=None, stop=None, step=None): 
    '''Creates one :class:`~mdpow.analysis.ensemble.Ensemble` for the MDPOW
       project and runs :class:`~mdpow.analysis.dihedral.DihedralAnalysis`
       for each dihedral atom group identified by the :code:`SMARTS`
       selection string.
       
       .. seealso:: :func:`~mdpow.workflows.dihedrals.automated_dihedral_analysis`,
                    :class:`~mdpow.analysis.dihedral.DihedralAnalysis`

       :keywords:

       *dirname*
           Molecule Simulation directory. Loads simulation files present in
           lambda directories into the new instance. With this method for
           generating an :class:`~mdpow.analysis.ensemble.Ensemble` the
           lambda directories are explored and
           :meth:`~mdpow.analysis.ensemble.Ensemble._load_universe_from_dir`
           searches for .gro, .gro.bz2, .gro.gz, and .tpr files for topology,
           and .xtc files for trajectory. It will default to using the tpr file
           available.

       *atom_group_indices*
           tuples of atom indices for dihedral atom groups
           .. seealso:: :func:`~mdpow.workflows.dihedrals.dihedral_indices`

       *solvents*
           The default solvents are documented under :data:`SOLVENTS_DEFAULT`

       *interactions*
           The default interactions are documented under :data:`INTERACTIONS_DEFAULT`

       *start, stop, step*
           .. seealso:: :class:`mdpow.analysis.ensemble.EnsembleAnalysis`

       :returns:
       
       *df*
           :class:`pandas.DataFrame` of :class:`~mdpow.analysis.dihedral.DihedralAnalysis`
           results, including all dihedral atom groups for molecule of current project

    '''

    dih_ens = mdpow.analysis.ensemble.Ensemble(dirname=dirname,
                                               solvents=solvents,
                                               interactions=interactions)
    indices = atom_group_indices
    all_dihedrals = [dih_ens.select_atoms(f'index {i[0]}',
                                          f'index {i[1]}',
                                          f'index {i[2]}',
                                          f'index {i[3]}' ) for i in indices]

    da = DihedralAnalysis(all_dihedrals)
    da.run(start=start, stop=stop, step=step)
    df = da.results

    return df

def save_df(df, df_save_dir, resname=None, molname=None):
    '''Takes a :class:`~pandas.DataFrame` of results from 
       :class:`~mdpow.analysis.dihedral.DihedralAnalysis`
       as input before padding the angles to optionaly save the raw
       data.

       Optionally saves results before padding the angles for periodicity
       and plotting dihedral angle frequencies as KDE violins
       with :func:`~mdpow.workflows.dihedrals.dihedral_violins`.
       Given a parent directory, creates subdirectory for molecule,
       saves fully sampled, unpadded results as compressed csv.

       :keywords:

        *df*
           :class:`pandas.DataFrame` of :class:`~mdpow.analysis.dihedral.DihedralAnalysis`
           results, including all dihedral atom groups for molecule of current project

       *df_save_dir*
           optional, path to the location to save results :ref:`pandas.DataFrame`

       *resname*
           :code:`resname` for the molecule as defined in 
           the topology and trajectory

       *molname*
           molecule name to be used for labelling
           plots, if different from :code:`resname`

       :returns:
       
       *.csv.bz2*
           saves the results :class:`pandas.DataFrame` as a compressed
           csv file, default: .bz2
    '''

    df = df.sort_values(by=["selection",
                            "solvent",
                            "interaction",
                            "lambda"]).reset_index(drop=True)

    if molname is None:
        molname = resname

    subdir = molname
    newdir = os.path.join(df_save_dir, subdir)
    os.mkdir(newdir)

    # time and compress level can be adjusted as kwargs
    df.to_csv(f'{newdir}/{molname}_full_df.csv.bz2',
              index=False, compression='bz2')

    logger.info(f'Results DataFrame saved as '
                 '{newdir}/{molname}_full_df.csv.bz2')

    #return

def periodic_angle(df, padding=45):
    '''Pads the angles from the results :class:`~pandas.DataFrame`
       to maintain periodicity in the violin plots.
    
       Takes a :class:`pandas.DataFrame` of results from
       :class:`~mdpow.analysis.dihedral.DihedralAnalysis`
       as input and pads the angles to maintain periodicity
       for properly plotting dihedral angle frequencies as KDE violins
       with :func:`~mdpow.workflows.dihedrals.dihedral_violins`.
       Creates two new :class:`pandas.DataFrame` based on the 
       cutoff value specified, adds to the angle values, concatenates
       all three :class:`pandas.DataFrame`, maintaining original data and
       adding padding, and returns new augmented :class:`pandas.DataFrame`.

       :keywords:

       *df*
           :class:`pandas.DataFrame` of :class:`~mdpow.analysis.dihedral.DihedralAnalysis`
           results, including all dihedral atom groups for molecule of current project

       *padding*
           value in degrees
           default: 45
           
       :returns:
       
       *df_aug*
           augmented results :class:`pandas.DataFrame` containing
           padded dihedral angles as specified by :code:`padding`

       .. rubric:: Example
       
       Typical Workflow::

           da = DihedralAnalysis(all_dihedrals)
           da.run(start=0, stop=100, step=10)
           df = da.results
           df_aug = periodic_angle(df, padding=45)
           plot = dihedral_violins(df_aug, width=0.9)

    '''

    # To-Do
    # warning for trying to set a value on a copy
    # of a slice of a DataFrame needs to be suppressed
    # is not relevant, double check again
    df1 = df[df.dihedral > 180 - padding]
    df1.dihedral -= 360
    df2 = df[df.dihedral < -180 + padding]
    df2.dihedral += 360
    df_aug = pd.concat([df1, df, df2]).reset_index()

    return df_aug

def dihedral_violins(df, width=0.9, solvents=SOLVENTS_DEFAULT):
    '''Plots distributions of dihedral angles for one dihedral atom group
       as violin plots, using as input the augmented :class:`pandas.DataFrame`
       from :func:`~mdpow.workflows.dihedrals.periodic_angle`.

       :keywords:

       *df*
           augmented results :class:`pandas.DataFrame` from
           :func:`~mdpow.workflows.dihedrals.periodic_angle`

       *width*
           width of the violin element (>1 overlaps)
           default: 0.9

       *solvents*
           The default solvents are documented under :data:`SOLVENTS_DEFAULT`
       
       :returns:

       *violin plot*
           Returns a :class:`seaborn.FacetGrid` object containing a violin plot of the
           kernel density estimations (KDE) of the dihedral angle frequencies for each
           relevant dihedral atom group in the molecule from the current MDPOW project

    '''

    # To-Do
    # make sure value sorting and index dropping/resetting
    # is consistent and makes sense throughout
    df['lambda'] = df['lambda'].astype('float') / 1000
    df = df.sort_values(by=["selection",
                            "solvent",
                            "interaction",
                            "lambda"]).reset_index(drop=True)

    # number of Coul windows + 1 / number of VDW windows
    # (+1 for additional space with axes)
    width_ratios = [len(df[df['interaction'] == "Coulomb"]["lambda"].unique()) + 1,
                    len(df[df['interaction'] == "VDW"]["lambda"].unique())]

    solvs = np.asarray(solvents)
    solv2 = 'octanol'
    if solvs.size > 1:
        solv2 = solvs[1]
    # To-Do
    # ^necessary when plotting only water
    # there is likely a better solution
    # review this before merge (1/30/23)
    # see lines 440-450
    
    g = sns.catplot(data=df, x="lambda", y="dihedral", hue="solvent", col="interaction",
                    kind="violin", split=True, width=width, inner=None, cut=0,
                    linewidth=0.5,
                    hue_order=[solvs[0], solv2], col_order=["Coulomb", "VDW"],
                    #hue_order=["water", "octanol"], col_order=["Coulomb", "VDW"],
                    #hue_order=[df.solvent.iloc[0], df.solvent.iloc[-1]], col_order=["Coulomb", "VDW"],
                    # causes problems and is not a good fix 
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

def plot_violins(df, resname, figdir=None, molname=None, width=0.9, solvents=SOLVENTS_DEFAULT):
    '''Coordinates plotting and optionally saving figures for all dihedral
       atom groups. Makes a subdirectory within the specified
       :code:`figdir` using :code:`resname` or :code:`molname` provided.
       See :func:`~mdpow.workflows.dihedrals.automated_dihedral_analysis`
       or :func:`~mdpow.workflows.dihedrals.dihedral_violins`
       for examples. Returns figures in the form of :class:`seaborn.FacetGrid`.
       
       :keywords:
       
       *df*
           augmented results :class:`pandas.DataFrame` from
           :func:`~mdpow.workflows.dihedrals.periodic_angle`

       *resname*
           :code:`resname` for the molecule as defined in 
           the topology and trajectory

       *figdir*
           optional, path to the location to save figures

       *molname*
           molecule name to be used for labelling
           plots, if different from :code:`resname`

       *width*
           width of the violin element (>1 overlaps)
           default: 0.9
           .. seealso:: :func:`~mdpow.workflows.dihedrals.dihedral_violins`

       *solvents*
           The default solvents are documented under :data:`SOLVENTS_DEFAULT`

       :returns:

       *violin plot*
           Returns a :class:`seaborn.FacetGrid` object containing a violin plot of the
           kernel density estimations (KDE) of the dihedral angle frequencies for each
           relevant dihedral atom group in the molecule from the current MDPOW project

    '''

    if molname is None:
        molname = resname

    if figdir is not None:
        subdir = molname
        newdir = os.path.join(figdir, subdir)
        os.mkdir(newdir)

    section = df.groupby(by='selection')

    for name in section:
        plot = dihedral_violins(name[1], width=width, solvents=solvents)
        # decide if 'set_titles' belongs here or in dihedral_violins
        # appears twice in new modified plots that add Mol image to plots
        plot.set_titles(f'{molname},{name[0]}, ''{col_name}')
        # ^belongs here bc of naming convention for iteration (1/30/23)

        if figdir is not None:
            figfile = pathlib.Path(newdir) / f"{molname}_{name[0]}_violins.pdf"
            plot.savefig(figfile)
            logger.info(f'Figure saved as {figfile}')

    return plot

def automated_dihedral_analysis(dirname=None, df_save_dir=None, figdir=None,
                                resname=None, molname=None,
                                SMARTS=SMARTS_DEFAULT,
                                dataframe=None, padding=45, width=0.9,
                                solvents=SOLVENTS_DEFAULT,
                                interactions=INTERACTIONS_DEFAULT,
                                start=None, stop=None, step=None):
    '''Runs :class:`~mdpow.analysis.dihedral.DihedralAnalysis` for a single MDPOW
       project and creates violin plots of dihedral angle frequencies for each
       relevant dihedral atom group.
    
       For one MDPOW project, automatically determines all relevant dihedral atom groups
       in the molecule, runs :class:`~mdpow.analysis.dihedral.DihedralAnalysis` for each group,
       pads the dihedral angles from analysis results for all groups to maintain periodicity,
       creates violin plots of dihedral angle frequencies for each group, separately, from the
       padded results.

       Optionally saves all pre-padded :class:`mdpow.analysis.dihedral.DihedralAnalysis` results
       as a single :class:`pandas.DataFrame`, and separate violin plots for each dihedral atom group
       in :code:`df_save_dir`, and :code:`figdir` directories provided, respectively.

       :keywords:

       *dirname*
           Molecule Simulation directory. Loads simulation files present in
           lambda directories into the new instance. With this method for
           generating an :class:`~mdpow.analysis.ensemble.Ensemble` the
           lambda directories are explored and
           :meth:`~mdpow.analysis.ensemble.Ensemble._load_universe_from_dir`
           searches for .gro, .gro.bz2, .gro.gz, and .tpr files for topology,
           and .xtc files for trajectory. It will default to using the tpr file
           available.

       *df_save_dir*
           optional, path to the location to save results :ref:`pandas.DataFrame`

       *figdir*
           optional, path to the location to save figures

       *resname*
           :code:`resname` for the molecule as defined in 
           the topology and trajectory

       *molname*
           molecule name to be used for labelling
           plots, if different from :code:`resname`
           
       *SMARTS*
           The default SMARTS string is described in detail under :data:`SMARTS_DEFAULT`.

       *dataframe*
           optional, if :class:`~mdpow.analysis.dihedral.DihedralAnalysis`
           was done prior, then results :class:`pandas.DataFrame` can be
           input to utilize angle padding and violin plotting functionality

       *padding*
           value in degrees
           default: 45
           .. seealso:: :func:`~mdpow.workflows.dihedrals.periodic_angle`

       *width*
           width of the violin element (>1 overlaps)
           default: 0.9
           .. seealso:: :func:`~mdpow.workflows.dihedrals.dihedral_violins`

       *solvents*
           The default solvents are documented under :data:`SOLVENTS_DEFAULT`

       *interactions*
           The default interactions are documented under :data:`INTERACTIONS_DEFAULT`

       *start, stop, step*
           .. seealso:: :class:`mdpow.analysis.ensemble.EnsembleAnalysis`

       :returns:

       *violin plot*
           Returns a :class:`seaborn.FacetGrid` object containing a violin plot of the
           kernel density estimations (KDE) of the dihedral angle frequencies for each
           relevant dihedral atom group in the molecule from the current MDPOW project

       .. rubric:: Example
       
       Typical Workflow::

           import automated_dihedral_analysis as ada

           ada.automated_dihedral_analysis(dirname='/foo/bar/MDPOW_project_data',
                                           figdir='/foo/bar/MDPOW_figure_directory',
                                           resname='UNK', molname='benzene',
                                           padding=45, width=0.9,
                                           solvents=('water','octanol'),
                                           interactions=('Coulomb','VDW'),
                                           start=0, stop=100, step=10)
    '''

    if dataframe is not None:

        df = dataframe
        logger.info(f'Proceeding with results DataFrame provided.')

    else:
        # new arrangement (1/31), temporary
        # this will need to change or the indices will need
        # to be saved in results DF for future changes
        atom_group_indices = dihedral_indices(dirname=dirname, resname=resname, SMARTS=SMARTS)
        df = dihedral_groups_ensemble(atom_group_indices=atom_group_indices, dirname=dirname,
                                      solvents=solvents, interactions=interactions,
                                      start=start, stop=stop, step=step)

    if df_save_dir is not None:
        save_df(df=df, df_save_dir=df_save_dir, resname=resname, molname=molname)

    df_aug = periodic_angle(df, padding=padding)

    # double check all kwargs for consistency/redundancy
    # create dictionary input option for downstream use
    # in workflows.base and workflows.registry
    return plot_violins(df_aug, resname=resname, figdir=figdir, molname=molname, width=width, solvents=solvents)
