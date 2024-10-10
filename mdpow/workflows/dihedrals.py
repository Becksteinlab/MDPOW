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

Atom indices obtained by :func:`get_atom_indices` are 0-based,
atom index labels on the molecule in the plots are 0-based,
but atom `names` in plots and file names are 1-based.

.. autodata:: SOLVENTS_DEFAULT
    :annotation: = ('water', 'octanol')
.. autodata:: INTERACTIONS_DEFAULT
    :annotation: = ('Coulomb', 'VDW')
.. autodata:: SMARTS_DEFAULT
    :annotation: = [!#1]~[!$(*#*)&!D1]-!@[!$(*#*)&!D1]~[!#1]
.. autodata:: PLOT_WIDTH_DEFAULT
    :annotation: = 190
.. autofunction:: automated_dihedral_analysis
.. autofunction:: build_universe
.. autofunction:: rdkit_conversion
.. autofunction:: get_atom_indices
.. autofunction:: get_bond_indices
.. autofunction:: get_dihedral_groups
.. autofunction:: get_paired_indices
.. autofunction:: dihedral_groups_ensemble
.. autofunction:: save_df
.. autofunction:: periodic_angle_padding
.. autofunction:: dihedral_violins
.. autofunction:: build_svg
.. autofunction:: plot_dihedral_violins

"""

import os
import logging
import pathlib

import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.topology.guessers import guess_types
from rdkit import Chem
from rdkit.Chem import rdCoordGen
from rdkit.Chem.Draw import rdMolDraw2D
import matplotlib.pyplot as plt
import seaborn as sns
import pypdf
import cairosvg
import svgutils
import svgutils.compose
import svgutils.transform

from .base import guess_elements
from ..analysis import ensemble, dihedral

logger = logging.getLogger("mdpow.workflows.dihedrals")

SOLVENTS_DEFAULT = ("water", "octanol")
"""Default solvents are water and octanol:

    * must match solvents used in project directory
    * one or two solvents can be specified
    * current solvents supported,

   .. seealso:: :mod:`mdpow.forcefields`

"""

INTERACTIONS_DEFAULT = ("Coulomb", "VDW")
"""Default interactions set to Coulomb and VDW:

    * default values should not be changed
    * order should not be changed

"""

SMARTS_DEFAULT = "[!#1]~[!$(*#*)&!D1]-!@[!$(*#*)&!D1]~[!#1]"
"""Default SMARTS string to identify relevant dihedral atom groups:

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

PLOT_WIDTH_DEFAULT = 190
"""Plot width (`plot_pdf_width`) should be provided in millimeters (mm),
   and is converted to pixels (px) for use with :mod:`cairosvg`.

   conversion factor: 1 mm = 3.7795275591 px
   default value: 190 mm = 718.110236229 pixels
"""


def build_universe(dirname, solvents=SOLVENTS_DEFAULT):
    """Builds :class:`~MDAnalysis.core.universe.Universe` from the
    ``./Coulomb/0000`` topology and trajectory of the project for
    the first solvent specified.

    Output used by :func:`~mdpow.workflows.dihedrals.rdkit_conversion`
    and :func:`~mdpow.workflows.dihedrals.get_atom_indices` to obtain atom indices
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

     *solvents*
        The default solvents are documented under :data:`SOLVENTS_DEFAULT`.
        Normally takes a two-tuple, but analysis is compatible with single solvent selections.
        Single solvent analyses will result in a figure with fully filled violins for the single solvent.

    :returns:

    *u*
        :class:`~MDAnalysis.core.universe.Universe` object

    """

    path = pathlib.Path(dirname)
    topology = path / f"FEP/{solvents[0]}/Coulomb/0000" / "md.tpr"
    trajectory = path / f"FEP/{solvents[0]}/Coulomb/0000" / "md.xtc"
    u = mda.Universe(str(topology), str(trajectory))

    return u


def rdkit_conversion(u, resname):
    """Converts the solute, `resname`, of the
    :class:`~MDAnalysis.core.universe.Universe` to :class:`rdkit.Chem.rdchem.Mol` object
    for use with a SMARTS selection string to identify dihedral atom groups.

    Accepts :class:`~MDAnalysis.core.universe.Universe` object made with
    :func:`~mdpow.workflows.dihedrals.build_universe` and a `resname` as input.
    Uses `resname` to select the solute for conversion by
    :class:`~MDAnalysis.converters.RDKit.RDKitConverter` to :class:`rdkit.Chem.rdchem.Mol`,
    and will add element attributes for Hydrogen if not listed in the topology,
    using :func:`MDAnalysis.topology.guessers.guess_atom_element`.

    :keywords:

    *u*
        :class:`~MDAnalysis.core.universe.Universe` object

    *resname*
        `resname` for the molecule as defined in
        the topology and trajectory

    :returns:

    *tuple(mol, solute)*
        function call returns tuple, see below

    *mol*
        :class:`rdkit.Chem.rdchem.Mol` object converted from `solute`

    *solute*
        the :any:`MDAnalysis` `AtomGroup` for the solute

    """

    try:
        solute = u.select_atoms(f"resname {resname}")
        mol = solute.convert_to("RDKIT")
    except AttributeError:
        guessed_elements = guess_elements(u.atoms)
        u.add_TopologyAttr("elements", guessed_elements)
        solute = u.select_atoms(f"resname {resname}")
        mol = solute.convert_to("RDKIT")

    rdCoordGen.AddCoords(mol)

    for atom in mol.GetAtoms():
        atom.SetProp("atomNote", str(atom.GetIdx()))

    return mol, solute


def get_atom_indices(mol, SMARTS=SMARTS_DEFAULT):
    """Uses a SMARTS selection string to identify atom indices
    for relevant dihedral atom groups.

    Requires a :class:`rdkit.Chem.rdchem.Mol` object as input
    for the :data:`SMARTS_DEFAULT` kwarg to match patterns to and
    identify relevant dihedral atom groups.

    :keywords:

    *mol*
        :class:`rdkit.Chem.rdchem.Mol` object converted from `solute`

    *SMARTS*
        The default SMARTS string is described in detail under :data:`SMARTS_DEFAULT`.

    :returns:

    *atom_indices*
        tuple of tuples of indices for each dihedral atom group

    """

    pattern = Chem.MolFromSmarts(SMARTS)
    atom_indices = mol.GetSubstructMatches(pattern)

    return atom_indices


def get_bond_indices(mol, atom_indices):
    """From the :class:`rdkit.Chem.rdchem.Mol` object, uses
    `atom_indices` to determine the indices of the bonds
    between those atoms for each dihedral atom group.

    :keywords:

    *mol*
        :class:`rdkit.Chem.rdchem.Mol` object converted from `solute`

    *atom_indices*
        tuple of tuples of indices for each dihedral atom group

    :returns:

    *bond_indices*
        tuple of tuples of indices for the bonds in each dihedral atom group

    """

    bonds = []

    for atom_index in atom_indices:
        x = mol.GetBondBetweenAtoms(atom_index[0], atom_index[1]).GetIdx()
        y = mol.GetBondBetweenAtoms(atom_index[1], atom_index[2]).GetIdx()
        z = mol.GetBondBetweenAtoms(atom_index[2], atom_index[3]).GetIdx()
        bix = (x, y, z)

        bonds.append(bix)

    bond_indices = tuple(bonds)

    return bond_indices


def get_dihedral_groups(solute, atom_indices):
    """Uses the 0-based `atom_indices` of the relevant dihedral atom groups
    determined by :func:`~mdpow.workflows.dihedrals.get_atom_indices`
    and returns the 1-based index names for each atom in each group.

    Requires the `atom_indices` from :func:`~mdpow.workflows.dihedrals.get_atom_indices`
    to index the `solute` specified by :func:`~MDAnalysis.core.groups.select_atoms`
    and return an array of the names of each atom within its respective
    dihedral atom group as identified by the SMARTS selection string.

    :keywords:

    *solute*
        the :any:`MDAnalysis` `AtomGroup` for the solute

    *atom_indices*
        tuple of tuples of indices for each dihedral atom group

    :returns:

    *dihedral_groups*
        list of :func:`numpy.array` for atom names in each dihedral atom group

    """
    # currently uses RDKit Mol object atom indices to retrieve
    # atom names from the MDAnalysis solute object
    # RDKit-MDAnalysis index consistency is currently tested
    dihedral_groups = [
        solute.atoms[list(atom_index)].names for atom_index in atom_indices
    ]

    return dihedral_groups


def get_paired_indices(atom_indices, bond_indices, dihedral_groups):
    """Combines `atom_indices` and `bond_indices` in tuples
    to be paired with their respective dihedral atom groups.

    A dictionary is created with key-value pairs as follows:
    `atom_indices` and `bond_indices` are joined in a tuple
    as the value, with the key being the respective member
    of `dihedral_groups` to facilitate highlighting the
    relevant dihedral atom group when generating violin plots.
    As an example, `'C1-N2-O3-S4': ((0, 1, 2, 3), (0, 1, 2))`,
    would be one key-value pair in the dictionary.

    :keywords:

    *atom_indices*
        tuple of tuples of indices for each dihedral atom group

    *bond_indices*
        tuple of tuples of indices for the bonds in each dihedral atom group

    *dihedral_groups*
        list of :func:`numpy.array` for atom names in each dihedral atom group

    :returns:

    *name_index_pairs*
        dictionary with key-value pair for dihedral atom group,
        atom indices, and bond indices

    """

    all_dgs = [f"{dg[0]}-{dg[1]}-{dg[2]}-{dg[3]}" for dg in dihedral_groups]

    name_index_pairs = {}
    name_index_pairs = {
        all_dgs[i]: (atom_indices[i], bond_indices[i]) for i in range(len(all_dgs))
    }

    return name_index_pairs


def dihedral_groups_ensemble(
    dirname,
    atom_indices,
    solvents=SOLVENTS_DEFAULT,
    interactions=INTERACTIONS_DEFAULT,
    start=None,
    stop=None,
    step=None,
):
    """Creates one :class:`~mdpow.analysis.ensemble.Ensemble` for the MDPOW
    project and runs :class:`~mdpow.analysis.dihedral.DihedralAnalysis`
    for each dihedral atom group identified by the SMARTS
    selection string.

    .. seealso::

       :func:`~mdpow.workflows.dihedrals.automated_dihedral_analysis`,
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

    *atom_indices*
        tuples of atom indices for dihedral atom groups

        .. seealso:: :func:`~mdpow.workflows.dihedrals.get_atom_indices`, :data:`SMARTS_DEFAULT`

    *solvents*
        The default solvents are documented under :data:`SOLVENTS_DEFAULT`.
        Normally takes a two-tuple, but analysis is compatible with single solvent selections.
        Single solvent analyses will result in a figure with fully filled violins for the single solvent.

    *interactions*
        The default interactions are documented under :data:`INTERACTIONS_DEFAULT`.

    *start, stop, step*
        arguments passed to :func:`~mdpow.analysis.ensemble.EnsembleAnalysis.run`,
        as parameters for iterating through the trajectories of the current ensemble

        .. seealso:: :class:`~mdpow.analysis.ensemble.EnsembleAnalysis`

    :returns:

    *df*
        :class:`pandas.DataFrame` of :class:`~mdpow.analysis.dihedral.DihedralAnalysis`
        results, including all dihedral atom groups for molecule of current project

    """

    dih_ens = ensemble.Ensemble(
        dirname=dirname, solvents=solvents, interactions=interactions
    )
    indices = atom_indices
    all_dihedrals = [
        dih_ens.select_atoms(
            f"index {i[0]}", f"index {i[1]}", f"index {i[2]}", f"index {i[3]}"
        )
        for i in indices
    ]

    da = dihedral.DihedralAnalysis(all_dihedrals)
    da.run(start=start, stop=stop, step=step)
    df = da.results

    return df


def save_df(df, df_save_dir, resname, molname=None):
    """Takes a :class:`pandas.DataFrame` of results from
    :class:`~mdpow.analysis.dihedral.DihedralAnalysis`
    as input before padding the angles to optionaly save the raw
    data.

    Optionally saves results before padding the angles for periodicity
    and plotting dihedral angle frequencies as KDE violins
    with :func:`~mdpow.workflows.dihedrals.dihedral_violins`.
    Given a parent directory, creates subdirectory for molecule,
    saves fully sampled, unpadded results :class:`pandas.DataFrame`
    as a compressed csv file, default: .csv.bz2.

    :keywords:

    *df*
        :class:`pandas.DataFrame` of :class:`~mdpow.analysis.dihedral.DihedralAnalysis`
        results, including all dihedral atom groups for molecule of current project

    *df_save_dir*
        optional, path to the location to save results :class:`pandas.DataFrame`

    *resname*
        `resname` for the molecule as defined in
        the topology and trajectory

    *molname*
        molecule name to be used for labelling
        plots, if different from `resname`

    """

    df = df.sort_values(
        by=["selection", "solvent", "interaction", "lambda"]
    ).reset_index(drop=True)

    if molname is None:
        molname = resname

    subdir = molname
    newdir = os.path.join(df_save_dir, subdir)
    os.mkdir(newdir)

    # time and compression level can be adjusted as kwargs
    df.to_csv(f"{newdir}/{molname}_full_df.csv.bz2", index=False, compression="bz2")

    logger.info(f"Results DataFrame saved as {newdir}/{molname}_full_df.csv.bz2")


def periodic_angle_padding(df, padding=45):
    """Pads the angles from the results :class:`~pandas.DataFrame`
    to maintain periodicity in the violin plots.

    Takes a :class:`pandas.DataFrame` of results from
    :class:`~mdpow.analysis.dihedral.DihedralAnalysis`
    or :func:`~mdpow.workflows.dihedrals.dihedral_groups_ensemble`
    as input and pads the angles to maintain periodicity
    for properly plotting dihedral angle frequencies as KDE violins
    with :func:`~mdpow.workflows.dihedrals.dihedral_violins` and
    :func:`~mdpow.workflows.dihedrals.plot_dihedral_violins`.
    Creates two new :class:`pandas.DataFrame` based on the
    `padding` value specified, pads the angle values, concatenates
    all three :class:`pandas.DataFrame`, maintaining original data and
    adding padded values, and returns new augmented :class:`pandas.DataFrame`.

    :keywords:

    *df*
        :class:`pandas.DataFrame` of :class:`~mdpow.analysis.dihedral.DihedralAnalysis`
        results, including all dihedral atom groups for molecule of current project

    *padding*
        value in degrees to specify angle augmentation threshold
        default: 45

    :returns:

    *df_aug*
        augmented results :class:`pandas.DataFrame` containing
        padded dihedral angles as specified by `padding`

    """

    df1 = df[df.dihedral > 180 - padding].copy(deep=True)
    df1.dihedral -= 360
    df2 = df[df.dihedral < -180 + padding].copy(deep=True)
    df2.dihedral += 360
    df_aug = pd.concat([df1, df, df2]).reset_index(drop=True)

    return df_aug


def dihedral_violins(df, width=0.9, solvents=SOLVENTS_DEFAULT, plot_title=None):
    """Plots kernel density estimates (KDE) of dihedral angle frequencies for
    one dihedral atom group as violin plots, using as input the augmented
    :class:`pandas.DataFrame` from
    :func:`~mdpow.workflows.dihedrals.periodic_angle_padding`.

    Output is converted to SVG by :func:`~mdpow.workflows.dihedrals.build_svg`
    and final output is saved as PDF by
    :func:`~mdpow.workflows.dihedrals.plot_dihedral_violins`

    :keywords:

    *df*
        augmented results :class:`pandas.DataFrame` from
        :func:`~mdpow.workflows.dihedrals.periodic_angle_padding`

    *width*
        width of the violin element (>1 overlaps); default: 0.9

    *solvents*
        The default solvents are documented under :data:`SOLVENTS_DEFAULT`.
        Normally takes a two-tuple, but analysis is compatible with single
        solvent selections.  Single solvent analyses will result in a figure
        with fully filled violins for the single solvent.

    *plot_title*
        generated by :func:`~mdpow.workflows.dihedrals.build_svg` using
        `molname`, `dihedral_groups`, `atom_indices`, and `interactions`
        in this order and format: f'{molname}, {name[0]} {a} | ''{col_name}'

    :returns:

    *g*
        returns a :class:`seaborn.FacetGrid` object containing a violin plot of the
        kernel density estimates (KDE) of the dihedral angle frequencies for each
        dihedral atom group identified by :data:`SMARTS_DEFAULT`

    """

    df["lambda"] = df["lambda"].astype("float") / 1000
    df = df.sort_values(
        by=["selection", "solvent", "interaction", "lambda"]
    ).reset_index(drop=True)

    width_ratios = [
        len(df[df["interaction"] == "Coulomb"]["lambda"].unique()) + 1,
        len(df[df["interaction"] == "VDW"]["lambda"].unique()),
        len(df[df["interaction"] == "Coulomb"]["lambda"].unique()) - 1,
    ]

    # Usage in Jupyter causes matplotlib figure object output, not the completed figure
    # Upcoming fix in issue #260
    assert (
        0 < len(solvents) < 3
    ), "one or two solvents must be specified, otherwise SOLVENTS_DEFAULT is used"
    split = len(solvents) > 1
    g = sns.catplot(
        data=df,
        x="lambda",
        y="dihedral",
        hue="solvent",
        col="interaction",
        kind="violin",
        split=split,
        width=width,
        inner=None,
        cut=0,
        linewidth=0.5,
        hue_order=list(solvents),
        col_order=["Coulomb", "VDW", "Structure"],
        sharex=False,
        sharey=True,
        height=3.0,
        aspect=2.0,
        facet_kws={
            "ylim": (-180, 180),
            "gridspec_kws": {
                "width_ratios": width_ratios,
            },
        },
    )

    g.set_xlabels(r"$\lambda$")
    g.set_ylabels(r"dihedral angle $\phi$")
    g.despine(offset=5)

    axC = g.axes_dict["Coulomb"]
    axC.yaxis.set_major_locator(plt.matplotlib.ticker.MultipleLocator(60))
    axC.yaxis.set_minor_locator(plt.matplotlib.ticker.MultipleLocator(30))
    axC.yaxis.set_major_formatter(
        plt.matplotlib.ticker.FormatStrFormatter(r"$%g^\circ$")
    )

    axV = g.axes_dict["VDW"]
    axV.yaxis.set_visible(False)
    axV.spines["left"].set_visible(False)

    axIM = g.axes_dict["Structure"]
    axIM.axis("off")

    g.set_titles(plot_title)

    return g


def build_svg(
    mol,
    molname,
    name_index_pairs,
    atom_group_selection,
    solvents=SOLVENTS_DEFAULT,
    width=0.9,
):
    """Converts and combines figure components into an SVG object to be
    converted and saved as a publication quality PDF.

    :keywords:

    *mol*
        :class:`rdkit.Chem.rdchem.Mol` object converted from `solute`

    *molname*
        molecule name to be used for labelling plots, if different from
        `resname` (in this case, carried over from an upstream decision between
        the two)

    *name_index_pairs*
        dictionary with key-value pair for dihedral atom group, atom indices,
        and bond indices

        .. seealso:: :func:`~mdpow.workflows.dihedrals.get_paired_indices`

    *atom_group_selection*
        `name` of each section in the `groupby` series of atom group selections

        .. seealso:: :func:`~mdpow.workflows.dihedrals.plot_dihedral_violins`

    *solvents*
        The default solvents are documented under :data:`SOLVENTS_DEFAULT`.
        Normally takes a two-tuple, but analysis is compatible with single
        solvent selections.  Single solvent analyses will result in a figure
        with fully filled violins for the single solvent.

    *width*
        width of the violin element (>1 overlaps); default: 0.9

    :returns:

    *fig*
        :mod:`svgutils` SVG figure object

    """

    atom_index = name_index_pairs[atom_group_selection[0]][0]
    bond_index = name_index_pairs[atom_group_selection[0]][1]

    drawer = rdMolDraw2D.MolDraw2DSVG(250, 250)
    drawer.DrawMolecule(mol=mol, highlightAtoms=atom_index, highlightBonds=bond_index)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText().replace("svg:", "")

    mol_svg = svgutils.transform.fromstring(svg)
    m = mol_svg.getroot()
    m.scale(0.0125).moveto(21.8, 0.35)

    plot_title = f"{molname}, {atom_group_selection[0]} {atom_index} | " "{col_name}"
    plot = dihedral_violins(
        atom_group_selection[1], width=width, solvents=solvents, plot_title=plot_title
    )

    plot_svg = svgutils.transform.from_mpl(plot)
    p = plot_svg.getroot()
    p.scale(0.02)

    # order matters here, plot down first, mol on top (p, m)
    fig = svgutils.compose.Figure("28cm", "4.2cm", p, m)

    return fig


def plot_dihedral_violins(
    df,
    resname,
    mol,
    name_index_pairs,
    figdir=None,
    molname=None,
    width=0.9,
    plot_pdf_width=PLOT_WIDTH_DEFAULT,
    solvents=SOLVENTS_DEFAULT,
):
    """Coordinates plotting and saving figures for all dihedral atom groups.

    Makes a subdirectory for the current project within the specified
    `figdir` using `resname` or `molname` as title and saves production
    quality PDFs for each dihedral atom group separately.

    .. seealso::

       :func:`~mdpow.workflows.dihedrals.automated_dihedral_analysis`,
       :func:`~mdpow.workflows.dihedrals.dihedral_violins`,
       :func:`~mdpow.workflows.dihedrals.build_svg`

    :keywords:

    *df*
        augmented results :class:`pandas.DataFrame` from
        :func:`~mdpow.workflows.dihedrals.periodic_angle_padding`

    *resname*
        `resname` for the molecule as defined in
        the topology and trajectory

    *mol*
        :class:`rdkit.Chem.rdchem.Mol` object converted from `solute`

    *name_index_pairs*
        dictionary with key-value pair for dihedral atom group,
        atom indices, and bond indices

        .. seealso:: :func:`~mdpow.workflows.dihedrals.get_paired_indices`

    *figdir*
        path to the location to save figures (REQUIRED but marked
        as a kwarg for technical reasons; will be changed in #244)

    *molname*
        molecule name to be used for labelling
        plots, if different from `resname`

    *width*
        width of the violin element (>1 overlaps)
        default: 0.9

        .. seealso:: :func:`~mdpow.workflows.dihedrals.dihedral_violins`

    *plot_pdf_width*
        The default value for width of plot output is described in detail under
        :data:`PLOT_WIDTH_DEFAULT`.

    *solvents*
        The default solvents are documented under :data:`SOLVENTS_DEFAULT`.
        Normally takes a two-tuple, but analysis is compatible with single solvent selections.
        Single solvent analyses will result in a figure with fully filled violins for the single solvent.

    """

    assert (
        figdir is not None
    ), "figdir MUST be set, even though it is a kwarg. Will be changed with #244"

    if molname is None:
        molname = resname

    subdir = molname
    newdir = os.path.join(figdir, subdir)
    os.mkdir(newdir)

    section = df.groupby(by="selection")

    plot_pdf_width_px = plot_pdf_width * 3.7795275591

    pdf_list = []

    for name in section:
        fig = build_svg(
            mol=mol,
            molname=molname,
            atom_group_selection=name,
            name_index_pairs=name_index_pairs,
            solvents=solvents,
            width=width,
        )

        figfile = pathlib.Path(newdir) / f"{molname}_{name[0]}_violins.pdf"
        if figdir is not None:
            plot_pdf = cairosvg.svg2pdf(
                bytestring=fig.tostr(),
                write_to=str(figfile),
                output_width=plot_pdf_width_px,
            )

            # add PDF for each dihedral atom group to all_PDFs list
        pdf_list.append(f"{figfile}")

        logger.info(f"Figure saved as {figfile}")

    logger.info(f"All figures generated and saved in {figdir}")

    merger = pypdf.PdfWriter()

    for pdf in pdf_list:
        merger.append(pdf)
    merger.write(f"{figdir}/{molname}/{molname}_all_figures.pdf")
    merger.close()
    logger.info(
        f"PDF of combined figures generated and saved as {figdir}/{molname}/{molname}_all_figures.pdf"
    )

    return None


def automated_dihedral_analysis(
    dirname,
    resname,
    figdir=None,
    # figdir is required and will cause issues if not specified
    # figdir=None is a temporary way to satisfy
    # workflows base tests until issue #244 is resolved
    # because it currently uses a **kwargs convention and the
    # positional argument figdir will not carry over nicely
    df_save_dir=None,
    molname=None,
    SMARTS=SMARTS_DEFAULT,
    plot_pdf_width=PLOT_WIDTH_DEFAULT,
    dataframe=None,
    padding=45,
    width=0.9,
    solvents=SOLVENTS_DEFAULT,
    interactions=INTERACTIONS_DEFAULT,
    start=None,
    stop=None,
    step=None,
):
    """Runs :class:`~mdpow.analysis.dihedral.DihedralAnalysis` for a single MDPOW
    project and creates violin plots of dihedral angle frequencies for each
    relevant dihedral atom group.

    For one MDPOW project, automatically determines all relevant dihedral atom groups
    in the molecule, runs :class:`~mdpow.analysis.dihedral.DihedralAnalysis` for each group,
    pads the dihedral angles to maintain periodicity, creates violin plots of dihedral angle
    frequencies (KDEs), and saves publication quality PDF figures for each group, separately.

    Optionally saves all pre-padded :class:`~mdpow.analysis.dihedral.DihedralAnalysis` results
    as a single :class:`pandas.DataFrame` in `df_save_dir` provided.

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

    *figdir*
        path to the location to save figures (REQUIRED but marked
        as a kwarg for technical reasons; will be changed in #244)

    *resname*
        `resname` for the molecule as defined in
        the topology and trajectory

    *df_save_dir*
        optional, path to the location to save results :class:`pandas.DataFrame`

    *molname*
        molecule name to be used for labelling
        plots, if different from `resname`

    *SMARTS*
        The default SMARTS string is described in detail under :data:`SMARTS_DEFAULT`.

    *plot_pdf_width*
        The default value for width of plot output is described in detail under
        :data:`PLOT_WIDTH_DEFAULT`.

    *dataframe*
        optional, if :class:`~mdpow.analysis.dihedral.DihedralAnalysis`
        was done prior, then results :class:`pandas.DataFrame` can be
        input to utilize angle padding and violin plotting functionality

    *padding*
        value in degrees
        default: 45

        .. seealso:: :func:`~mdpow.workflows.dihedrals.periodic_angle_padding`

    *width*
        width of the violin element (>1 overlaps)
        default: 0.9

        .. seealso:: :func:`~mdpow.workflows.dihedrals.dihedral_violins`

    *solvents*
        The default solvents are documented under :data:`SOLVENTS_DEFAULT`.
        Normally takes a two-tuple, but analysis is compatible with single solvent selections.
        Single solvent analyses will result in a figure with fully filled violins for the single solvent.

    *interactions*
        The default interactions are documented under :data:`INTERACTIONS_DEFAULT`.

    *start, stop, step*
        arguments passed to :func:`~mdpow.analysis.ensemble.EnsembleAnalysis.run`,
        as parameters for iterating through the trajectories of the current ensemble

        .. seealso:: :class:`~mdpow.analysis.ensemble.EnsembleAnalysis`

    .. rubric:: Example

    Typical Workflow::

        import dihedrals

        dihedrals.automated_dihedral_analysis(dirname='/foo/bar/MDPOW_project_data',
                                        figdir='/foo/bar/MDPOW_figure_directory',
                                        resname='UNK', molname='benzene',
                                        padding=45, width=0.9,
                                        solvents=('water','octanol'),
                                        interactions=('Coulomb','VDW'),
                                        start=0, stop=100, step=10)
    """

    u = build_universe(dirname=dirname, solvents=solvents)
    mol, solute = rdkit_conversion(u=u, resname=resname)
    atom_indices = get_atom_indices(mol=mol, SMARTS=SMARTS)
    bond_indices = get_bond_indices(mol=mol, atom_indices=atom_indices)
    dihedral_groups = get_dihedral_groups(solute=solute, atom_indices=atom_indices)
    name_index_pairs = get_paired_indices(
        atom_indices=atom_indices,
        bond_indices=bond_indices,
        dihedral_groups=dihedral_groups,
    )

    if dataframe is not None:
        df = dataframe
        logger.info(f"Proceeding with results DataFrame provided.")

    else:
        df = dihedral_groups_ensemble(
            dirname=dirname,
            atom_indices=atom_indices,
            solvents=solvents,
            interactions=interactions,
            start=start,
            stop=stop,
            step=step,
        )

    if df_save_dir is not None:
        save_df(df=df, df_save_dir=df_save_dir, resname=resname, molname=molname)

    df_aug = periodic_angle_padding(df, padding=padding)

    plot_dihedral_violins(
        df=df_aug,
        resname=resname,
        mol=mol,
        name_index_pairs=name_index_pairs,
        figdir=figdir,
        molname=molname,
        width=width,
        plot_pdf_width=plot_pdf_width,
        solvents=solvents,
    )

    logger.info(f"DihedralAnalysis completed for all projects in {dirname}")

    return
