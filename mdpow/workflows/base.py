# MDPOW: base.py
# 2022 Cade Duckworth

"""
:mod:`mdpow.workflows.base` --- Automated workflow base functions
=================================================================

To analyze multiple MDPOW projects, provide :func:`project_paths`
with the top-level directory containing all MDPOW projects' simulation data
to obtain a :class:`pandas.DataFrame` containing the project information
and paths. Then, :func:`automated_project_analysis` takes as input the
aforementioned :class:`pandas.DataFrame` and runs the specified
:class:`~mdpow.analysis.ensemble.EnsembleAnalysis` for all MDPOW projects
under the top-level directory provided to :func:`project_paths`.

.. seealso:: :mod:`~mdpow.workflows.registry`

.. autofunction:: project_paths
.. autofunction:: automated_project_analysis
.. autofunction:: guess_elements

"""

import os
import re
import logging

import numpy as np
import pandas as pd
from MDAnalysis.topology import guessers, tables

logger = logging.getLogger("mdpow.workflows.base")


def project_paths(parent_directory=None, csv=None, csv_save_dir=None):
    """Takes a top directory containing MDPOW projects and determines
    the molname, resname, and path, of each MDPOW project within.

    Optionally takes a .csv file containing `molname`, `resname`, and
    `paths`, in that order.

    :keywords:

    *parent_directory*
        the path for the location of the top directory
        under which the subdirectories of MDPOW simulation
        data exist, additionally creates a 'project_paths.csv' file
        for user manipulation of metadata and for future reference

    *csv*
        .csv file containing the molecule names, resnames,
        and paths, in that order, for the MDPOW simulation
        data to be iterated over must contain header of the
        form: `molecule,resname,path`

    *csv_save_dir*
        optionally provided directory to save .csv file, otherwise,
        data will be saved in current working directory

    :returns:

    *project_paths*
        :class:`pandas.DataFrame` containing MDPOW project metadata

    .. rubric:: Example

    Typical Workflow::

        project_paths = project_paths(parent_directory='/foo/bar/MDPOW_projects')
        automated_project_analysis(project_paths)

    or::

        project_paths = project_paths(csv='/foo/bar/MDPOW.csv')
        automated_project_analysis(project_paths)

    """

    if parent_directory is not None:
        locations = []

        reg_compile = re.compile("FEP")
        for dirpath, dirnames, filenames in os.walk(parent_directory):
            result = [
                dirpath.strip() for dirname in dirnames if reg_compile.match(dirname)
            ]
            if result:
                locations.append(result[0])

        resnames = []

        for loc in locations:
            res_temp = loc.strip().split("/")
            resnames.append(res_temp[-1])

        project_paths = (
            pd.DataFrame({"molecule": resnames, "resname": resnames, "path": locations})
            .sort_values(by=["molecule", "resname", "path"])
            .reset_index(drop=True)
        )
        if csv_save_dir is not None:
            project_paths.to_csv(f"{csv_save_dir}/project_paths.csv", index=False)
            logger.info(f"project_paths saved under {csv_save_dir}")
        else:
            current_directory = os.getcwd()
            project_paths.to_csv("project_paths.csv", index=False)
            logger.info(f"project_paths saved under {current_directory}")

    elif csv is not None:
        locations = pd.read_csv(csv)
        project_paths = locations.sort_values(
            by=["molecule", "resname", "path"]
        ).reset_index(drop=True)

    return project_paths


def automated_project_analysis(project_paths, ensemble_analysis, **kwargs):
    """Takes a :class:`pandas.DataFrame` created by :func:`~mdpow.workflows.base.project_paths`
    and iteratively runs the specified :class:`~mdpow.analysis.ensemble.EnsembleAnalysis`
    for each of the projects by running the associated automated workflow
    in each project directory returned by :func:`~mdpow.workflows.base.project_paths`.

    Compatibility with more automated analyses in development.

    :keywords:

    *project_paths*
        :class:`pandas.DataFrame` that provides paths to MDPOW projects

    *ensemble_analysis*
        name of the :class:`~mdpow.analysis.ensemble.EnsembleAnalysis`
        that corresponds to the desired automated workflow module

    *kwargs*
        keyword arguments for the supported automated workflows,
        see the :mod:`~mdpow.workflows.registry` for all available
        workflows and their call signatures

    .. rubric:: Example

    A typical workflow is the automated dihedral analysis from
    :mod:`mdpow.workflows.dihedrals`, which applies the *ensemble analysis*
    :class:`~mdpow.analysis.dihedral.DihedralAnalysis` to each project.
    The :data:`~mdpow.workflows.registry.registry` contains this automated
    workflow under the key *"DihedralAnalysis"* and so the automated execution
    for all `project_paths` (obtained via :func:`project_paths`) is performed by
    passing the specific key to :func:`automated_project_analysis`::

        project_paths = project_paths(parent_directory='/foo/bar/MDPOW_projects')
        automated_project_analysis(project_paths, ensemble_analysis='DihedralAnalysis', **kwargs)

    """
    # import inside function to avoid circular imports
    from .registry import registry

    for row in project_paths.itertuples():
        molname = row.molecule
        resname = row.resname
        dirname = row.path

        logger.info(f"starting {molname}")

        try:
            registry[ensemble_analysis](
                dirname=dirname, resname=resname, molname=molname, **kwargs
            )
            logger.info(f"{molname} completed")
        except KeyError as err:
            msg = (
                f"Invalid ensemble_analysis {err}. An EnsembleAnalysis type that corresponds "
                "to an existing automated workflow module must be input as a kwarg. "
                "ex: ensemble_analysis='DihedralAnalysis'"
            )
            logger.error(f"{err} is an invalid selection")
            raise KeyError(msg)
        except TypeError as err:
            msg = (
                f"Invalid ensemble_analysis {ensemble_analysis}. An EnsembleAnalysis type that "
                "corresponds to an existing automated workflow module must be input as a kwarg. "
                "ex: ensemble_analysis='DihedralAnalysis'"
            )
            logger.error(f"workflow module for {ensemble_analysis} does not exist yet")
            raise TypeError(msg)

    logger.info("all analyses completed")
    return


def guess_elements(atoms, rtol=1e-3):
    """guess elements for atoms from masses

    Given masses, we perform a reverse lookup on
    :data:`MDAnalysis.topology.tables.masses` to find the corresponding
    element. Only atoms where the standard MDAnalysis guesser finds elements
    with masses contradicting the topology masses are corrected.

    .. Note:: This function *requires* correct masses to be present.
              No sanity checks because MDPOW always uses TPR files that
              contain correct masses.

    :arguments:

    *atoms*
         MDAnalysis AtomGroup *with masses defined*

    :keywords:

    *rtol*
         relative tolerance for a match (as used in :func:`numpy.isclose`);
         atol=1e-6 is at a fixed value, which means that "zero" is only
         recognized for values =< 1e-6

         .. note:: In order to reliably match GROMACS masses, *rtol* should
                   be at least 1e-3.

    :returns:

    *elements*
         array of guessed element symbols, in same order as `atoms`

    .. rubric:: Example

    As an example we guess masses and then set the elements for all atoms::

       elements = guess_elements(atoms)
       atoms.add_TopologyAttr("elements", elements)

    """
    ATOL = 1e-6

    names = atoms.names
    masses = atoms.masses

    mda_elements = np.fromiter(tables.masses.keys(), dtype="U5")
    mda_masses = np.fromiter(tables.masses.values(), dtype=np.float64)

    guessed_elements = guessers.guess_types(names)
    guessed_masses = np.array([guessers.get_atom_mass(n) for n in guessed_elements])
    problems = np.logical_not(np.isclose(masses, guessed_masses, atol=ATOL, rtol=rtol))

    # match only problematic  masses against the MDA reference masses
    iproblem, ielem = np.nonzero(
        np.isclose(masses[problems, np.newaxis], mda_masses, atol=ATOL, rtol=rtol)
    )
    # We should normally find a match for each problem but just in case, assert and
    # give some useful information for debugging.
    assert len(ielem) == sum(problems), (
        "Not all masses could be assigned an element, "
        f"missing names {set(names[problems]) - set(names[problems][iproblem])}"
    )

    guessed_elements[problems] = mda_elements[ielem]

    # manually fix some dummies that are labelled "D": set ALL zero masses to DUMMY
    guessed_elements[np.isclose(masses, 0, atol=ATOL)] = "DUMMY"

    return guessed_elements
