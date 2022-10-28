import automated_dihedral_analysis as ada



def directory_iteration(dirpaths, figdir=None, padding=45, width=0.9,
                        solvents=('water','octanol'), interactions=('Coulomb','VDW'),
                        start=None, stop=None, step=None):
    #fix to update how I end up handling arguments

    for row in dirpaths.itertuples():
            molname = row.molecule
            resname = row.resname
            datadir = row.path

            ada.automated_dihedral_analysis(datadir=datadir, figdir=figdir,
                                            molname=molname, resname=resname,
                                            padding=padding, width=width,
                                            solvents=solvents, interactions=interactions,
                                            start=start, stop=stop, step=step)

    return
