import os
import re
import pandas as pd

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