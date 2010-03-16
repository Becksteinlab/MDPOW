# setuptools installation of POW
# Copyright (c) 2010 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)


from ez_setup import use_setuptools
use_setuptools()
from setuptools import setup, find_packages

setup(name="POW",
      version="0.0.2",
      description="A library for computing octanol/water partitioning coefficients",
      long_description="""The POW module simplifies the setup and
execution of free energy calculations of small molecules in water and
octanol. It uses Gromacs (http://www.gromacs.org) for the molecular dynamics
(MD) simulations and relies on GromacsWrapper
(http://sbcb.bioch.ox.ac.uk/oliver/software/GromacsWrapper).
""",
      author="Oliver Beckstein",
      author_email="orbeckst@gmail.com",
      license="GPLv3",
      url="http://sbcb.bioch.ox.ac.uk/oliver/software/POW/",
      download_url="http://sbcb.bioch.ox.ac.uk/oliver/download/Python/",
      keywords="science Gromacs analysis 'molecular dynamics'",
      packages=find_packages(exclude=['examples']),
      package_data={'mdpow': ['top/*', 'templates/*'], },
      install_requires = ['numpy>=1.0',
                          'GromacsWrapper>=0.1.1'],
      dependency_links = ["http://sbcb.bioch.ox.ac.uk/oliver/download/Python/"],
      zip_safe = True,
)

      
