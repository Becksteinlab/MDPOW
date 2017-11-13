# setuptools installation of POW
# Copyright (c) 2010 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)


from setuptools import setup, find_packages

# Dynamically calculate the version based on VERSION.
from importlib import import_module
import sys
sys.path.insert(0, "./src")
try:
    version = import_module('mdpow.version').get_version()
finally:
    sys.path.pop()

setup(name="MDPOW",
      version=version,
      description="A library for computing solvation/water partitioning coefficients using molecular dynamics simulations",
      long_description=open("README.rst").read(),
      author="Oliver Beckstein",
      author_email="orbeckst@gmail.com",
      license="GPLv3",
      url="https://github.com/Becksteinlab/MDPOW",
      keywords="science Gromacs analysis 'molecular dynamics'",
      classifiers=[
          "Development Status :: 4 - Beta",
          "Environment :: Console",
          "Intended Audience :: Science/Research",
          "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
          "Operating System :: POSIX",
          "Programming Language :: Python :: 2.7",
          "Topic :: Scientific/Engineering :: Chemistry",
          "Topic :: Scientific/Engineering :: Physics",
      ],
      packages=find_packages('src', exclude=['examples']),
      package_dir={'': 'src'},
      scripts=['src/scripts/mdpow-pow',
               'src/scripts/mdpow-pcw',
               'src/scripts/mdpow-ghyd',
               'src/scripts/mdpow-check',
               'src/scripts/mdpow-rebuild-fep',
               'src/scripts/mdpow-rebuild-simulation',
               'src/scripts/mdpow-equilibrium',
               'src/scripts/mdpow-fep',
               'src/scripts/mdpow-cfg2yaml.py',
               'src/scripts/mdpow-solvationenergy',
               'src/scripts/mdpow-get-runinput'
      ],
      package_data={'mdpow': ['top/*.dat', 'top/*.gro', 'top/*.itp',
                              'top/oplsaa.ff/*',
                              'templates/*'], },
      install_requires=['numpy>=1.6', 'scipy',
                        'pyyaml',
                        'GromacsWrapper>=0.5.1',
                        'numkit',
                        'six',
      ],
      setup_requires=['pytest-runner',],
      tests_require=['pytest', 'pybol', 'py'],
      zip_safe=True,
)


