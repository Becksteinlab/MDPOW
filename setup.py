# setuptools installation of POW
# Copyright (c) 2010 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)


from setuptools import setup, find_packages

# Dynamically calculate the version based on VERSION.
version = __import__('mdpow.version').get_version()

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
      packages=find_packages(exclude=['examples']),
      scripts=['scripts/mdpow-pow',
               'scripts/mdpow-pcw',
               'scripts/mdpow-ghyd',
               'scripts/mdpow-check',
               'scripts/mdpow-rebuild-fep',
               'scripts/mdpow-rebuild-simulation',
               'scripts/mdpow-equilibrium',
               'scripts/mdpow-fep',
               'scripts/mdpow-cfg2yaml.py',
               'scripts/mdpow-solvationenergy',
               'scripts/mdpow-get-runinput'
      ],
      package_data={'mdpow': ['top/*.dat', 'top/*.gro', 'top/*.itp',
                              'top/oplsaa.ff/*',
                              'top/charmm36-mar2019.ff/*',
                              'top/amber99sb.ff/*',
                              'templates/*'], },
      install_requires=['numpy>=1.6', 'scipy',
                        'pyyaml',
                        'GromacsWrapper>=0.5.1',
                        'numkit',
                        'six',
                        'mdanalysis',
                        'alchemlyb',
                        'pandas',
                        'pymbar',
      ],
      #setup_requires=['pytest-runner',],
      tests_require=['pytest', 'pybol', 'py'],
      zip_safe=True,
)
