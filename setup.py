# setuptools installation of POW
# Copyright (c) 2010 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)


from setuptools import setup, find_packages

# Dynamically calculate the version based on VERSION.
version = __import__('mdpow.version').get_version()

setup(name="POW",
      version=version,
      description="A library for computing octanol/water partitioning coefficients",
      long_description="""The POW module simplifies the setup and
execution of free energy calculations of small molecules in water and
and other solvents such as octanol and cyclohexane.
It uses Gromacs (http://www.gromacs.org) for the molecular dynamics
(MD) simulations and relies on GromacsWrapper
(https://github.com/Becksteinlab/GromacsWrapper)
""",
      author="Oliver Beckstein",
      author_email="orbeckst@gmail.com",
      license="GPLv3",
      url="https://github.com/Becksteinlab/MDPOW",
      keywords="science Gromacs analysis 'molecular dynamics'",
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
                              'templates/*'], },
      install_requires=['numpy>=1.6', 'scipy',
                        'pyyaml',
                        'GromacsWrapper>=0.5.1',
                        'numkit',
                        'six',
      ],
      setup_requires=['pytest-runner',],
      # alchemtest does not have a release yet
      #dependency_links = ["https://github.com/alchemistry/alchemtest/tarball/master#egg=alchemtest"],
      tests_require=['pytest', 'pybol', 'py'],
      zip_safe=True,
)


