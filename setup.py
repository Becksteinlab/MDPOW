# setuptools installation of POW
# Copyright (c) 2010 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)

from setuptools import setup, find_packages
import versioneer

setup(
    name="MDPOW",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="A library for computing solvation/water partitioning coefficients using molecular dynamics simulations",
    long_description=open("README.rst").read(),
    long_description_content_type="text/x-rst",
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
        "Operating System :: MacOS :: MacOS X",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    packages=find_packages(exclude=["examples"]),
    scripts=[
        "scripts/mdpow-pow",
        "scripts/mdpow-pcw",
        "scripts/mdpow-ptw",
        "scripts/mdpow-check",
        "scripts/mdpow-rebuild-fep",
        "scripts/mdpow-rebuild-simulation",
        "scripts/mdpow-equilibrium",
        "scripts/mdpow-fep",
        "scripts/mdpow-cfg2yaml.py",
        "scripts/mdpow-solvationenergy",
        "scripts/mdpow-get-runinput",
    ],
    # exclude large data sets in tests/testing_resources/*
    include_package_data=True,
    package_data={
        "mdpow": [
            "top/*.dat",
            "top/*.gro",
            "top/*.itp",
            "top/oplsaa.ff/*",
            "top/charmm36-mar2019.ff/*",
            "top/amber99sb.ff/*",
            "templates/*",
        ],
    },
    exclude_package_data={
        "mdpow.tests": [
            "testing_resources/*",
        ]
    },
    install_requires=[
        "numpy>=1.6",
        "scipy>=1.11.0",
        "pyyaml",
        "GromacsWrapper>=0.5.1",
        "numkit",
        "six",
        "mdanalysis>=2",
        "alchemlyb>=2",
        "pandas",
        "pymbar>=4",
        "matplotlib",
        "seaborn",
        "rdkit",
        "svgutils",
        "cairosvg",
        "pypdf",
    ],
    # setup_requires=['pytest-runner',],
    # tests_require=["pytest", "pybol", "py"],
    zip_safe=True,
)
