=============================================
 Quick installation instructions for *MDPOW*
=============================================

**Only Python 2.7 is supported** (Python 3 is *not* supported, see
`#84 <https://github.com/Becksteinlab/MDPOW/issues/84>`_).  Python 2.7
is rock-stable and frozen but not officially supported anymore by the
Python developers.

We recommend that you install MDPOW in a virtual environment.


Installation for users
======================

**Most users should follow these instructions.**

We use the Anaconda distribution with the ``conda`` command to manage
dependencies (see `installing the Anaconda distribution
<https://docs.anaconda.com/anaconda/install/>`_ for details on how to
set up ``conda``).

To run MDPOW you will also need to install a compatible version of
GROMACS_.

.. _GROMACS: http://www.gromacs.org



Conda environment with pre-requisites
-------------------------------------

We make a conda environment with the latest packages for Python 2.7
and name it *mdpow*; this installs the larger dependencies that are
pre-requisites for MDPOW::

 conda create -c conda-forge -n mdpow python=2.7 numpy scipy 'matplotlib<3.3' 'mdanalysis<2' 'mdanalysistests<2' pyyaml
 conda activate mdpow  
 pip install gromacswrapper

Installation from source
------------------------
 
Install from source::

 conda activate mdpow
 git clone https://github.com/Becksteinlab/MDPOW.git
 pip install ./MDPOW

(Older releases are available but outdated; use the latest source for now.)


Checking that the installation worked
-------------------------------------

Check that you can run the ``mdpow-*`` commandline tools::

  mdpow-equilibrium --help

should print a whole bunch of messages. If your GROMACS_ installation
cannot be found, make sure you ``source GMXRC`` or load the appopriate
modules or whatever else you have to do so that the ``gmx`` command is
found. See https://gromacswrapper.readthedocs.io for more details.


Check that you can import the module::

  python
  >>> import mdpow
  >>> help(mdpow)

In case of problems  file an issue at
https://github.com/Becksteinlab/MDPOW/issues




Developer installation
======================

A development install is useful while hacking away on the code::

 cd MDPOW
 pip install -e  .

  
