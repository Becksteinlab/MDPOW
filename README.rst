=================== 
 README for MDPOW
=================== 

|build| |cov|

.. |P_ow| replace:: *P*\ :sub:`OW`
.. |P_cw| replace:: *P*\ :sub:`CW`

*MDPOW* is a python package that automates the calculation of
solvation free energies via molecular dynamics (MD) simulations. In
particular, it facilitates the computation of partition
coeffcients. Currently implemented:

- *water-octanol* partition coefficient (|P_ow|)
- *water-cyclohexane* partition coefficient (|P_cw|)

Calculations are performed with the Gromacs_ MD software package
[#GromacsWrapperNote]_. Currently, *OPLS-AA force field* parameters are
supported.

As *input*, the user only needs to provide a structure file (PDB or
GRO) and a Gromacs ITP file containing the parametrization of the
small molecule (e.g. from LigandBook_ or ParamChem_).

.. _Gromacs: http://www.gromacs.org
.. _GromacsWrapper: http://gromacswrapper.readthedocs.org/en/latest/
.. _LigandBook: http://ligandbook.icsn.cnrs-gif.fr/
.. _ParamChem: https://cgenff.paramchem.org/



Installation
------------

Install from the checked out source::

  git clone https://github.com/Becksteinlab/MDPOW.git  
  pip install MDPOW/

This will also install additional dependencies such as GromacsWrapper_. You
will also need `Gromacs`_ (currently tested with versions 4.6.7,
5.1.2, and Gromacs 2016).


Source code
-----------

*MDPOW* is open source and published under the `GNU General Public License
v3`_. Source code is available at https://github.com/Becksteinlab/MDPOW .

.. _`GNU General Public License v3`: 
   http://www.gnu.org/licenses/gpl-3.0.html

Footnotes
---------

.. [#GromacsWrapperNote] The package is built on top of the GromacsWrapper_
                         framework (which is automatically installed).

.. |build| image:: https://travis-ci.org/Becksteinlab/MDPOW.svg?branch=develop
   :alt: Build Status
   :target: https://travis-ci.org/Becksteinlab/MDPOW

.. |cov| image:: https://codecov.io/github/Becksteinlab/MDPOW/coverage.svg?branch=develop
   :alt: Coverage Status
   :target: https://codecov.io/github/Becksteinlab/MDPOW?branch=develop

   
