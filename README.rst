===================
 README for MDPOW
===================

|build| |cov| |docs| |black| |zenodo|

.. |P_ow| replace:: *P*\ :sub:`OW`
.. |P_cw| replace:: *P*\ :sub:`CW`
.. |P_tw| replace:: *P*\ :sub:`TW`  

*MDPOW* is a python package that automates the calculation of
solvation free energies via molecular dynamics (MD) simulations. In
particular, it facilitates the computation of partition
coefficients. Currently implemented:

- *water-octanol* partition coefficient (|P_ow|)
- *water-cyclohexane* partition coefficient (|P_cw|)
- *water-toluene* partition coefficient (|P_tw|)
  
Calculations are performed with the Gromacs_ MD software package
[#GromacsWrapperNote]_. Currently, *OPLS-AA*, *CHARMM/CGENFF*, and
*AMBER/GAFF* parameters are supported.

As *input*, the user only needs to provide a structure file (PDB or
GRO) and a Gromacs ITP file containing the parametrization of the
small molecule (e.g. from LigandBook_ or ParamChem_).

.. _Gromacs: http://www.gromacs.org
.. _GromacsWrapper: http://gromacswrapper.readthedocs.org/en/latest/
.. _LigandBook: http://ligandbook.org/
.. _ParamChem: https://cgenff.paramchem.org/


Documentation
-------------

* https://mdpow.readthedocs.io
* `Tutorial`_ : computing the octanol-water partition coefficient of
  benzene (uses the `example files`_)


.. _Tutorial: http://mdpow.readthedocs.io/en/latest/init.html#tutorial-using-the-mdpow-scripts-to-compute-logpow-of-benzene
.. _example files: https://github.com/Becksteinlab/MDPOW/tree/develop/doc/examples

Installation
------------

See `INSTALL`_ for detailed instructions. MDPOW currently supports and
is tested with Python 3.10 to 3.12.

You will also need `Gromacs`_ (currently tested with versions 4.6.5,
2018, 2020, 2021, 2022, 2023, 2024 but 2016 and 2019 should also work).


Development version
~~~~~~~~~~~~~~~~~~~

If you want to install the development version, get the sources from
GitHub (the development branch) ::

  git clone https://github.com/Becksteinlab/MDPOW.git

and Install from the checked out source::

  pip install MDPOW/

(Note the trailing slash ``/`` to indicate the directory.)



Source code
-----------

*MDPOW* is open source and published under the `GNU General Public License
v3`_. Source code is available at https://github.com/Becksteinlab/MDPOW .

We use `black`_ for uniform code formatting.

.. _`GNU General Public License v3`:
   http://www.gnu.org/licenses/gpl-3.0.html

.. _`black`: https://github.com/psf/black


Footnotes
---------

.. [#GromacsWrapperNote] The package is built on top of the GromacsWrapper_
                         framework (which is automatically installed).

.. |build| image:: https://github.com/Becksteinlab/MDPOW/actions/workflows/ci.yaml/badge.svg?branch=develop
   :alt: Build Status
   :target: https://github.com/Becksteinlab/MDPOW/actions/workflows/ci.yaml

.. |cov| image:: https://codecov.io/github/Becksteinlab/MDPOW/coverage.svg?branch=develop
   :alt: Coverage Status
   :target: https://codecov.io/github/Becksteinlab/MDPOW?branch=develop

.. |docs| image:: https://readthedocs.org/projects/mdpow/badge/?version=latest
   :target: http://mdpow.readthedocs.org/en/latest/?badge=latest
   :alt: Documentation
   
.. |zenodo| image:: https://zenodo.org/badge/44999898.svg
   :target: https://zenodo.org/badge/latestdoi/44999898
   :alt: Zenodo

.. |black| image:: https://img.shields.io/badge/code%20style-black-000000.svg
   :target: https://github.com/psf/black	 
   :alt: black   

.. _INSTALL: INSTALL.rst
