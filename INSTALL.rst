=============================================
 Quick installation instructions for *MDPOW*
=============================================

Releases
--------

The latest release can be installed from the internet with ::

  pip install --upgrade MDPOW


Standard installation from source 
---------------------------------

Example for installation as a user from the checked out source (by
default, uses the development branch)::

  git clone https://github.com/Becksteinlab/MDPOW.git  
  pip install --user MDPOW/

Check that you can import the module::

  python
  >>> import mdpow
  >>> help(mdpow)

In case of problems  file an issue at
https://github.com/Becksteinlab/MDPOW/issues


Developer installation
----------------------

A development install is useful while hacking away on the code::

 cd MDPOW
 pip install -e --user .

  
