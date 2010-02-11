# config.py
# Copyright (c) 2010 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""
:mod:`pow.config` -- Configuration for POW
==========================================================

The config module provides configurable options for the whole package;
eventually it might grow into a sophisticated configuration system such as
matplotlib's rc system but right now it mostly serves to define which gromacs
tools and other scripts are offered in the package and where template files are
located. If the user wants to change anything they will still have to do it
here in source until a better mechanism with rc files has been implemented.


Logging
-------

Gromacs commands log their invocation to a log file; typically at
loglevel *INFO* (see the python `logging module`_ for details).

.. autodata:: logfilename
.. autodata:: loglevel_console
.. autodata:: loglevel_file

.. _logging module: http://docs.python.org/library/logging.html



Location of template files
--------------------------

*Template variables* list files in the package that can be used as
templates such as run input files. Because the package can be a zipped
egg we actually have to unwrap these files at this stage but this is
completely transparent to the user.

.. autodata:: templates
.. autodata:: sge_template


Functions
---------

The following functions can be used to access configuration data.

.. autofunction:: get_template
.. autofunction:: get_templates

"""

import os
from pkg_resources import resource_filename

import utilities



# Location of template files
# --------------------------

# TODO: This is becoming unwieldy: should be more abstract/automatic.
#       Maybe simply use filename as key? This would require some cleaning up
#       in gromacs.setup but would not be too bad.
templates = {
    'md_OPLSAA_mdp': resource_filename(__name__, 'templates/md_OPLSAA.mdp'),
    }
"""Templates have to be extracted from the egg because they are used
by external code. All template filenames are stored in
:data:`gromacs.config.templates`.

**Gromacs mdp templates**

   These are supplied as examples and there is *NO GUARANTEE THAT THEY
   PRODUCE SENSIBLE OUTPUT* --- check for yourself!  Note that only
   existing parameter names can be modified with
   :func:`gromacs.cbook.edit_mdp` at the moment; if in doubt add the
   parameter with its gromacs default value (or empty values) and
   modify later with :func:`~gromacs.cbook.edit_mdp`.

"""

# Functions to access configuration data
# --------------------------------------

def get_template(t):
    """Find template file *t* and return its real path.

    *t* can be a single string or a list of strings. A string
    should be one of

    1. a relative or absolute path,
    2. a filename in the package template directory (defined in the template dictionary
       :data:`gromacs.config.templates`) or
    3. a key into :data:`~gromacs.config.templates`.

    The first match (in this order) is returned. If the argument is a
    single string then a single string is returned, otherwise a list
    of strings.

    :Arguments: *t* : template file or key (string or list of strings)
    :Returns:   os.path.realpath(*t*) (or a list thereof)
    :Raises:    :exc:`ValueError` if no file can be located.
       
    """
    templates = [_get_template(s) for s in asiterable(t)]
    if len(templates) == 1:
         return templates[0]
    return templates

def get_templates(t):
    """Find template file(s) *t* and return their real paths.

    *t* can be a single string or a list of strings. A string should
    be one of

    1. a relative or absolute path,
    2. a filename in the package template directory (defined in the template dictionary
       :data:`gromacs.config.templates`) or
    3. a key into :data:`~gromacs.config.templates`.

    The first match (in this order) is returned for each input argument.

    :Arguments: *t* : template file or key (string or list of strings)
    :Returns:   list of os.path.realpath(*t*) 
    :Raises:    :exc:`ValueError` if no file can be located.
       
    """
    return [_get_template(s) for s in utilities.asiterable(t)]

def _get_template(t):
    """Return a single template *t*."""
    if os.path.exists(t):           # 1) Is it an accessible file?
        pass
    else:                           # 2) check the packaged template files
        _t = os.path.basename(t)
        _t_found = False
        for p in templates.values():
            if _t == os.path.basename(p):
                t = p
                _t_found = True     # NOTE: in principle this could match multiple
                break               #       times if more than one template dir existed.
        if not _t_found:            # 3) try it as a key into templates
            try:
                t = templates[t]
            except KeyError:
                pass
            else:
                _t_found = True
        if not _t_found:            # 4) nothing else to try... or introduce a PATH?
            raise ValueError("Failed to locate the template file %(t)r." % vars())
    return os.path.realpath(t)


# utility functions  (from gromacs.utilities)
# Copied so that config does not have a dependency on gromacs.utilities

def iterable(obj):
    """Returns ``True`` if *obj* can be iterated over and is *not* a  string."""
    if type(obj) is str:
        return False    # avoid iterating over characters of a string

    if hasattr(obj, 'next'):
        return True    # any iterator will do 
    try: 
        len(obj)       # anything else that might work
    except TypeError: 
        return False
    return True

def asiterable(obj):
    """Returns obj so that it can be iterated over; a string is *not* treated as iterable"""
    if not iterable(obj):
        obj = [obj]
    return obj
