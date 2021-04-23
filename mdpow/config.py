# POW: config.py
# Copyright (c) 2010-2011 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""
:mod:`mdpow.config` -- Configuration for MDPOW
==============================================

The config module provides configurable options for the whole package;
eventually it might grow into a more sophisticated configuration system but
right now it mostly serves to define which gromacs tools and other scripts are
offered in the package and where template files are located. If the user wants
to change anything they will still have to do it here in source until a better
mechanism with a global configuration file has been implemented.


Force field
-----------

By default, MDPOW uses a collection of OPLS/AA force field files based on the
Gromacs 4.5.3 distribution, with the following differences:

* For ions we use the new alkali and halide ion parameters from Table 2 in
  [Jensen2006]_ which had shown some small improvements in the paper. They
  should only be used with the TIP4P water model.

* OPLS/AA parameters for 1-octanol were added. These parameters were validated
  against experimental data by computing the density (neat), hydration free
  energy and logP (the latter being a self consistency check).

  .. TODO add the results of the checks

The force field files are found in the directory pointed to by the environment
variable :envvar:`GMXLIB`. By default, :mod:`mdpow.config` sets
:envvar:`GMXLIB` to :data:`includedir` unless :envvar:`GMXLIB` has already been
set. This mechanism allows the user to override the choice of location of force
field.

At the moment, only OPLS/AA is tested with MDPOW although in principle it is
possible to use other force fields by supplying appropriately customized
template files.

.. rubric:: References

.. [Jensen2006] K.P. Jensen and W.L. Jorgensen, *J Comp Theor Comput* **2**
                (2006), 1499.  doi:`10.1021/ct600252r`_

.. _`10.1021/ct600252r`: http://dx.doi.org/10.1021/ct600252r



Location of template files
--------------------------

*Template variables* list files in the package that can be used as
templates such as run input files. Because the package can be a zipped
egg we actually have to unwrap these files at this stage but this is
completely transparent to the user.

.. autodata:: templates
.. autodata:: topfiles
.. autodata:: includedir
.. autodata:: defaults

Functions
---------

The following functions can be used to access configuration data.

.. autofunction:: get_template
.. autofunction:: get_templates
.. autofunction:: get_configuration

.. rubric:: Developer note

Templates have to be extracted from the egg because they are used by external
code. All template filenames are stored in :data:`config.templates` or
:data:`config.topfiles`.

Sub-directories are extracted (see `Resource extraction`_) but the file names
themselves will not appear in the template dict. Thus, only store files in
subdirectories that don't have to be explicitly found by the package (e.g. the
Gromacs force field files are ok).

.. _Resource extraction:
   http://packages.python.org/distribute/pkg_resources.html#resource-extraction

.. autofunction:: _generate_template_dict

"""

from __future__ import absolute_import

import os, errno
from pkg_resources import resource_filename, resource_listdir
import yaml

import numpy as np

import logging
logger = logging.getLogger("mdpow.config")


# Reading of configuration files
# ------------------------------

#: Locations of default run input files and configurations.
defaults = {
    "runinput": resource_filename(__name__, "templates/runinput.yml"),
    }

def merge_dicts(user, default):
    """Merge two dictionaries recursively.

    Based on https://stackoverflow.com/a/823240/334357
    """
    if isinstance(user, dict) and isinstance(default, dict):
        for k, v in default.iteritems():
            if k not in user:
                user[k] = v
            else:
                user[k] = merge_dicts(user[k], v)
    return user


class POWConfigParser(object):
    """Parse YAML config file."""

    def __init__(self):
        self.conf = None

    def readfp(self, fn):
        """Read YAML from open stream ``fn``.

        Overwrites everything.
        """
        self.conf = yaml.safe_load(fn)
        return True

    def merge(self, fn):
        """Load YAML from open stream ``fn`` and merge into :attr:`conf`.

        Data from this file will overwrite data from the existing
        configuration; anything not set in this file will be taken
        from the loaded configuration. Arrays are overwritten and not
        appended/merged.
        """
        user = yaml.safe_load(fn)
        self.conf = merge_dicts(user, self.conf)
        return self.conf

    def write(self, filename):
        with open(filename, 'w') as f:
            f.write(yaml.dump(self.conf))

    def get(self, section, option):
        """Return option, unless its "None" --> ``None``,

        Conversion to basic python types str, float, int, boolean is
        carried out automatically (unless it was None).

        .. Note:: "none" remains a string, which is essential, see
                  `Issue 20 <https://github.com/Becksteinlab/MDPOW/issues/20>`_

        .. versionchanged:: 0.6.0
           Prior versions would convert case-insensitively (e.g. "NONE"
           and "none")
        """
        value = self.conf[section][option]
        return value if value != "None" else None

    # TODO:
    # The YAML parser does automatic conversion: the following
    # methods are for backward compatibility with the old ini parser
    # and should be cleaned up. --- orbeckst 2016-01-18
    getstr = get
    getfloat = get
    getint = get
    getboolean = get

    def getpath(self, section, option):
        """Return option as an expanded path."""
        return os.path.expanduser(os.path.expandvars(self.get(section, option)))

    def findfile(self, section, option):
        """Return location of a file ``option``.

        Uses :func:`mdpow.config.get_template`.
        """
        return get_template(self.getpath(section, option))

    # TODO: Change input file format to use yaml lists and make this method superfluous
    def getlist(self, section, option):
        """Return option as a list of strings.

        *option* must be comma-separated; leading/trailing whitespace
        is stripped and quotes are treated verbatim.
        """
        return [x.strip() for x in str(self.get(section, option)).split(",")]

    def getarray(self, section, option):
        """Return option as a numpy array of floats.

        *option* must be comma-separated; leading/trailing whitespace
        is stripped and quotes are treated verbatim.
        """
        return np.asarray(self.getlist(section, option), dtype=np.float)

    def getintarray(self, section, option):
        """Return option as a numpy array of integers.

        *option* must be comma-separated; leading/trailing whitespace
        is stripped and quotes are treated verbatim.
        """
        return np.asarray(self.getlist(section, option), dtype=np.int)

def get_configuration(filename=None):
    """Reads and parses a run input config file.

    Uses the package-bundled defaults as a basis.
    """
    cfg = POWConfigParser()
    cfg.readfp(open(defaults["runinput"]))
    logger.debug("Loaded runinput defaults from %r", defaults["runinput"])
    if filename is not None:
        cfg.merge(open(filename))   # override package defaults
        logger.debug("Loaded user runinput from %r (replacing defaults)", filename)
    else:
        logger.warning("Running with package defaults for the run; you should supply a runinput file!")
    return cfg

def modify_gromacs_environment(name, value):
    from gromacs.environment import flags
    if flags[name] != value:
        logger.warn("Changing GromacsWrapper environment: flags[%(name)r] = %(value)r", vars())
        flags[name] = value

def set_gromacsoutput(cfg):
    # maybe allow setting this on a per protocol basis?
    modify_gromacs_environment('capture_output', not cfg.getboolean('setup', 'gromacsoutput'))


# Functions to locate template files
# ----------------------------------

def _generate_template_dict(dirname):
    """Generate a list of included top-level files *and* extract them to a temp space.

    Only lists files and directories at the *top level* of the *dirname*;
    however, all directories are extracted recursively and will be available.
    """
    return dict((resource_basename(fn), resource_filename(__name__, dirname+'/'+fn))
                for fn in resource_listdir(__name__, dirname)
                if not fn.endswith('~'))

def resource_basename(resource):
     """Last component of a resource (which always uses '/' as sep)."""
     if resource.endswith('/'):
          resource = resource[:-1]
     parts = resource.split('/')
     return parts[-1]


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
    if isinstance(obj, basestring):
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


# Setting up configuration variables and paths
#---------------------------------------------

templates = _generate_template_dict('templates')
"""*POW* comes with a number of templates for run input files
and queuing system scripts. They are provided as a convenience and
examples but **WITHOUT ANY GUARANTEE FOR CORRECTNESS OR SUITABILITY FOR
ANY PURPOSE**.

All template filenames are stored in
:data:`gromacs.config.templates`. Templates have to be extracted from
the GromacsWrapper python egg file because they are used by external
code: find the actual file locations from this variable.

**Gromacs mdp templates**

   These are supplied as examples and there is *NO GUARANTEE THAT THEY
   PRODUCE SENSIBLE OUTPUT* --- check for yourself!  Note that only
   existing parameter names can be modified with
   :func:`gromacs.cbook.edit_mdp` at the moment; if in doubt add the
   parameter with its gromacs default value (or empty values) and
   modify later with :func:`~gromacs.cbook.edit_mdp`.


   The safest bet is to use one of the ``mdout.mdp`` files produced by
   :func:`gromacs.grompp` as a template as this mdp contains all
   parameters that are legal in the current version of Gromacs.
"""

#: List of all topology files that are included in the package.
#: (includes force field files under ``top/oplsaa.ff``)
topfiles = _generate_template_dict('top')
topfiles.update(_generate_template_dict('top/oplsaa.ff'))  # added manually!

# Find the top include dir by looking for an important file 'ffoplsaa.itp'.
# Since Gromacs 4.5.x, force fields are ONLY found in
# 1) the current directory
# 2) the directory pointed to by the environment variable GMXLIB
# 3) or the default Gromacs installation top directory
try:
    #: The package's include directory for :func:`gromacs.grompp`; the
    #: environment variable :envvar:`GMXLIB` is set to :data:`includedir`
    #: so that the bundled version of the force field is picked up.
    includedir = os.path.dirname(topfiles['ffoplsaa.itp'])
except KeyError:
    errmsg = "Missing required data files (ffoplsaa.itp). Check your installation."
    logger.fatal(errmsg)
    raise ImportError(errmsg)

if not 'GMXLIB' in os.environ:
    if not os.path.exists(includedir):
        errmsg = "Likely installation problem: cannot access the package GMXLIB " \
            "directory (try re-installing): "
        logger.fatal(errmsg + includedir)
        raise OSError(errno.ENOENT, errmsg, includedir)
    os.environ['GMXLIB'] = includedir
    logger.info("Using the bundled force fields from GMXLIB=%(includedir)r.", vars())
    logger.info("If required, override this behaviour by setting the environment variable GMXLIB yourself.")
else:
    logger.warn("Using user-supplied environment variable GMXLIB=%r to find force fields", os.environ['GMXLIB'])
    logger.info("(You can use the MDPOW default by executing 'unset GMXLIB' in your shell before running MDPOW.)")

