# POW package __init__.py
# Copyright (c) 2010 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""\
MDPOW version information
=========================

MDPOW uses `semantic versioning`_ with the release number consisting
of a triplet *MAJOR.MINOR.PATCH*. *PATCH* releases are bug fixes or
updates to docs or meta data only and do not introduce new features or
change the API. Within a *MAJOR* release, the user API is stable
except during the development cycles with MAJOR = 0 where the API may
also change (rarely) between MINOR releases. *MINOR* releases can
introduce new functionality or deprecate old ones.

Development versions will have the suffix *-dev* after the version
string.

.. _semantic versioning: http://semver.org

Accessing release information
-----------------------------

User code should use :func:`get_version` or `get_version_tuple`.

.. autodata:: VERSION
.. autofunction:: get_version
.. autofunction:: get_version_tuple

"""

#: Package version; this is the only place where it is set.
VERSION = 0,7,0
#: Set to ``True`` for a release. If set to ``False`` then the patch level
#: will have the suffix "-dev".
RELEASE = False
if not RELEASE:
    VERSION = VERSION[:2] + (str(VERSION[2]) + '-dev',)

def get_version():
    """Return current package version as a string."""
    return ".".join(map(str,VERSION))

def get_version_tuple():
    """Return current package version as a tuple (*MAJOR*, *MINOR*, *PATCHLEVEL*)."""
    return tuple(map(str,VERSION))
