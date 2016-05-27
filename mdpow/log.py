# log.py
# Copyright (c) 2010 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""
:mod:`mdpow.log` --- Configure logging for POW analysis
=======================================================

Import this module if logging is desired in application code and
create the logger in ``__init__.py``::

  import log
  logger = log.create(logname, logfile)

In modules simply use::

  import logging
  logger = logging.getLogger(logname)
"""

from __future__ import absolute_import

# log is the only package that is imported in __init__ so it MAY NOT
# HAVE ANY DEPENDENCIES except standard library; in particular, DO NOT
# 'from . import config' or similar because it will break pip install
# (because in setup.py we import mdpow.version (and thus __init__.py)
# for the dynamic version information)

import logging

def create(logname, logfile):
    """Create a top level logger.

    - The file logger logs everything (including DEBUG).
    - The console logger only logs INFO and above.
    """

    logger = logging.getLogger(logname)

    logger.setLevel(logging.DEBUG)

    logfile = logging.FileHandler(logfile)
    logfile_formatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
    logfile.setFormatter(logfile_formatter)
    logger.addHandler(logfile)

    # define a Handler which writes INFO messages or higher to the sys.stderr
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    # set a format which is simpler for console use
    formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
    console.setFormatter(formatter)
    logger.addHandler(console)

    return logger

def clear_handlers(logger):
    """clean out handlers in the library top level logger

    (only important for reload/debug cycles...)
    """
    for h in logger.handlers:
        logger.removeHandler(h)


