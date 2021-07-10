# POW package __init__.py
# Copyright (c) 2010 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.


from __future__ import absolute_import

from .version import VERSION, get_version, get_version_tuple
from . import log

__all__ = ['fep', 'equil']


def create_logger(logfile="mdpow.log"):
    """Create the default logger.

    Channels the output from :mod:`mdpow`, :mod:`gromacs`, and
    :mod:`numkit` into the file *logfile*.
    """
    logger = log.create('mdpow', logfile)
    log.create('numkit', logfile)   # capture numkit messages to same file
    log.create('gromacs', logfile)  # and the GromacsWrapper messages
    return logger

def log_banner():
    """Log program name and licence at INFO level."""
    logger.info("MDPOW %s starting.", get_version())
    logger.info("Copyright (c) 2010-2021 Shujie Fan, Ian Kenney, Bogdan Iorga, and Oliver Beckstein")
    logger.info("Released under the GNU Public Licence, version 3.")
    logger.info("For bug reports and help: https://github.com/Becksteinlab/MDPOW/issues")

logger = create_logger()
log_banner()

# AVOID IMPORTS OF OTHER PACKAGES IN __init__.py; only standard
# library and mdpow.log are allowed. Anything else will break 'pip
# install' because in setup.py we import mdpow.version (and thus
# __init__.py) for the dynamic version information BEFORE pip has a
# chance to install dependencies.

# XXX move to tables XXX
#: Avogadro's constant |NA| in mol^-1 (`NA NIST value`_).
N_AVOGADRO = 6.02214179e23
#: Boltzmann's constant |kB| in kJ mol^-1 (`kB NIST value`_).
kBOLTZ = 1.3806504e-23 *1e-3 * N_AVOGADRO

