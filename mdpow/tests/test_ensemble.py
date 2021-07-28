from __future__ import absolute_import

from . import tempdir

import os.path
import sys

import pytest
import py.path

import yaml
import pybol

from six.moves import cPickle as pickle

from mdpow.analysis import Ensemble
import mdpow.fep

from pkg_resources import resource_filename
RESOURCES = py.path.local(resource_filename(__name__, 'testing_resources'))
MANIFEST = RESOURCES.join("manifest.yml")

ensemble_keys = [('water', 'Coulomb', '0000'),
                 ('water', 'Coulomb', '0250'),
                 ('water', 'Coulomb', '0500'),
                 ('water', 'Coulomb', '0750'),
                 ('water', 'Coulomb', '1000'),
                 ('water', 'VDW', '0050'),
                 ('water', 'VDW', '0100'),
                 ('water', 'VDW', '0200'),
                 ('water', 'VDW', '0300'),
                 ('water', 'VDW', '0400'),
                 ('water', 'VDW', '0500'),
                 ('water', 'VDW', '0600'),
                 ('water', 'VDW', '0650'),
                 ('water', 'VDW', '0750'),
                 ('water', 'VDW', '0800'),
                 ('water', 'VDW', '0850'),
                 ('water', 'VDW', '0900'),
                 ('water', 'VDW', '0950'),
                 ('water', 'VDW', '1000')]


def fix_manifest(topdir):
    """Create a temporary manifest with a custom `path`.

    Fix manifest in a local temporary copy in existing dir topdir
    where the `path` is an absolute path to our "states" directory. We
    use `pkg_resources.resource_filename` to anchor the path.

    Arguments
    ---------
    topdir : py.path.local
        existing temporary directory (as provided by, for instance,
        `pytest.tmpdir`)

    Returns
    -------
    new_manifest : py.path.local
        Path to the new manifest.

    Example
    -------
    Use as ::

        new_manifest = fix_manifest(tmpdir)
        m = pybol.Manifest(new_manifest.strpath)

    """
    manifest = yaml.safe_load(MANIFEST.open())
    # simple heuristic: last element of the recorded manifest::path is the name
    # of the states directory, typically 'states' (from .../testing_resources/states)
    manifest['path'] = RESOURCES.join(os.path.basename(manifest['path'])).strpath
    new_manifest = topdir.join("local_manifest.yml")
    yaml.dump(manifest, stream=new_manifest.open("w"))
    return new_manifest


# session scope if read-only use

@pytest.fixture(scope="function")
def fep_benzene_directory(tmpdir_factory):
    topdir = tmpdir_factory.mktemp('ensemble')
    m = pybol.Manifest(fix_manifest(topdir).strpath)
    m.assemble('FEP', topdir.strpath)
    return topdir.join("benzene")


class TestAnalyze(object):
    def get_Gsolv(self, pth):
        gsolv = pth.join("FEP", "water", "Gsolv.fep")
        # Needed to load old pickle files in python 3
        if sys.version_info.major >= 3:
            with open(gsolv, 'rb') as f:
                G = pickle.load(f, encoding='latin1')
                # patch paths
        elif sys.version_info.major == 2:
            G = pickle.load(gsolv.open())
        G.basedir = pth.strpath
        G.filename = gsolv.strpath
        return G

    def test_build_ensemble(self):
        G = self.get_Gsolv(fep_benzene_directory)
        G.convert_edr()
        bnz = Ensemble(dirname=fep_benzene_directory, solvents=['water'])
        assert bnz.get_keys() == ensemble_keys
