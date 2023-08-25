# tests for MDPOW

import py.path
from pkg_resources import resource_filename


RESOURCES = py.path.local(resource_filename(__name__, "testing_resources"))

MANIFEST = RESOURCES / "manifest.yml"

MOLECULES = {
    "benzene": RESOURCES.join("molecules", "benzene"),
}
STATES = {
    "FEP": RESOURCES.join("states", "FEP"),
    "base": RESOURCES.join("states", "base"),
    "md_npt": RESOURCES.join("states", "FEP"),
    "workflows": RESOURCES.join("states", "workflows"),
}
CONFIGURATIONS = RESOURCES.join("test_configurations")
