import pytest

import mdpow

@pytest.fixture(scope="module")
def version():
    return mdpow.__version__

def test_version_string(version):
    assert isinstance(version, six.string_types)

def test_version(version):
    # generic non-empty check because versioneer can provide different
    # answers depending on VCS status
    assert version


