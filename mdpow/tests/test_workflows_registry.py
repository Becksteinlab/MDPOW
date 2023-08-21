import pytest

from mdpow.workflows import registry


def test_registry():
    assert list(registry.registry.keys()) == ["DihedralAnalysis"]
