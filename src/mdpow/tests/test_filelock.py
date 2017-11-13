import os.path

import pytest

import mdpow.filelock

def test_FileLock_acquire(tmpdir, filename="test.txt"):
    with tmpdir.as_cwd():
        with mdpow.filelock.FileLock(filename, timeout=2) as lock:
            with open(filename, "w") as f:
                f.write("Humpty Dumpty sat on a wall")
        assert os.path.exists(filename)

def test_FileLock_lock(filename="test.txt"):
    with mdpow.filelock.FileLock(filename, timeout=2) as lock:
        with pytest.raises(mdpow.filelock.FileLockException):
            with mdpow.filelock.FileLock(filename, timeout=0.1) as lock2:
                pass
