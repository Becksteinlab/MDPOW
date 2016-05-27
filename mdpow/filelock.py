"""
Cross-platform filelocking
==========================

:Author: Evan Fosmark
:Year: 2009
:License: BSD (as mentioned in Evan's comment dating 22 April 2009 1:04pm)
:URL: http://www.evanfosmark.com/2009/01/cross-platform-file-locking-support-in-python/

On occasion, one requires the need to lock a file. Now, this is relatively easy
if you're targeting a specific platform because there is often a function in
the library to do it for you. But what if you want to target a larger set of
platforms? The following is a solution I wrote up today. It's lockfile creation
is an atomic operation and thus doesn't suffer from any race conditions. It
should work in both Windows and Unix environments.

The above class is best used in a context manager fashion through the with statement like in the example below::

  with FileLock("test.txt", timeout=2) as lock:
      print("Lock acquired.")
      # Do something with the locked file

The largest downside of this is that the directory the file is located in must
be writable. I hope this code helps you. Of course, if you have a better
recipe, please share it in the comments. ;)

"""

from __future__ import absolute_import

import os
import time
import errno

class FileLockException(Exception):
    pass

class FileLock(object):
    """ A file locking mechanism that has context-manager support so
        you can use it in a with statement. This should be relatively cross
        compatible as it doesn't rely on msvcrt or fcntl for the locking.
    """

    def __init__(self, file_name, timeout=10, delay=.05):
        """ Prepare the file locker. Specify the file to lock and optionally
            the maximum timeout and the delay between each attempt to lock.
        """
        self.is_locked = False
        self.lockfile = os.path.join(os.getcwd(), "%s.lock" % file_name)
        self.file_name = file_name
        self.timeout = timeout
        self.delay = delay


    def acquire(self):
        """ Acquire the lock, if possible. If the lock is in use, it check again
            every `wait` seconds. It does this until it either gets the lock or
            exceeds `timeout` number of seconds, in which case it throws
            an exception.
        """
        start_time = time.time()
        while True:
            try:
                self.fd = os.open(self.lockfile, os.O_CREAT|os.O_EXCL|os.O_RDWR)
                break;
            except OSError as e:
                if e.errno != errno.EEXIST:
                    raise
                if (time.time() - start_time) >= self.timeout:
                    raise FileLockException("Timeout occured.")
                time.sleep(self.delay)
        self.is_locked = True


    def release(self):
        """ Get rid of the lock by deleting the lockfile.
            When working in a `with` statement, this gets automatically
            called at the end.
        """
        if self.is_locked:
            os.close(self.fd)
            os.unlink(self.lockfile)
            self.is_locked = False


    def __enter__(self):
        """ Activated when used in the with statement.
            Should automatically acquire a lock to be used in the with block.
        """
        if not self.is_locked:
            self.acquire()
        return self


    def __exit__(self, type, value, traceback):
        """ Activated at the end of the with statement.
            It automatically releases the lock if it isn't locked.
        """
        if self.is_locked:
            self.release()


    def __del__(self):
        """ Make sure that the FileLock instance doesn't leave a lockfile
            lying around.
        """
        self.release()
