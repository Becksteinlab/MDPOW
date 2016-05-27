# restart.py
# Copyright (c) 2010-2011 Oliver Beckstein

"""
:mod:`mdpow.restart` --- Restarting and checkpointing
=====================================================

The module provides classes and functions to keep track of which stages of a
simulation protocol have been completed. It uses a :class:`Journal` class for
the book-keeping. Together with saving the current state of a protocol to a
checkpoint file (using :meth:`Journalled.save`) it is possible to implement
restartable simulation protocols (for example :program:`mdpow-equilibrium`).

.. autoexception:: JournalSequenceError

.. autoclass:: Journal
   :members:

.. autoclass:: Journalled
   :members:

.. autofunction:: checkpoint
"""
from __future__ import absolute_import

import os
import errno
import cPickle

import logging
logger = logging.getLogger('mdpow.checkpoint')

def checkpoint(name, sim, filename):
    """Execute the :meth:`Journalled.save` method and log the event."""
    logger.info("checkpoint: %(name)s", vars())
    sim.save(filename)

class JournalSequenceError(Exception):
    """Raised when a stage is started without another one having been completed."""

class Journal(object):
    """Class that keeps track of the stage in a protocol.

    Transaction blocks have to be bracketed by calls to :meth:`~Journal.start`
    and :meth:`~Journal.completed`. If a block is started before completion, a
    :exc:`JournalSequenceError` will be raised.

    Other methods such as :meth:`~Journal.has_completed` and
    :meth:`~Journal.has_not_completed` can be used to query the status. The
    attribute :attr:`~Journal.incomplete` flags the state of the current stage
    (:attr:`~Journal.current`).

    All completed stages are recorded in the attribute
    :attr:`~Journal.history`.

    The current (incomplete) stage can be reset to its initial state with
    :meth:`Journal.clear`.

    Example::

      J = Journal(['pre', 'main', 'post'])
      J.start('pre')
      ...
      J.completed('pre')
      J.start('main')
      ...
      # main does not finish properly
      print J.incomplete
      # --> 'main'
      J.start('post')
      # raises JournalSequenceError

    """
    def __init__(self, stages):
        """Initialise the journal that keeps track of stages.

        :Arguments:
          *stages*
              list of the stage identifiers, in the order that they
              should per performed. Stage identifiers are checked
              against this list before they are accepted as arguments
              to most methods.
        """
        self.stages = stages  # list of stage identifiers
        self.__current = None
        self.__history = []
        self.__incomplete = None

    @property
    def current(self):
        """Current stage identifier"""
        return self.__current

    @current.setter
    def current(self, stage):
        if not stage in self.stages:
            raise ValueError("Can only assign a registered stage from %r, not %r" %
                             (self.stages, stage))
        self.__current = stage

    @current.deleter
    def current(self):
        self.__current = None

    @property
    def incomplete(self):
        """This last stage was not completed."""
        return self.__incomplete

    @incomplete.setter
    def incomplete(self, stage):
        if not stage in self.stages:
            raise ValueError("can only assign a registered stage from %(stages)r" % vars(self))
        self.__incomplete = stage

    @incomplete.deleter
    def incomplete(self):
        self.__incomplete = None

    @property
    def history(self):
        """List of stages completed"""
        return self.__history

    @history.deleter
    def history(self):
        self.__history = []

    def completed(self, stage):
        """Record completed stage and reset :attr:`Journal.current`"""
        assert stage == self.current, "Program logic error: can only complete the current stage"
        self.__history.append(self.current)
        del self.current

    def start(self, stage):
        """Record that *stage* is starting."""
        if self.current is not None:
            errmsg = "Cannot start stage %s because previously started stage %s " \
                "has not been completed." % (stage, self.current)
            logger.error(errmsg)
            raise JournalSequenceError(errmsg)
        self.current = stage

    def has_completed(self, stage):
        """Returns ``True`` if the *stage* has been started and completed at any time."""
        return stage in self.history

    def has_not_completed(self, stage):
        """Returns ``True`` if the *stage* had been started but not completed yet."""
        return self.current is None and not self.has_completed(stage)

    def clear(self):
        """Reset incomplete status and current stage"""
        del self.incomplete
        del self.current

    def __repr__(self):
        return "%s(%r)" % (self.__class__.__name__, self.stages)

class Journalled(object):
    """A base class providing methods for journalling and restarts.

    It installs an instance of :class:`Journal` in the attribute
    :attr:`Journalled.journal` if it does not exist already.
    """
    #: Class-attribute that contains the names of computation protocols
    #: supported by the class. These are either method names or dummy names,
    #: whose logic is provided by an external callback function.
    #: The method :meth:`get_protocol` raises a :exc:`ValueError` if a
    #: protocol is not listed in :attr:`~Journalled.protocols`.
    protocols = []

    def __init__(self, *args, **kwargs):
        # add journal unless we are starting from a save file that already
        # contains the journal
        try:
            len(self.journal.history)
        except AttributeError:
            self.journal = Journal(self.protocols)
        super(Journalled, self).__init__(*args, **kwargs)

    def get_protocol(self, protocol):
        """Return method for *protocol*.

        - If *protocol* is a real method of the class then the method is
          returned.

        - If *protocol* is a registered protocol name but no method of
          the name exists (i.e. *protocol* is a "dummy protocol") then
          a wrapper function is returned. The wrapper has the
          signature

          .. function:: dummy_protocol(func, *args, **kwargs)

             Runs *func* with the arguments and keywords between calls
             to :meth:`Journal.start` and :meth:`Journal.completed`,
             with the stage set to *protocol*.

        - Raises a :exc:`ValueError` if the *protocol* is not
          registered (i.e. not found in :attr:`Journalled.protocols`).

        """
        if protocol not in self.protocols:
            raise ValueError("%r: protocol must be one of %r" % (protocol, self.protocols))
        try:
            return self.__getattribute__(protocol)
        except AttributeError:
            # catch *_run dummy protocols and have the user provide the function
            return self._journalled_func(protocol)

    def _journalled_func(self, protocol):
        def dummy_protocol(*args, **kwargs):
            """Wrap call to func(args) in journaling."""
            assert len(args) > 0, "f(func, *args, **kwargs) --> func(*args,**kwargs)"
            func = args[0]
            self.journal.start(protocol)
            success = func(*args[1:], **kwargs)
            if success:
                self.journal.completed(protocol)
            return success
        return dummy_protocol

    def save(self, filename=None):
        """Save instance to a pickle file.

        The default filename is the name of the file that was last loaded from
        or saved to. Also sets the attribute :attr:`~Journalled.filename` to
        the absolute path of the saved file.
        """
        if filename is None:
            try:
                if self.filename is not None:
                    filename = self.filename
                else:
                    raise AttributeError
            except AttributeError:
                errmsg = "Neither filename nor default filename provided to save to."
                logger.error(errmsg)
                raise ValueError(errmsg)
        else:
            self.filename = os.path.abspath(filename)
        with open(self.filename, 'wb') as f:
            cPickle.dump(self, f, protocol=cPickle.HIGHEST_PROTOCOL)
        logger.debug("Instance pickled to %(filename)r" % vars(self))

    def load(self, filename=None):
        """Re-instantiate class from pickled file.

        If no *filename* was supplied then the filename is taken from the
        attribute :attr:`~Journalled.filename`.
        """
        if filename is None:
            try:
                if self.filename is not None:
                    filename = self.filename
                else:
                    raise AttributeError
            except AttributeError:
                errmsg = "Neither filename nor default filename provided to load from."
                logger.error(errmsg)
                raise ValueError(errmsg)
        with open(filename, 'rb') as f:
            instance = cPickle.load(f)
        self.__dict__.update(instance.__dict__)
        logger.debug("Instance loaded from %(filename)r" % vars())
