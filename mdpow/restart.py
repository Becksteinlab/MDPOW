# restart.py
# Copyright (c) 2010-2011 Oliver Beckstein

"""Support for restarting protocols and checkpointing"""

import os, errno

import logging
logger = logging.getLogger('mdpow.checkpoint')

def checkpoint(name, sim, filename):
    logger.info("checkpoint: %(name)s", vars())
    sim.save(filename)

class JournalSequenceError(Exception):
    """Raised when a stage is started without another one having been completed."""

class Journal(object):
    """Class that keeps track of the stage in protocol.

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

    def current():
        doc = """Current stage identifier"""
        def fget(self):
            return self.__current
        def fset(self, stage):
            if not stage in self.stages:
                raise ValueError("Can only assign a registered stage from %r, not %r" %
                                 (self.stages, stage))
            self.__current = stage
        def fdel(self):
            self.__current = None
        return locals()
    current = property(**current())

    def incomplete():
        doc = """This last stage was not completed."""
        def fget(self):
            return self.__incomplete
        def fset(self, stage):
            if not stage in self.stages:
                raise ValueError("can only assign a registered stage from %(stages)r" % vars(self))
            self.__incomplete = stage
        def fdel(self):
            self.__incomplete = None
        return locals()
    incomplete = property(**incomplete())

    # maybe make this a ringbuffer
    def history():
        doc = """List of stages completed"""
        def fget(self):
            return self.__history
        def fdel(self):
            self.__history = []
        return locals()
    history = property(**history())

    def completed(self, stage):
        """Record completed stage and reset current"""
        assert stage == self.current, "Program logic error: can only complete the current stage"
        self.__history.append(self.current)
        del self.current

    def start(self, stage):
        """Record that we are starting a stage"""
        if not self.current is None:
            errmsg = "Cannot start stage %s because previously started stage %s " \
                "has not been completed." % (stage, self.current)
            logger.error(errmsg)
            raise JournalSequenceError(errmsg)
        self.current = stage

    def has_completed(self, stage):
        return stage in self.history

    def has_not_completed(self, stage):
        return self.current is None and not self.has_completed(stage)

    def clear(self):
        """Reset incomplete status and current stage"""
        del self.incomplete
        del self.current
