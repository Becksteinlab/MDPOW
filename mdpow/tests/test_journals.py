import pytest

from mdpow import restart


@pytest.fixture
def journal():
    return restart.Journal(["pre", "main", "post"])


class TestJournal:
    def test_full_sequence(self, journal):
        journal.start("pre")
        assert journal.current == "pre"
        journal.completed("pre")

        journal.start("main")
        assert journal.current == "main"
        journal.completed("main")

        journal.start("post")
        assert journal.current == "post"
        journal.completed("post")

    def test_set_wrong_stage_ValueError(self, journal):
        with pytest.raises(ValueError, match="Can only assign a registered stage"):
            journal.start("BEGIN !")

    def test_JournalSequenceError_no_completion(self, journal):
        with pytest.raises(restart.JournalSequenceError, match="Cannot start stage"):
            journal.start("pre")
            assert journal.current == "pre"

            journal.start("main")

    @pytest.mark.xfail
    def test_JournalSequenceError_skip_stage(self, journal):
        # Currently allows skipping a stage and does not enforce ALL previous
        # stages to have completed.
        with pytest.raises(restart.JournalSequenceError, match="Cannot start stage"):
            journal.start("pre")
            assert journal.current == "pre"
            journal.completed("pre")

            journal.start("post")

    def test_start_idempotent(self, journal):
        # test that start() can be called multiple time (#278)
        journal.start("pre")
        journal.start("pre")
        assert journal.current == "pre"

    def test_incomplete_known_stage(self, journal):
        journal.incomplete = "main"
        assert journal.incomplete == "main"

    def test_incomplete_unknown_stage_ValueError(self, journal):
        with pytest.raises(ValueError, match="Can only assign a registered stage from"):
            journal.incomplete = "BEGIN !"

    def test_clear(self, journal):
        journal.start("pre")
        journal.completed("pre")
        journal.start("main")
        # manually setting incomplete
        journal.incomplete = journal.current

        assert journal.current == "main"
        assert journal.incomplete == journal.current

        journal.clear()
        assert journal.current is None
        assert journal.incomplete is None

    def test_history(self, journal):
        journal.start("pre")
        journal.completed("pre")
        journal.start("main")
        journal.completed("main")
        journal.start("post")

        # completed stages
        assert journal.history == ["pre", "main"]

    def test_history_del(self, journal):
        journal.start("pre")
        journal.completed("pre")
        journal.start("main")
        journal.completed("main")
        assert journal.history

        del journal.history
        assert journal.history == []

    def test_has_completed(self, journal):
        journal.start("pre")
        journal.completed("pre")

        assert journal.has_completed("pre")
        assert not journal.has_completed("main")

    def test_has_not_completed(self, journal):
        journal.start("pre")
        journal.completed("pre")
        journal.start("main")
        # simulate crash/restart
        del journal.current

        assert journal.has_not_completed("main")
        assert not journal.has_not_completed("pre")


# need a real class so that it can be pickled later
class JournalledMemory(restart.Journalled):
    # divide is a dummy protocol
    protocols = ["divide", "multiply"]

    def __init__(self):
        self.memory = 1
        super().__init__()

    def multiply(self, x):
        self.journal.start("multiply")
        self.memory *= x
        self.journal.completed("multiply")


@pytest.fixture
def journalled():
    return JournalledMemory()


class TestJournalled:
    @staticmethod
    def divide(m, x):
        return m.memory / x

    def test_get_protocol_of_class(self, journalled):
        f = journalled.get_protocol("multiply")
        f(10)
        assert journalled.memory == 10
        assert journalled.journal.has_completed("multiply")

    def test_get_protocol_dummy(self, journalled):
        dummy_protocol = journalled.get_protocol("divide")
        result = dummy_protocol(self.divide, journalled, 10)

        assert result == 1 / 10
        assert journalled.journal.has_completed("divide")

    def test_get_protocol_dummy_incomplete(self, journalled):
        dummy_protocol = journalled.get_protocol("divide")
        with pytest.raises(ZeroDivisionError):
            result = dummy_protocol(self.divide, journalled, 0)
        assert not journalled.journal.has_completed("divide")

    def test_save_load(self, tmp_path):
        # instantiate a class that can be pickled (without pytest magic)
        journalled = JournalledMemory()
        f = journalled.get_protocol("multiply")
        f(10)
        assert journalled.memory == 10

        pickle = tmp_path / "memory.pkl"
        journalled.save(pickle)

        assert pickle.exists()

        # change instance
        f(99)
        assert journalled.memory == 10 * 99

        # reload previous state
        journalled.load(pickle)
        assert journalled.memory == 10
