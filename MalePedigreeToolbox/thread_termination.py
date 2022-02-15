import os
import sys
from pathlib import Path
import logging
import uuid
import tempfile


# in case runnign from the exe --> might return different paths depending on the way of executing
if getattr(sys, 'frozen', False):
    static_path = Path(os.path.dirname(sys.executable))
elif __file__:
    static_path = Path(os.path.dirname(__file__))
else:
    raise SystemExit("Can not find application path.")


# ensure unique semaphore for different instances of the program
_TEMP_FILE = tempfile.NamedTemporaryFile(delete=False)
_SEMAPHORE_FILE = Path(_TEMP_FILE.name)
_DISABLED = False

LOG = logging.getLogger("mpt")


def clear_semaphore():
    """Make sure that the semaphore is cleared and present. The file should have been deleted, but you never know"""
    open(_SEMAPHORE_FILE, "w").close()


def disable_thread_termination():
    """Disable thread termination, usefull when not working with threads and wanting to prevent the semaphore opening
    overhead"""
    global _DISABLED
    _DISABLED = True


def enable_thread_termination():
    """enable thread termination, used when starting the gui from the command line"""
    global _DISABLED
    _DISABLED = False


def request_exit():
    with open(_SEMAPHORE_FILE, "w") as f:
        f.write("1")


class ThreadTerminable:
    """Simply checks if a semaphore file has a value in it. If that is the case terminate execution. This allows to
    signal a thread to terminate pre-maturely

    Use this as little as possible since it relies on file reading which can slow down an application slightly
    """

    def __init__(self, function):
        # ensure this file is present since it is integral to working of the wrapper
        # clearing the file can not be done here since it can not be guaranteed that a written value will than always
        # cause an abort
        self._function = function
        if not _SEMAPHORE_FILE.exists():
            open(_SEMAPHORE_FILE, "a").close()

    def __call__(self, *args, **kwargs):
        # make sure that semaphore checkign can be skipped
        if _DISABLED:
            return self._function(*args, **kwargs)
        if not _SEMAPHORE_FILE.exists():
            raise SystemExit("'semaphore' is missing. Can not ensure pre-mature thread termintation")
        do_exit = False
        with open(_SEMAPHORE_FILE) as f:
            if f.read() == "1":
                do_exit = True
        if do_exit:
            LOG.info("Process has been succesfully termintated.")
            _TEMP_FILE.close()
            os.unlink(_TEMP_FILE.name)
            raise SystemExit("User requested abort")
        else:
            return self._function(*args, **kwargs)

    def __del__(self):
        # remove the file on application exit
        try:
            if _SEMAPHORE_FILE.exists():
                _TEMP_FILE.close()
                os.unlink(_TEMP_FILE.name)
        except FileNotFoundError:
            pass
