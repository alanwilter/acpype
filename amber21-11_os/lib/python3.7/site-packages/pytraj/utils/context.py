import os
from contextlib import contextmanager
import tempfile
from shutil import rmtree

from ..externals.wurlitzer import pipes
capture_stdout = pipes


@contextmanager
def tempfolder():
    my_temp = tempfile.mkdtemp()
    cwd = os.getcwd()
    os.chdir(my_temp)
    yield
    os.chdir(cwd)
    rmtree(my_temp)
