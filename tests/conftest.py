"""
Shared pytest fixtures.

To run:
pytest
pytest --cov=acpype --cov-report=term-missing:skip-covered
"""

import os
import shutil
import tempfile
from pathlib import Path
from shutil import rmtree

import pytest

from acpype import __version__ as version

TESTS_DIR = Path(__file__).parent


def pytest_report_header(config):
    return f">>>\tVersion: {version}\n"


def _data_files():
    """Yield the input files the tests read, i.e. everything but the test modules."""
    for path in sorted(TESTS_DIR.iterdir()):
        if path.is_file() and path.suffix != ".py":
            yield path


@pytest.fixture(scope="session")
def data_dir(tmp_path_factory):
    """One pristine copy of the test inputs, shared by the whole session.

    Tests are linked against this copy rather than the files in the repository, so a
    test that writes to an input name cannot modify the checked-in data.
    """
    dest = tmp_path_factory.mktemp("acpype-data")
    for src in _data_files():
        shutil.copy2(src, dest / src.name)
    return dest


@pytest.fixture
def janitor(tmp_path, monkeypatch, data_dir):
    """Run a test in its own directory, with the test inputs available by name.

    ACPYPE and antechamber both write their scratch files into the working directory.
    Running every test in the shared ``tests/`` directory meant those files collided
    between tests, left artefacts in the repository, and made two concurrent pytest
    runs fail unpredictably. Each test now gets an empty directory of its own.

    Yields:
        list: paths to remove after the test; kept for the tests that append to it,
        though anything written inside the temporary directory is discarded anyway.
    """
    for src in data_dir.iterdir():
        (tmp_path / src.name).symlink_to(src)
    monkeypatch.chdir(tmp_path)
    # acpype.utils.parmMerge caches merged parameter files in the temp directory.
    # Pointing that at tmp_path keeps a warm cache on the developer's machine from
    # skipping the merge, which otherwise makes coverage vary by several percent.
    monkeypatch.setattr(tempfile, "tempdir", str(tmp_path))

    to_delete = []
    yield to_delete
    for item in to_delete:
        if os.path.exists(item):
            rmtree(item)
