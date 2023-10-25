"""
Requirements:
* pytest
* pytest-html

To run:
pytest -s --html=report.html
pytest --cov=acpype --cov-report=term-missing:skip-covered
"""

import os
from shutil import rmtree

import pytest

from acpype import __version__ as version


def pytest_report_header(config):
    return f">>>\tVersion: {version}\n"


@pytest.fixture
def janitor():
    os.chdir(os.path.dirname(__file__))
    to_delete = []
    yield to_delete
    for item in to_delete:
        if os.path.exists(item):
            rmtree(item)
