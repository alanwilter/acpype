"""
Requirements:
* pytest
* pytest-html

To run:
pytest -s --html=report.html
pytest --cov=acpype --cov-report=term-missing:skip-covered
"""

from acpype import __version__ as version


def pytest_report_header(config):
    return f">>>\tVersion: {version}\n"
