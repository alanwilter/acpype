from acpype_lib import __version__ as version


def pytest_report_header(config):
    return f">>>\tVersion: {version}\n"
