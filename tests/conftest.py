from acpype_lib import version


def pytest_report_header(config):
    return f">>>\tVersion: {version}\n"
