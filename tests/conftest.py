from acpype_lib.acpype import version


def pytest_report_header(config):
    return f">>>\tVersion: {version}\n"
