from subprocess import Popen, STDOUT, PIPE

try:
    try:
        from importlib.metadata import version

        version = str(version("acpype"))
    except Exception:
        import pkg_resources

        version = str(pkg_resources.get_distribution("acpype").version)
except Exception:
    version = (
        Popen("git describe --tags --always", shell=True, stderr=STDOUT, stdout=PIPE)
        .communicate()[0][:-1]
        .decode()[0:10]
    )
