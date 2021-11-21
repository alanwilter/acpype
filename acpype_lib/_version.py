from subprocess import Popen, STDOUT, PIPE
import re


version = (
    Popen("git describe --tags --always", shell=True, stderr=STDOUT, stdout=PIPE).communicate()[0][:-1].decode()[0:10]
)

if not re.match("^\d{4}\.(0[1-9]|1[012])\.(0[1-9]|[12][0-9]|3[01])$", version):
    try:
        from importlib.metadata import version as ver

        version = ver("acpype")
    except Exception:
        import pkg_resources

        version = str(pkg_resources.get_distribution("acpype").version)
