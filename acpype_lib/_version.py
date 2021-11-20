from subprocess import Popen, STDOUT, PIPE
from pbr.version import VersionInfo

try:
    version = VersionInfo("acpype").version_string()
except Exception:
    version = (
        Popen("git describe --tags --always", shell=True, stderr=STDOUT, stdout=PIPE)
        .communicate()[0][:-1]
        .decode()[0:10])
