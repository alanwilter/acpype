from subprocess import Popen, STDOUT, PIPE
from pbr.version import VersionInfo

try:
    __version__ = Popen("git describe --tags --always", shell=True, stderr=STDOUT, stdout=PIPE).communicate()[0][:-1].decode()[0:10]
except Exception:
    __version__ = VersionInfo("acpype").version_string()
