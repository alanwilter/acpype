import os
from subprocess import run, STDOUT, PIPE


def run_git():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    return run("git describe --tags --always", shell=True, stderr=STDOUT, stdout=PIPE)


out = run_git()

if out.returncode == 0:

    version = out.stdout.decode()[0:10]
    version = out.stdout.decode().rsplit("-", 1)[0]
else:
    try:
        from importlib.metadata import version as ver

        version = ver("acpype")
    except Exception:
        import pkg_resources

        version = str(pkg_resources.get_distribution("acpype").version)
