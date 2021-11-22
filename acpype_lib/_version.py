import os
from subprocess import run, STDOUT, PIPE


def run_git():
    cur_dir = os.getcwd()
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    ver = run("git describe --tags --always", shell=True, stderr=STDOUT, stdout=PIPE)
    os.chdir(cur_dir)
    return ver


out = run_git()
breakpoint()

if out.returncode == 0:

    version = out.stdout.decode().rsplit("-", 1)[0].replace("-", ".").strip()
else:
    try:
        from importlib.metadata import version as ver

        version = ver("acpype")
    except Exception:
        import pkg_resources

        version = str(pkg_resources.get_distribution("acpype").version)
