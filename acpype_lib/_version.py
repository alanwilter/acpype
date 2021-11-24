import os
import re
from subprocess import run, STDOUT, PIPE


def run_git():
    cur_dir = os.getcwd()
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    ver = run("git describe --tags --always", shell=True, stderr=STDOUT, stdout=PIPE)
    os.chdir(cur_dir)
    m = re.match(r"(\d{4}.\d{2}.\d{2}(-\d+)?)(\-.*)", ver.stdout.decode())
    if m:
        ver = m.group(1).replace("-", ".")
        status = 0
    else:
        ver = "0"
        status = 1
    return ver, status


version, status = run_git()

if status != 0:
    try:
        from importlib.metadata import version as ver

        version = ver("acpype")
    except Exception:
        version = "0"
