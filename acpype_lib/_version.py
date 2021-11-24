import os
import re
from subprocess import run, STDOUT, PIPE
from importlib.metadata import version as ver
from datetime import datetime

today = datetime.today().strftime("%Y.%m.%d")


def parse_ver(ver):
    m = re.match(r"(\d{4}.\d{2}.\d{2}(-\d+)?)(\-.*)?", ver)
    if m:
        version = m.group(1).replace("-", ".")
        status = 0
    else:
        version = "0"
        status = 1
    return version, status


def run_git():
    cur_dir = os.getcwd()
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    ver = run("gits describe --tags --always", shell=True, stderr=STDOUT, stdout=PIPE)
    os.chdir(cur_dir)
    return parse_ver(ver.stdout.decode())


version, status = run_git()

if status != 0:
    try:
        version, status = parse_ver(ver("acpype"))
    except Exception:
        version = today
        status = 0
if status != 0:
    version = today
