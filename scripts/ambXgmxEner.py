"""
Script to get energies from AMBER and GROMACS
for a given system and compare them
"""

import re
import sys

from acpype.utils import _getoutput


def error(v1, v2):
    """percentage relative error"""
    if v1 == v2:
        return 0
    return abs(v1 - v2) / max(abs(v1), abs(v2)) * 100


ff = 4.1840

vv = ["ANGLE", "BOND", "DIHED", "EPTOT", "VDW14", "VDW", "QQ14", "QQ"]
norm_gmx = {"POTENTIAL": "EPTOT", "LJ-14": "VDW14", "LJ_SR": "VDW", "COULOMB-14": "QQ14", "COULOMB_SR": "QQ"}
norm_amb = {"1-4_NB": "VDW14", "VDWAALS": "VDW", "1-4_EEL": "QQ14", "EELEC": "QQ"}

(e_amb, e_gmx) = sys.argv[1:]

cmd_amb = f"cat {e_amb}"
cmd_gmx = f"echo 1 2 3 4 5 6 7 8 9 | gmx energy -f {e_gmx}.edr"

amb_out = {
    y[0].upper(): float(y[1]) * ff
    for y in [
        x.split("=") for x in re.sub(r"\s+=\s+", "=", _getoutput(cmd_amb)).strip().replace("1-4 ", "1-4_").split()
    ]
}

for k in list(amb_out.keys())[:]:
    v = norm_amb.get(k)
    if v:
        amb_out[v] = amb_out[k]

gmx_out = {
    y[0].upper(): float(y[1])
    for y in [
        x.split()[:2]
        for x in _getoutput(cmd_gmx)
        .split("-" * 3)[-1][2:]
        .replace(" Dih.", "_Dih")
        .replace(" (SR)", "_SR")
        .splitlines()
    ]
}
gmx_out["DIHED"] = gmx_out.get("PROPER_DIH", 0) + gmx_out.get("IMPROPER_DIH", 0) + gmx_out.get("RYCKAERT-BELL.", 0)
for k in list(gmx_out.keys())[:]:
    v = norm_gmx.get(k)
    if v:
        gmx_out[v] = gmx_out[k]

for ii in vv:
    print(f"{ii:<10}: {error(gmx_out[ii],amb_out[ii]):-2.5f}%")

# DIHED = PROPER_DIH + IMPROPER_DIH + RYCKAERT-BELL.
# VDW14 : LJ-14, 1-4_NB
# VDW   : LJ_SR, VDWAALS
# QQ14  : COULOMB-14, 1-4_EEL
# QQ    : COULOMB_SR, EELEC
