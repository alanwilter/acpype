#!/usr/bin/env python3
"""Compare ACPYPE topologies and energies against GROMACS/AMBER99SB references.

This is the comparative test described in ``NOTE.txt``. It builds the 22 tripeptides
(one per usual amino acid, plus the three histidine protonation states), converts each
one with both ``pdb2gmx -ff amber99sb`` and ACPYPE, and reports where the atom types,
bonded parameters and single-point potential energies disagree.

It is a development harness, not part of the test suite: it needs external tools, takes
minutes to run, and its output is meant to be read rather than asserted on.

Requirements:
    * GROMACS on ``PATH`` (or ``$GMX``), providing ``gmx`` and its ``amber99sb.ff``
      data. ``brew install gromacs`` or ``conda install -c conda-forge gromacs``.
    * AmberTools for ``tleap``. The copy bundled with ACPYPE is used automatically.
    * ``sander`` and ``ambpdb``, to energy minimise each tripeptide. These are *not*
      part of the AmberTools bundled with ACPYPE, so a full install is needed:
      ``conda install -c conda-forge ambertools``. The script runs without them, using
      the geometry ``tleap savepdb`` writes, but that geometry is strained enough that
      antechamber mis-perceives aromatic rings -- PHE's CD2 comes out as type ``DU``,
      which Gasteiger has no parameter for -- and the energies drift far from the
      figures quoted in ``NOTE.txt``. Treat a run without them as indicative only.
    * Optionally ``pymol``, an alternative way to build the tripeptides (``--pymol``).

Usage:
    uv run python scripts/check_acpype.py [--ff amber|gaff] [--charge gas|bcc]
                                          [--dir DIR] [--pymol] [--keep]
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from glob import glob

from acpype.topol import ACTopol
from acpype.utils import _getoutput


def _parse_args():
    """Read the command line, keeping the historical defaults."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--ff",
        default="amber",
        choices=("amber", "amber2", "gaff", "gaff2"),
        help="ACPYPE atom types to compare (only 'amber'/'amber2' are like-for-like with the AMBER reference)",
    )
    parser.add_argument("--charge", default="gas", choices=("gas", "bcc"), help="ACPYPE charge method")
    parser.add_argument(
        "--protein_ff", default="ff14SB", choices=("ff14SB", "ff99SB"), help="protein force field for -a amber"
    )
    parser.add_argument(
        "--gmxff",
        default=None,
        choices=("amber14sb", "amber99sb"),
        help="GROMACS reference force field (default: the one matching --protein_ff)",
    )
    parser.add_argument("--dir", default=None, help="working directory (default: a fresh temporary one)")
    parser.add_argument("--pymol", action="store_true", help="build the tripeptides with pymol instead of tleap")
    parser.add_argument("--keep", action="store_true", help="do not delete the working directory afterwards")
    return parser.parse_args()


def _which(name, *extraDirs):
    """Find an executable next to the bundled AmberTools first, then on PATH."""
    for directory in extraDirs:
        candidate = os.path.join(directory, name)
        if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            return candidate
    return shutil.which(name)


def _gmxDataDir(gmxExe):
    """Locate the GROMACS ``share`` tree holding ``gromacs/top/amber99sb.ff``."""
    if os.environ.get("GMXDATA"):
        return os.environ["GMXDATA"]
    # gmx is usually a symlink into a versioned prefix, e.g. Homebrew's Cellar.
    prefix = os.path.dirname(os.path.dirname(os.path.realpath(gmxExe)))
    return os.path.join(prefix, "share")


args = _parse_args()

usePymol = args.pymol
ffType = args.ff
cType = args.charge
proteinFF = args.protein_ff
gmxFF = args.gmxff or {"ff14SB": "amber14sb", "ff99SB": "amber99sb"}[proteinFF]
debug = True

water = " -water none"

print(f"usePymol: {usePymol}, ffType: {ffType}, cType: {cType}, proteinFF: {proteinFF}, reference: {gmxFF}")

tmpDir = args.dir or tempfile.mkdtemp(prefix="acpype-check-")

delList = ["topol.top", "posre.itp"]

# create Dummy PDB, only so ACTopol resolves the AmberTools binaries for us
tpdb = os.path.join(tempfile.gettempdir(), "acpype-check-dummy.pdb")
dummyLine = "ATOM      1  N   ALA A   1      -1.188  -0.094   0.463  1.00  0.00           N\n"
open(tpdb, "w").writelines(dummyLine)
tempObj = ACTopol(tpdb, chargeVal=0, verbose=False)

# Binaries for ambertools
acExe = tempObj.acExe
tleapExe = tempObj.tleapExe
amberBin = os.path.dirname(acExe)
# sander and ambpdb are not part of the AmberTools bundled with acpype, and are only
# needed to minimise the tripeptides. Without them tleap's savepdb geometry is used.
sanderExe = _which("sander", amberBin)
ambpdbExe = _which("ambpdb", amberBin)

exePymol = _which("pymol")

# Binaries for gromacs
gmxExe = os.environ.get("GMX") or _which("gmx") or _which("gmx_mpi")
if not gmxExe:
    sys.exit("ERROR: no 'gmx' found. Install GROMACS or point $GMX at it.")
gmxTopDir = _gmxDataDir(gmxExe)
if not os.path.isdir(os.path.join(gmxTopDir, "gromacs", "top", f"{gmxFF}.ff")):
    sys.exit(f"ERROR: no {gmxFF}.ff under {gmxTopDir}/gromacs/top. Set $GMXDATA to the GROMACS share directory.")
pdb2gmx = f"{gmxExe} pdb2gmx"
grompp = f"{gmxExe} grompp"
mdrun = f"{gmxExe} mdrun"
editconf = f"{gmxExe} editconf"
g_energy = f"{gmxExe} energy"
gmxdump = f"{gmxExe} dump"

print(
    f"gmx: {gmxExe}\ntleap: {tleapExe}\nminimise: {'yes' if sanderExe and ambpdbExe else 'no (sander/ambpdb absent)'}"
)

amberff = "oldff/leaprc.ff99SB"

genPdbTemplate = """
source %(amberff)s
mol = sequence {N%(res1)s %(res)s C%(res1)s}
savepdb mol %(aai3)s.pdb
saveamberparm mol prmtop inpcrd
quit
"""

# Single point energy: one force evaluation. GROMACS dropped the group cutoff scheme
# in 2020, and Verlet supports neither infinite cutoffs nor pbc=no, so the molecule now
# goes in a large box (editconf -d 3.0) with cutoffs well beyond its own size. Both
# sides of the comparison use these same settings, so the difference is unaffected.
# This used to carry its shell heredoc wrapper ("cat << EOF >| SPE.mdp" ... "EOF")
# into the file itself, which modern grompp rejects as three malformed parameters.
spe_mdp = """integrator               = md
nsteps                   = 0
dt                       = 0.001
constraints              = none
cutoff-scheme            = Verlet
nstlist                  = 1
verlet-buffer-tolerance  = -1
rlist                    = 2.0
pbc                      = xyz
coulombtype              = Cut-off
rcoulomb                 = 2.0
vdwtype                  = Cut-off
rvdw                     = 2.0
Tcoupl                   = no
Pcoupl                   = no
gen_vel                  = no
nstcomm                  = 1
comm_mode                = Linear
nstlog                   = 1
nstenergy                = 1
nstcalcenergy            = 1
nstxout                  = 1
continuation             = yes
"""

aa_dict = {
    "A": "ala",
    "C": "cys",
    "D": "asp",
    "E": "glu",
    "F": "phe",
    "G": "gly",
    "H": "hie",
    "J": "hip",
    "O": "hid",
    "I": "ile",
    "K": "lys",
    "L": "leu",
    "M": "met",
    "N": "asn",
    "P": "pro",
    "Q": "gln",
    "R": "arg",
    "S": "ser",
    "T": "thr",
    "V": "val",
    "W": "trp",
    "Y": "tyr",
}

hisDictConv = {"HHH": ["HIE", "HISE"], "JJJ": ["HIP", "HISH"], "OOO": ["HID", "HISD"]}

pymolScript = (
    "aa_dict = %s" % str(aa_dict)
    + """
def build3(object_name, sequence, first_residue = "1"):
    if len(sequence):
        codeNT, code, codeCT = sequence
        cmd.fragment('nt_'+aa_dict[codeNT],object_name)
        cmd.alter(object_name,'resi="%s"'%first_residue)
        cmd.edit(object_name+" and name C")
        for code in sequence[1:-1]:
            editor.attach_amino_acid("pk1",aa_dict[code])
        editor.attach_amino_acid("pk1",'ct_'+aa_dict[codeCT])
        cmd.edit()

seqi=aa_dict.keys()

count=0

cmd.delete('all')

for aai in seqi:
    aai3=3*aai
    build3(aai3,aai3)
    cmd.alter(aai3,"chain='%s'"%aai3)
    cmd.translate([15*count,0,0],aai3)
    cmd.sculpt_activate(aai3)
    for n in range(100):
        cmd.sculpt_iterate(aai3)
    count=count+1
    aafile = aai3+'.pdb'
    cmd.set('pdb_conect_all')
    cmd.save(aafile,'all')
    cmd.delete('all')
"""
)


def pdebug(msg):
    if debug:
        print("DEBUG: %s" % msg)


def perror(msg):
    if debug:
        print("ERROR: %s" % msg)


def createAmberPdb(fpdb, opt):
    """Pdb -> apdb (amber)
    opt = 0 : pymol pdb to apdb
    opt = 1 : ogpdb to apdb
    """
    fLines = open(fpdb).readlines()
    fName = "a" + fpdb
    famb = open(fName, "w")
    for line in fLines:
        if "ATOM  " in line:
            if opt == 0:
                if "AAA" not in fpdb:
                    line = line.replace(" HB3 ", " HB1 ")
                line = line.replace("CD1 ILE", "CD  ILE")
                line = line.replace(" HG3 ", " HG1 ")
                line = line.replace(" HA3 GLY", " HA1 GLY")
                line = line.replace("HG13 ILE", "HG11 ILE")
                line = line.replace("HG13 ILE", "HG11 ILE")
                if line[22:26] == "   3":
                    line = line.replace("  O  ", "  OC1").replace("  OXT", "  OC2")
            elif opt == 1:
                if line[22:26] == "   3":
                    line = line.replace(" O1 ", " OC1").replace(" O2 ", " OC2")
            if line[22:26] == "   1":
                res = line[17:20]
                line = line.replace("%s " % res, "N%s" % res)
            if line[22:26] == "   3":
                res = line[17:20]
                line = line.replace("%s " % res, "C%s" % res)
        famb.write(line)  # pay attention here!
    famb.close()
    return fName


def createAmberPdb3(fpdb):
    """Add N and C XXX --> NXXX, CXXX"""
    fLines = open(fpdb).readlines()
    fName = fpdb.replace("_", "a")
    famb = open(fName, "w")
    for line in fLines:
        if "ATOM  " in line:
            if line[22:26] == "   1":
                res = line[17:20]
                line = line.replace("%s " % res, "N%s" % res)
            if line[22:26] == "   3":
                res = line[17:20]
                line = line.replace("%s " % res, "C%s" % res)
        famb.write(line)  # pay attention here!
    famb.close()
    return fName


def createAmberPdb2(fpdb, opt):
    """Use formatConverter
    Add N and C XXX --> NXXX, CXXX and OC1 & OC2
    """
    projectName = "testImport"
    if os.path.exists(projectName):
        shutil.rmtree(projectName)
    _fpdb = "_%s" % fpdb

    fLines = open(_fpdb).readlines()
    fName = "a" + fpdb
    famb = open(fName, "w")
    for line in fLines:
        if "ATOM  " in line:
            line = line.replace("CYS", "CYN")
            line = line.replace("LYS", "LYP")
            if opt == 0:
                if line[22:26] == "   3":
                    line = line.replace("  O   ", "  OC1 ").replace("  OXT ", "  OC2 ")
            elif opt == 1:
                if line[22:26] == "   3":
                    line = line.replace(" O1 ", " OC1").replace(" O2 ", " OC2")
            if line[22:26] == "   1":
                res = line[17:20]
                line = line.replace("%s " % res, "N%s" % res)
            if line[22:26] == "   3":
                res = line[17:20]
                line = line.replace("%s " % res, "C%s" % res)
        famb.write(line)  # pay attention here!
    famb.close()
    return fName


def createOldPdb2(fpdb):
    """Using my own dict for e.g. 2HB = HB2 -> HB1"""
    defHB1 = [" HB2", "2HB ", " HB1"]
    defHB2 = [" HB3", "3HB ", " HB2"]
    defHG1 = [" HG2", "2HG ", " HG1"]
    defHG2 = [" HG3", "3HG ", " HG2"]
    defHD1 = [" HD2", "2HD ", " HD1"]
    defHD2 = [" HD3", "3HD ", " HD2"]
    def1HG2 = ["HG21", "1HG2"]
    def2HG2 = ["HG22", "2HG2"]
    def3HG2 = ["HG23", "3HG2"]

    dictPdb2GmxAtomNames = {
        "ALA": (),
        "VAL": (),
        "CYS": (defHB1, defHB2),
        "ASP": (defHB1, defHB2),
        "GLU": (defHB1, defHB2, defHG1, defHG2),
        "PHE": (defHB1, defHB2),
        "GLY": ([" HA2 ", " HA  ", " HA1 "], [" HA3", "3HA ", " HA2"]),
        "HIE": (defHB1, defHB2),
        "ILE": (
            ["CD1", "CD "],
            ["HG12", "2HG1", "1HG1"],
            ["HG13", "3HG1", "2HG1"],
            def1HG2,
            def2HG2,
            def3HG2,
            ["HD11", "1HD1", " HD1"],
            ["HD12", "2HD1", " HD2"],
            ["HD13", "3HD1", " HD3"],
        ),
        "HIP": (defHB1, defHB2),
        "HID": (defHB1, defHB2),
        "LYS": (
            defHB1,
            defHB2,
            defHG1,
            defHG2,
            defHD1,
            defHD2,
            [" HE2", "2HE ", " HE1"],
            [" HE3", "3HE ", " HE2"],
        ),
        "LEU": (defHB1, defHB2),
        "MET": (defHB1, defHB2, defHG1, defHG2),
        "ASN": (defHB1, defHB2),
        "PRO": (
            defHB1,
            defHB2,
            defHG1,
            defHG2,
            defHD1,
            defHD2,
            [" H3 ", "3H  ", " H1 "],
        ),
        "GLN": (defHB1, defHB2, defHG1, defHG2, ["HE21", "1HE2"], ["HE22", "2HE2"]),
        "ARG": (defHB1, defHB2, defHG1, defHG2, defHD1, defHD2),
        "SER": (defHB1, defHB2),
        "THR": (def1HG2, def2HG2, def3HG2),
        "TRP": (defHB1, defHB2),
        "TYR": (defHB1, defHB2),
    }
    fName = "_" + fpdb
    fLines = open(fpdb).readlines()
    aLines = []
    for line in fLines:
        if "ATOM  " in line:
            aLines.append(line)
    nRes = int(aLines[-1][22:26])
    nLines = []
    for line in aLines:
        res = line[17:20]
        nResCur = int(line[22:26])
        for item in dictPdb2GmxAtomNames[res]:
            atAim = item[-1]
            for at in item[:-1]:
                # print(at, atAim)
                # print(line)
                line = line.replace("%s" % at, "%s" % atAim)
                # print(line)
        if nResCur == nRes:
            line = line.replace("O   %s" % res, "OC1 %s" % res)
            line = line.replace("OXT %s" % res, "OC2 %s" % res)
        nLines.append(line)
    open(fName, "w").writelines(nLines)
    return fName


def parseTopFile(lista):
    """Lista = top file in list format
    return a (dict,dict) with fields
    """
    defDict = {}
    parDict = {}
    for line in lista:
        line = line.split(";")[0]
        line = line.strip()
        if line.startswith("#def"):
            vals = line.split()
            defDict[vals[1]] = map(eval, vals[2:])
        elif line and not line.startswith("#"):
            if line.startswith("[ "):
                flag = line.split()[1]
                if flag not in parDict:
                    parDict[flag] = []
            else:
                u = []
                t = line.split()
                for i in t:
                    try:
                        v = eval(i)
                    except Exception:
                        v = i
                    u.append(v)
                parDict[flag].append(u)
    return parDict, defDict


def nbDict(lista):
    tdict = {}
    for line in lista:
        line = line.split(";")[0]
        line = line.strip()
        if line and not line.startswith("#") and not line.startswith("["):
            name, type_ = line.split()[0:2]
            tdict[name] = type_
    return tdict


def roundAllFloats(lista, ll):
    """Round to 3 decimals"""
    nlista = []
    for ii in lista:
        tt = ii[: ll + 1]
        for jj in ii[ll + 1 :]:
            if jj > 100.0:
                jj = round(jj, -1)
            nn = round(jj, 3)
            tt.append(nn)
        nlista.append(tt)
    return nlista


def checkTopAcpype(res):
    """Compare acpype gmx itp against amber pdb2gmx results"""
    os.chdir(tmpDir)

    def addParam(ll, item):
        """Ll : index for bonded types
        item : set of atomtypes
        """
        dict_ = {2: "bondtypes", 3: "angletypes", 4: "dihedraltypes"}
        dType = {}  # dict_
        for type_ in ffBon[0][dict_[ll]]:
            i = type_[: ll + 1]
            j = type_[ll + 1 :]
            dType[str(i)] = j  # dict_ {[atomtypes,funct] : parameters}
        entries = []
        lNum = item[ll]  # funct
        ent = [ffDictAtom[x] for x in item[:ll]]  # convert atomtypes ids to atomtypes names
        rent = ent[:]
        rent.reverse()
        entries.append([*ent, lNum])
        entries.append([*rent, lNum])
        if ll == 4:
            if len(item) == 6:
                par = ffBon[1][item[5]]
                return par
            tent = ent[:]
            rtent = rent[:]
            if lNum in [3, 9]:  # dih proper
                tent[0] = "X"
                entries.append([*tent, lNum])
                tent.reverse()
                entries.append([*tent, lNum])
                tent[0] = "X"
                entries.append([*tent, lNum])
                tent.reverse()
                entries.append([*tent, lNum])
                rtent[0] = "X"
                entries.append([*rtent, lNum])
                rtent.reverse()
                entries.append([*rtent, lNum])
            if lNum in [1, 4]:  # dih improp
                tent[0] = "X"
                entries.append([*tent, lNum])
                tent[1] = "X"
                entries.append([*tent, lNum])
                rtent[0] = "X"
                entries.append([*rtent, lNum])
                rtent[1] = "X"
                entries.append([*rtent, lNum])
        found = False
        for e in entries:
            try:
                par = dType[str(e)]
                found = True
                break
            except Exception:
                pass
        if not found:
            print(f"{dict_[ll]} {ent} {item} not found in {ffType} Bon")
        # print(item, e, par)
        return par

    compareParameters = True

    agRes = parseTopFile(open("ag%s.top" % (res)).readlines())
    acRes = parseTopFile(open(f"ag{res}.acpype/ag{res}_GMX.itp").readlines())
    ffBon = aBon
    ffgRes = agRes

    atError = False
    print("    ==> Comparing atomtypes AC x AMBER")

    ffDictAtom = {}  # dict link res atom numbers to amber atom types
    for item in ffgRes[0]["atoms"]:
        i, j = item[:2]  # id, at
        ambat = j
        ffDictAtom[i] = ambat
        # compare atom types AC x Amber
        acat = acRes[0]["atoms"][i - 1][1]
        if ambat != acat:
            print("    %i %s AC: %s   x   AMB: %s" % (i, item[4], acat, ambat))
            atError = True

    if not atError:
        print("        atomtypes OK")

    acDictAtom = {}
    for item in acRes[0]["atoms"]:
        i, j = item[:2]
        acDictAtom[i] = j  # dict for atom id -> atom type from acpype itp file

    flags = [
        ("pairs", 2),
        ("bonds", 2),
        ("angles", 3),
        ("dihedrals", 4),
    ]  # , ['dihedraltypes', 'angletypes', 'bondtypes']

    for flag, ll in flags:
        print("    ==> Comparing %s" % flag)
        if flag != flags[0][0] and compareParameters:  # not 'pairs'
            agres = []
            tAgRes = ffgRes[0][flag]  # e.g. dic 'bonds' from gmx top file
            for item in tAgRes:
                if flag == flags[1][0] or flag == flags[2][0]:  # 'bonds'
                    par = addParam(ll, item)
                elif flag == flags[3][0]:  # 'dihedrals', e.g. item = [2, 1, 5, 6, 9]
                    par = addParam(ll, item)
                    if len(item) == 6:
                        item.pop()
                agres.append(item + par)
            if compareParameters:
                acres = acRes[0][flag]
        else:
            if flag == flags[3][0]:
                agres = [x[: ll + 1] for x in ffgRes[0][flag]]
            else:
                agres = ffgRes[0][flag]
            if compareParameters:
                acres = acRes[0][flag]
            else:
                acres = [x[: ll + 1] for x in acRes[0][flag]]

        act = acres[:]
        agt = agres[:]

        if compareParameters:
            if flag != flags[0][0]:
                # round parameters
                act = roundAllFloats(act, ll)
                agt = roundAllFloats(agt, ll)
                act2 = act[:]
                agt2 = agt[:]

        if not compareParameters or flag == flags[0][0]:
            act2 = act[:]
            agt2 = agt[:]

        for ac in act:
            if str(ac) in str(agt):
                act2.remove(ac)
                agt2.remove(ac)
            else:
                t = ac[:-1]
                t.reverse()
                acr = [*t, ac[-1]]
                if str(acr) in str(agt):
                    act2.remove(ac)
                    agt2.remove(acr)

        act3 = act2[:]
        agt3 = agt2[:]

        if act2 and agt2:
            # specially for dih since it may need to resort indexes
            agl = {}
            for ag in agt3:
                t = ag[:-1]
                t.sort()
                ags = [*t, ag[-1]]
                agl[str(ags)] = ag
            for ac in act2:
                t = ac[:-1]
                t.sort()
                acs = [*t, ac[-1]]
                if str(acs) in str(agl.keys()):
                    act3.remove(ac)
                    agt3.remove(agl[str(acs)])

        act4 = []
        agt4 = []

        for ac in act3:
            if ac[:5] not in tAgRes:
                if flag == flags[3][0]:
                    if ac[6]:
                        act4.append(ac)
                else:
                    act4.append(ac)

        tAcRes = [x[:5] for x in acres]
        for ac in agt3:
            if ag[:5] not in tAcRes:
                act4.append(ac)

        if flag == flags[3][0]:
            for ac in act4:
                if not ac[6]:
                    act4.remove(ac)
            for ag in agt4:
                if not ag[6]:
                    agt4.remove(ac)

        if act4:
            act4.sort()
            print("    ac: ", act4)
        if agt4:
            agt4.sort()
            print("    %sg: " % ff, agt4)

        if not act4 and not agt4:
            print("        %s OK" % flag)


def parseOut(out):
    lines = out.splitlines()
    # for line in lines:
    count = 0
    while count < len(lines):
        line = lines[count]
        if "WARNING" in line.upper():
            if (
                "will be determined based" in line
                or "Problems reading a PDB file" in line
                or "Open Babel Warning" in line
                or "no charge value given" in line
                or "audit log messages" in line
                or "all CONECT" in line
                or "with zero occupancy" in line
                or "check:  Warnings:" in line
            ):
                pass
            else:
                print(line)
        elif "Total charge" in line:
            print("    ", line)
        elif "ERROR" in line.upper():
            if "Fatal error" in line:
                print(line, lines[count + 1])
            else:
                print(line)
        count += 1


def fixRes4Acpype(fpdb):
    code = fpdb[2]
    fLines = open(fpdb).readlines()
    famb = open(fpdb, "w")
    for line in fLines:
        if "ATOM  " in line:
            line = line[:17] + aa_dict[code].upper() + line[20:]
        famb.write(line)
    famb.close()


def build_residues_tleap():
    """Build residues tripeptides with tleap and minimise with sander"""
    mdin = """Minimization\n&cntrl\nimin=1, maxcyc=200, ntmin=2, ntb=0, igb=0,cut=999,/\n"""
    open("mdin", "w").writelines(mdin)
    seqi = aa_dict.keys()
    for aai in seqi:
        # if aai != 'H': continue
        aai3 = 3 * aai
        res = aa_dict[aai].upper()
        res1 = res
        leapDict = {"amberff": amberff, "res": res, "aai3": aai3, "res1": res1}
        tleapin = genPdbTemplate % leapDict
        open("tleap.in", "w").writelines(tleapin)
        cmd = "%s -f tleap.in" % (tleapExe)
        _getoutput(cmd)
        # cmd = "%s -O; %s < restrt > %s.pdb; mv mdinfo %s.mdinfo" % (sanderExe, ambpdbExe, aai3, aai3)
        # -i mdin -o mdout -p prmtop -c inpcrd" % (sanderExe)
        # tleap has already written aai3.pdb via savepdb; sander/ambpdb only replace
        # it with a minimised geometry, so skip that step when they are not installed.
        if sanderExe and ambpdbExe:
            cmd = f"{sanderExe} -O; {ambpdbExe} -c restrt > {aai3}.pdb"
            _getoutput(cmd)
    _getoutput("rm -f mdout mdinfo mdin restrt tleap.in prmtop inpcrd leap.log")


def error(v1, v2):
    """Percentage relative error"""
    if v1 == v2:
        return 0
    return abs(v1 - v2) / max(abs(v1), abs(v2)) * 100


def calcGmxPotEnerDiff(res):
    def getEnergies(template):
        cmd = template % cmdDict
        pdebug(cmd)
        out = _getoutput(cmd)
        out = out.split("\n")[-nEner:]
        dictEner = {}
        for item in out:
            # A gmx energy summary row is "<name words...> <average> <err.est> <rmsd>
            # <tot-drift> (<unit>)". Gluing the name with a fixed list of replacements
            # broke as soon as GROMACS renamed "Improper Dih." to "Per. Imp. Dih.", a
            # three word name, so the columns are counted from the right instead.
            tokens = item.split()
            if len(tokens) < 6:
                continue
            k = "_".join(tokens[:-5])
            v = float(tokens[-5])
            dictEner[k] = v
            if "Dih." in k or "Bell." in k:
                if "Dihedral P+I" in list(dictEner.keys()):
                    dictEner["Dihedral P+I"] = dictEner["Dihedral P+I"] + v
                else:
                    dictEner["Dihedral P+I"] = v
            if "Coulomb" in k or "LJ" in k:
                if "TotalNonBonded" in list(dictEner.keys()):
                    dictEner["TotalNonBonded"] = dictEner["TotalNonBonded"] + v
                else:
                    dictEner["TotalNonBonded"] = v
        dictEner["Total_Bonded"] = dictEner["Potential"] - dictEner["TotalNonBonded"]
        return dictEner

    os.chdir(tmpDir)
    nEner = 9  # number of energy entries from g_energy
    tEner = " ".join([str(x) for x in range(1, nEner + 1)])
    open("SPE.mdp", "w").writelines(spe_mdp)
    cmdDict = {
        "pdb2gmx": pdb2gmx,
        "grompp": grompp,
        "mdrun": mdrun,
        "res": res,
        "g_energy": g_energy,
        "tEner": tEner,
        "water": water,
        "gmxdump": gmxdump,
        "editconf": editconf,
        "gmxff": gmxFF,
    }

    # calc Pot Ener for aXXX.pdb (AMB_GMX)
    template = """%(pdb2gmx)s -ff %(gmxff)s -f a%(res)s.pdb -o a%(res)s_.pdb -p a%(res)s.top %(water)s
    %(editconf)s -f a%(res)s_.pdb -o a%(res)s_box.pdb -bt cubic -d 3.0
    %(grompp)s -c a%(res)s_box.pdb -p a%(res)s.top -f SPE.mdp -o a%(res)s.tpr -pp a%(res)sp.top
    %(mdrun)s -v -deffnm a%(res)s
    echo %(tEner)s | %(g_energy)s -f a%(res)s.edr
    """
    dictEnerAMB = getEnergies(template)
    # print(dictEnerAMB)

    # calc Pot Ener for agXXX.acpype/agXXX.pdb (ACPYPE_GMX)
    template = """
    %(editconf)s -f ag%(res)s.acpype/ag%(res)s_NEW.pdb -o ag%(res)s_box.pdb -bt cubic -d 3.0
    %(grompp)s -c ag%(res)s_box.pdb -p ag%(res)s.acpype/ag%(res)s_GMX.top -f SPE.mdp \
    -o ag%(res)s.tpr -pp ag%(res)sp.top
    %(mdrun)s -v -deffnm ag%(res)s
    echo %(tEner)s | %(g_energy)s -f ag%(res)s.edr
    """
    dictEnerACPYPE = getEnergies(template)
    # print(dictEnerACPYPE)

    order = [
        "LJ-14",
        "LJ_(SR)",
        "Coulomb-14",
        "Coulomb_(SR)",
        "TotalNonBonded",
        "Potential",
        "Angle",
        "Bond",
        "Proper_Dih.",
        "Improper_Dih.",
        "Dihedral P+I",
        "Total_Bonded",
    ]  # sorted(dictEnerAMB.items())
    for k in order:
        v = dictEnerAMB.get(k)
        v2 = dictEnerACPYPE.get(k)
        rerror = error(v2, v)
        if rerror > 0.1:
            print("%15s   %.3f   %5.6f x %5.6f" % (k, rerror, v2, v))
        else:
            print("%15s   %.3f" % (k, rerror))

    cmd = "{gmxdump} -s a{res}.tpr".format(**cmdDict)
    amb = _getoutput(cmd)
    cmd = "{gmxdump} -s ag{res}.tpr".format(**cmdDict)
    acp = _getoutput(cmd)

    # dihAmb = [x.split(']=')[-1] for x in amb.splitlines() if ('PIDIHS' in x or 'PDIHS' in x) and 'functype' in x ]
    # dihAcp = [x.split(']=')[-1] for x in acp.splitlines() if ('PIDIHS' in x or 'PDIHS' in x) and 'functype' in x ]

    dihAmb = [x.split("PIDIHS")[-1][2:] for x in amb.splitlines() if ("PIDIHS" in x)]
    dihAcp = [x.split("PIDIHS")[-1][2:] for x in acp.splitlines() if ("PIDIHS" in x)]

    dAcp = set(dihAcp).difference(set(dihAmb))
    dAmb = set(dihAmb).difference(set(dihAcp))
    rAcp = [" ".join(reversed(x.split())) for x in dAcp]
    rAmb = [" ".join(reversed(x.split())) for x in dAmb]
    ddAmb = sorted(dAmb.difference(rAcp))
    ddAcp = sorted(dAcp.difference(rAmb))
    if ddAmb:
        print("IDIH: Amb diff Acp", ddAmb)
    if ddAcp:
        print("IDIH: Acp diff Amb", ddAcp)

    return dihAmb, dihAcp


if __name__ == "__main__":
    """order: (tleap/EM or pymol) AAA.pdb -f-> _AAA.pdb -f-> aAAA.pdb
    --> (pdb2gmx) agAAA.pdb -f-> agAAA.pdb --> acpype"""
    aNb = nbDict(open(gmxTopDir + f"/gromacs/top/{gmxFF}.ff/ffnonbonded.itp").readlines())
    aBon = parseTopFile(open(gmxTopDir + f"/gromacs/top/{gmxFF}.ff/ffbonded.itp").readlines())

    tmpFile = "tempScript.py"
    if not os.path.exists(tmpDir):
        os.mkdir(tmpDir)
    os.chdir(tmpDir)
    subprocess.run(r"rm -fr \#* *.acpype", shell=True, check=False)
    # create res.pdb
    if usePymol:
        ff = open(tmpFile, "w")
        ff.writelines(pymolScript)
        ff.close()
        cmd = f"{exePymol} -qc {tmpFile}"
        subprocess.run(cmd, shell=True, check=False)
    else:
        build_residues_tleap()
    listRes = os.listdir(".")  # list all res.pdb
    listRes = glob("???.pdb")
    listRes.sort()
    failures = []
    for resFile in listRes:
        res, ext = os.path.splitext(resFile)  # eg. res = 'AAA'
        # if res != 'RRR': continue
        if len(resFile) == 7 and ext == ".pdb" and resFile[:3].isupper():
            print(f"\nFile {resFile} : residue {aa_dict[res[0]].upper()}")
            # ACPYPE and calcGmxPotEnerDiff both chdir, so a residue that fails part
            # way through would otherwise strand the loop in its output directory and
            # take every later residue down with it.
            os.chdir(tmpDir)
            try:
                _pdb = createOldPdb2(resFile)  # using my own dict
                apdb = createAmberPdb3(_pdb)  # create file aAAA.pdb with NXXX, CXXX

                # from ogpdb to amber pdb and top
                agpdb = "ag%s.pdb" % res  # output name
                agtop = "ag%s.top" % res
                cmd = f" {pdb2gmx} -f {apdb} -o {agpdb} -p {agtop} -ff {gmxFF} {water}"
                pdebug(cmd)
                out = _getoutput(cmd)
                # print(out)
                # parseOut(out)

                # acpype on agpdb file
                fixRes4Acpype(agpdb)
                cv = None
                # cmd = "%s -dfi %s -c %s -a %s" % (acpypeExe, agpdb, cType, ffType)
                if res in ["JJJ", "RRR"] and cType == "bcc":
                    cv = 3  # cmd += ' -n 3' # acpype failed to get correct charge
                mol = ACTopol(
                    agpdb,
                    chargeType=cType,
                    atomType=ffType,
                    protein_ff=proteinFF,
                    chargeVal=cv,
                    debug=False,
                    verbose=debug,
                )
                mol.createACTopol()
                mol.createMolTopol()

                # out = commands.getstatusoutput(cmd)
                # parseOut(out[1])
                print(f"Compare ACPYPE x GMX {gmxFF.upper()} topol & param")
                checkTopAcpype(res)
                # calc Pot Energy
                print(f"Compare ACPYPE x GMX {gmxFF.upper()} Pot. Energy (|ERROR|%)")
                # calc gmx amber energies
                dihAmb, dihAcp = calcGmxPotEnerDiff(res)
            # One bad residue must not end a run that takes minutes.
            except Exception as exc:
                failures.append((res, f"{type(exc).__name__}: {exc}"))
                perror(f"{res} FAILED: {type(exc).__name__}: {exc}")
            # calc gmx acpype energies

    print(f"\n=== {len(listRes) - len(failures)}/{len(listRes)} residues compared ===")
    if failures:
        print("Failed:")
        for res, why in failures:
            print(f"  {res}: {why}")
        if not (sanderExe and ambpdbExe):
            print(
                "\nNote: sander/ambpdb were not found, so the tripeptides were not energy\n"
                "minimised. tleap's raw geometry is often poor enough that antechamber\n"
                "mis-perceives aromatic rings (atom type DU), which is the usual cause of\n"
                "the failures above. Install a full AmberTools for a faithful comparison:\n"
                "    conda install -c conda-forge ambertools"
            )
    print(f"\nWorking directory: {tmpDir}")

    subprocess.run(r"rm -f %s/\#* posre.itp tempScript.py" % tmpDir, shell=True, check=False)
    subprocess.run(
        "find . -name 'ag*GMX*.itp' | xargs grep -v 'created by acpype on' > standard_ag_itp.txt",
        shell=True,
        check=False,
    )
