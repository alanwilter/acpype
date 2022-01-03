import io
import json
import os
import shutil
import sys
import time
import traceback

from acpype.params import MAXTIME
from acpype.topol import ACTopol, header
from acpype.utils import elapsedTime, while_replace

em_mdp = io.StringIO()
AC_frcmod = io.StringIO()
AC_inpcrd = io.StringIO()
AC_lib = io.StringIO()
AC_prmtop = io.StringIO()
mol2 = io.StringIO()
CHARMM_inp = io.StringIO()
CHARMM_prm = io.StringIO()
CHARMM_rtf = io.StringIO()
CNS_inp = io.StringIO()
CNS_par = io.StringIO()
CNS_top = io.StringIO()
GMX_OPLS_itp = io.StringIO()
GMX_OPLS_top = io.StringIO()
GMX_gro = io.StringIO()
GMX_itp = io.StringIO()
GMX_top = io.StringIO()
NEW_pdb = io.StringIO()
md_mdp = io.StringIO()

filesInMemory = [
    (em_mdp, "em.mdp"),
    (AC_frcmod, "_AC.frcmod"),
    (AC_inpcrd, "_AC.inpcrd"),
    (AC_lib, "_AC.lib"),
    (AC_prmtop, "_AC.prmtop"),
    (mol2, ".mol2"),
    (CHARMM_inp, "_CHARMM.inp"),
    (CHARMM_prm, "_CHARMM.prm"),
    (CHARMM_rtf, "_CHARMM.rtf"),
    (CNS_inp, "_CNS.inp"),
    (CNS_par, "_CNS.par"),
    (CNS_top, "_CNS.top"),
    (GMX_OPLS_itp, "_GMX_OPLS.itp"),
    (GMX_OPLS_top, "_GMX_OPLS.top"),
    (GMX_gro, "_GMX.gro"),
    (GMX_itp, "_GMX.itp"),
    (GMX_top, "_GMX.top"),
    (NEW_pdb, "_NEW.pdb"),
    (md_mdp, "md.mdp"),
]


def clearFileInMemory():
    for files in filesInMemory:
        files[0].seek(0)
        files[0].truncate(0)


def readFiles(basename, chargeType, atomType):
    for files in filesInMemory:
        if files[1] == "em.mdp" or files[1] == "md.mdp":
            filename = files[1]
        elif files[1] == ".mol2":
            filename = basename + "_" + chargeType + "_" + atomType + files[1]
        else:
            filename = basename + files[1]
        readfile = tuple(open(filename))
        for line in readfile:
            files[0].write(line)


def acpype_api(
    inputFile,
    chargeType="bcc",
    chargeVal=None,
    multiplicity="1",
    atomType="gaff2",
    force=False,
    basename=None,
    debug=False,
    outTopol="all",
    engine="tleap",
    allhdg=False,
    timeTol=MAXTIME,
    qprog="sqm",
    ekFlag=None,
    verbose=True,
    gmx4=False,
    merge=False,
    direct=False,
    is_sorted=False,
    chiral=False,
    is_smiles=False,
):

    at0 = time.time()
    print(header)

    if debug:
        texta = "Python Version %s" % sys.version
        print("DEBUG: %s" % while_replace(texta))
    try:
        molecule = ACTopol(
            inputFile=inputFile,
            chargeType=chargeType,
            chargeVal=chargeVal,
            debug=debug,
            multiplicity=multiplicity,
            atomType=atomType,
            force=force,
            outTopol=outTopol,
            allhdg=allhdg,
            basename=basename,
            timeTol=timeTol,
            qprog=qprog,
            ekFlag=ekFlag,
            verbose=verbose,
            gmx4=gmx4,
            merge=merge,
            direct=direct,
            is_sorted=is_sorted,
            chiral=chiral,
        )

        molecule.createACTopol()
        molecule.createMolTopol()

        "Output in JSON format"
        os.chdir(molecule.absHomeDir)
        readFiles(molecule.baseName, chargeType, atomType)
        output = {
            "file_name": molecule.baseName,
            "em_mdp": em_mdp.getvalue(),
            "AC_frcmod": AC_frcmod.getvalue(),
            "AC_inpcrd": AC_inpcrd.getvalue(),
            "AC_lib": AC_lib.getvalue(),
            "AC_prmtop": AC_prmtop.getvalue(),
            "mol2": mol2.getvalue(),
            "CHARMM_inp": CHARMM_inp.getvalue(),
            "CHARMM_prm": CHARMM_prm.getvalue(),
            "CHARMM_rtf": CHARMM_rtf.getvalue(),
            "CNS_inp": CNS_inp.getvalue(),
            "CNS_par": CNS_par.getvalue(),
            "CNS_top": CNS_top.getvalue(),
            "GMX_OPLS_itp": GMX_OPLS_itp.getvalue(),
            "GMX_OPLS_top": GMX_OPLS_top.getvalue(),
            "GMX_gro": GMX_gro.getvalue(),
            "GMX_itp": GMX_itp.getvalue(),
            "GMX_top": GMX_top.getvalue(),
            "NEW_pdb": NEW_pdb.getvalue(),
            "md_mdp": md_mdp.getvalue(),
        }

    except Exception:
        _exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
        print("ACPYPE FAILED: %s" % exceptionValue)
        if debug:
            traceback.print_tb(exceptionTraceback, file=sys.stdout)
            output = {"file_name": f"ERROR: {str(exceptionValue)}"}

    execTime = int(round(time.time() - at0))
    if execTime == 0:
        amsg = "less than a second"
    else:
        amsg = elapsedTime(execTime)
    print("Total time of execution: %s" % amsg)
    clearFileInMemory()
    try:
        shutil.rmtree(molecule.absHomeDir)
    except Exception:
        print("DEBUG: No folder left to be removed")
    return json.dumps(output)
