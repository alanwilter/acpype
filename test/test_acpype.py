#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Requirements:
* pytest
* pytest-html

To run:
pytest -s --html=report.html

"""

import sys
import os
import traceback
import time
from acpype_lib.acpype import ACTopol, MAXTIME, while_replace, header, elapsedTime


def run_acpype(
    inputFile,
    chargeType="gas",
    chargeVal=None,
    multiplicity="1",
    atomType="gaff",
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
    disam=False,
    direct=False,
    is_sorted=False,
    chiral=False,
    InMemory=True,
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
            engine=engine,
            allhdg=allhdg,
            basename=basename,
            timeTol=timeTol,
            qprog=qprog,
            ekFlag=ekFlag,
            verbose=verbose,
            gmx4=gmx4,
            disam=disam,
            direct=direct,
            is_sorted=is_sorted,
            chiral=chiral,
            InMemory=InMemory,
        )

        molecule.createACTopol()
        molecule.createMolTopol()
        # if not basename:
        #     file_name = "smiles_molecule"
        # else:
        #     file_name = basename

    except Exception:
        _exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
        print("ACPYPE FAILED: %s" % exceptionValue)
        if debug:
            traceback.print_tb(exceptionTraceback, file=sys.stdout)
        return 1

    execTime = int(round(time.time() - at0))
    if execTime == 0:
        amsg = "less than a second"
    else:
        amsg = elapsedTime(execTime)
    print("Total time of execution: %s" % amsg)
    return 0


def acpype_test(input, name):
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    exec_acpype = run_acpype(inputFile=input, basename=name)
    assert exec_acpype == 0


def test_mol2(input="AAA.mol2"):
    acpype_test(input, name=None)


def test_pdb(input="AAA.pdb"):
    acpype_test(input, name=None)


def test_chiral(input="thalidomide.mol2"):
    acpype_test(input, name=None)


def test_smiles(input="c1ccc2c(c1)C(=O)N(C2=O)C3CCC(=NC3=O)O"):
    acpype_test(input, name="thalidomide_smiles")
    try:
        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        os.remove("smiles_molecule.mol2")
    except Exception as ex:
        print(ex)
