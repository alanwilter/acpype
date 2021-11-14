#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Requirements:
* pytest
* pytest-html

To run:
pytest -s --html=report.html

"""
import os
import shutil
from acpype_lib.acpype import ACTopol, MolTopol


def test_mol2():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    molecule = ACTopol("AAA.mol2", chargeType="gas", debug=True)
    molecule.createACTopol()
    molecule.createMolTopol()
    assert molecule
    assert len(molecule.molTopol.atoms) == 33
    assert len(molecule.molTopol.properDihedrals) == 95
    assert len(molecule.molTopol.improperDihedrals) == 5
    assert molecule.molTopol.totalCharge == 0
    shutil.rmtree(molecule.absHomeDir)


def test_pdb():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    molecule = ACTopol("FFF.mol2", chargeType="gas", debug=True)
    molecule.createACTopol()
    molecule.createMolTopol()
    assert molecule
    assert len(molecule.molTopol.atoms) == 63
    assert len(molecule.molTopol.properDihedrals) == 185
    assert len(molecule.molTopol.improperDihedrals) == 23
    shutil.rmtree(molecule.absHomeDir)


def test_charges():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    molecule = ACTopol("KKK.mol2", chargeType="gas", debug=True)
    molecule.createACTopol()
    molecule.createMolTopol()
    assert molecule
    assert len(molecule.molTopol.atoms) == 69
    assert len(molecule.molTopol.properDihedrals) == 215
    assert len(molecule.molTopol.improperDihedrals) == 5
    assert len(molecule.molTopol.chiralGroups) == 3
    assert molecule.chargeVal == "3"
    assert molecule.molTopol.totalCharge == 3
    shutil.rmtree(molecule.absHomeDir)


def test_smiles():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    smiles = "c1ccc2c(c1)C(=O)N(C2=O)C3CCC(=NC3=O)O"
    molecule = ACTopol(smiles, basename="thalidomide_smiles", chargeType="gas", debug=True)
    molecule.createACTopol()
    molecule.createMolTopol()
    assert molecule
    assert len(molecule.molTopol.atoms) == 29
    shutil.rmtree(molecule.absHomeDir)
    os.remove(molecule.absInputFile)


def test_glycam():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    molecule = MolTopol(acFileTop="glycam_exe.prmtop", acFileXyz="glycam_exe.inpcrd", debug=True, amb2gmx=True)
    molecule.writeGromacsTopolFiles()
    assert molecule
    assert molecule.topo14Data.hasNondefault14()
    assert len(molecule.topo14Data.scnb_scale_factor) == 31
    shutil.rmtree(molecule.absHomeDir)
