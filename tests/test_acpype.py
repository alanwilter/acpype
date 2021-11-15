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
import pytest
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
    assert molecule.molTopol.atoms[-1].__repr__() == "<Atom id=33, name=H8, <AtomType=hc>>"
    # check sorted
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    molecule = ACTopol("AAA.mol2", chargeType="gas", debug=True, force=True, is_sorted=True)
    molecule.createACTopol()
    molecule.createMolTopol()
    assert molecule
    assert molecule.molTopol.atoms[-1].__repr__() == "<Atom id=33, name=OXT, <AtomType=o>>"
    shutil.rmtree(molecule.absHomeDir)


def test_pdb():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    molecule = ACTopol("FFF.pdb", chargeType="gas", debug=True)
    molecule.createACTopol()
    molecule.createMolTopol()
    assert molecule
    assert len(molecule.molTopol.atoms) == 63
    assert len(molecule.molTopol.properDihedrals) == 185
    assert len(molecule.molTopol.improperDihedrals) == 23
    assert molecule.molTopol.atoms[0].__repr__() == "<Atom id=1, name=N, <AtomType=n4>>"
    # check gaff2 and force
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    molecule = ACTopol("FFF.pdb", chargeType="gas", debug=True, atomType="gaff2", force=True)
    molecule.createACTopol()
    molecule.createMolTopol()
    assert molecule
    assert len(molecule.molTopol.atoms) == 63
    assert len(molecule.molTopol.properDihedrals) == 188
    assert len(molecule.molTopol.improperDihedrals) == 23
    assert molecule.molTopol.atoms[0].__repr__() == "<Atom id=1, name=N, <AtomType=nz>>"
    # check for already present
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    molecule = ACTopol("FFF.mol2", chargeType="gas", debug=True)
    molecule.createACTopol()
    molecule.createMolTopol()
    assert molecule
    shutil.rmtree(molecule.absHomeDir)


def test_amber():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    molecule = ACTopol("FFF.mol2", chargeType="gas", debug=True, atomType="amber")
    molecule.createACTopol()
    molecule.createMolTopol()
    assert molecule
    assert len(molecule.molTopol.atoms) == 63
    assert len(molecule.molTopol.properDihedrals) == 189
    assert len(molecule.molTopol.improperDihedrals) == 23
    assert molecule.molTopol.atoms[0].__repr__() == "<Atom id=1, name=N, <AtomType=N3>>"
    # check amber2
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    molecule = ACTopol("FFF.mol2", chargeType="gas", debug=True, atomType="amber2", force=True)
    molecule.createACTopol()
    molecule.createMolTopol()
    assert molecule
    assert len(molecule.molTopol.atoms) == 63
    assert len(molecule.molTopol.properDihedrals) == 187
    assert len(molecule.molTopol.improperDihedrals) == 23
    assert molecule.molTopol.atoms[0].__repr__() == "<Atom id=1, name=N, <AtomType=N3>>"
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


def test_amb2gmx():
    # oct box with water and ions
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    molecule = MolTopol(acFileTop="RAMP1_ion.prmtop", acFileXyz="RAMP1_ion.inpcrd", debug=True, amb2gmx=True)
    molecule.writeGromacsTopolFiles()
    assert molecule
    assert not molecule.topo14Data.hasNondefault14()
    assert molecule.atoms[1300].__repr__() == "<Atom id=1301, name=NA+, <AtomType=Na+>>"
    assert molecule.atoms[1310].__repr__() == "<Atom id=1311, name=CL-, <AtomType=Cl->>"
    assert len(molecule.atoms) == 18618
    shutil.rmtree(molecule.absHomeDir)


def test_sqm_tleap():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    molecule = ACTopol("benzene.pdb", chargeType="bcc", debug=True)
    molecule.createACTopol()
    molecule.createMolTopol()
    assert molecule
    assert len(molecule.molTopol.atoms) == 12
    assert len(molecule.molTopol.properDihedrals) == 24
    assert len(molecule.molTopol.improperDihedrals) == 6
    shutil.rmtree(molecule.absHomeDir)


def test_time_limit():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    molecule = ACTopol("KKK.pdb", chargeType="bcc", debug=True, timeTol=2)
    with pytest.raises(Exception) as e_info:
        molecule.createACTopol()
    assert e_info.value.args[0] == 'Semi-QM taking too long to finish... aborting!'
    shutil.rmtree(molecule.absHomeDir)