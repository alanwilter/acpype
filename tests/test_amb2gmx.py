import os
import shutil
import pytest
from acpype.topol import MolTopol


def test_glycam():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    molecule = MolTopol(acFileTop="glycam_exe.prmtop", acFileXyz="glycam_exe.inpcrd", debug=True, amb2gmx=True)
    molecule.writeGromacsTopolFiles()
    assert molecule
    assert molecule.topo14Data.hasNondefault14()
    assert len(molecule.topo14Data.scnb_scale_factor) == 31
    shutil.rmtree(molecule.absHomeDir)


@pytest.mark.parametrize(
    ("dd", "g4", "ntext"),
    [(False, False, 14193), (True, False, 31516), (False, True, 12124), (True, True, 29447)],
)
def test_amb2gmx(dd, g4, ntext):
    # oct box with water and ions
    # modified from https://ambermd.org/tutorials/basic/tutorial7/index.php
    # using addIonsRand separated for each ion and TIP3PBOX
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    molecule = MolTopol(acFileTop="RAMP1_ion.prmtop", acFileXyz="RAMP1_ion.inpcrd", debug=True, direct=dd, gmx4=g4)
    molecule.writeGromacsTopolFiles()
    assert molecule
    assert len(molecule.topText) == ntext
    assert not molecule.topo14Data.hasNondefault14()
    assert molecule.atoms[1300].__repr__() == "<Atom id=1301, name=NA+, <AtomType=Na+>>"
    assert molecule.atoms[1310].__repr__() == "<Atom id=1311, name=CL-, <AtomType=Cl->>"
    assert len(molecule.atoms) == 18618
    shutil.rmtree(molecule.absHomeDir)


@pytest.mark.parametrize(
    ("merge", "gaff", "n_at", "acoef", "bcoef", "msg"),
    [
        (False, "1", 50, 379876.399, 564.885984, "<AtomType=o>"),
        (False, "2", 50, 376435.47, 469.350655, "<AtomType=o>"),
        (True, "1", 41, 361397.723, 495.732238, "<AtomType=os>"),
        (True, "2", 50, 376435.47, 469.350655, "<AtomType=o_>"),
    ],
)
def test_merge(merge, gaff, n_at, acoef, bcoef, msg):
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    molecule = MolTopol(acFileTop=f"ComplexG{gaff}.prmtop", acFileXyz=f"ComplexG{gaff}.inpcrd", debug=True, merge=merge)
    molecule.writeGromacsTopolFiles()
    assert molecule
    assert len(molecule.atomTypesGromacs) == n_at
    assert molecule.atomTypesGromacs[31].ACOEF == acoef
    assert molecule.atomTypesGromacs[31].BCOEF == bcoef
    assert molecule.atomTypesGromacs[31].__repr__() == msg
    shutil.rmtree(molecule.absHomeDir)


@pytest.mark.parametrize(
    ("mol", "n1", "n2", "n3", "n4", "n5", "msg"),
    [("ILDN", 24931, 12230, 577, 0, 47, "<AtomType=C6>"), ("Base", 24931, 12044, 577, 0, 43, "<AtomType=HS>")],
)
def test_ildn(mol, n1, n2, n3, n4, n5, msg):
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    molecule = MolTopol(acFileTop=f"{mol}.prmtop", acFileXyz=f"{mol}.inpcrd", debug=True)
    molecule.writeGromacsTopolFiles()
    assert molecule
    assert len(molecule.atoms) == n1
    assert len(molecule.properDihedrals) == n2
    assert len(molecule.improperDihedrals) == n3
    assert molecule.totalCharge == n4
    assert len(molecule.atomTypes) == n5
    assert molecule.atomTypes[23].__repr__() == msg
    shutil.rmtree(molecule.absHomeDir)


def test_ildn_gmx4_fail():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    molecule = MolTopol(acFileTop="ILDN.prmtop", acFileXyz="ILDN.inpcrd", debug=True, gmx4=True)
    with pytest.raises(Exception) as e_info:
        molecule.writeGromacsTopolFiles()
    assert e_info.value.args[0] == "Likely trying to convert ILDN to RB, DO NOT use option '-z'"
    shutil.rmtree(molecule.absHomeDir)


def test_wrong_input():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    with pytest.raises(Exception) as e_info:
        MolTopol(acFileTop="nope.prmtop", acFileXyz="ILDN.inpcrd", debug=True, gmx4=True)
    assert e_info.value.args == (2, "No such file or directory")
