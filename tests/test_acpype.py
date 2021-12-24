import os
import shutil
import pytest
from pytest import approx
from acpype import __version__ as version
from acpype.cli import init_main
from acpype.topol import ACTopol
from acpype.utils import _getoutput


@pytest.mark.parametrize(
    ("issorted", "charge", "msg"),
    [
        (False, 0.240324, "<Atom id=33, name=H8, <AtomType=hc>>"),
        (True, 0.240324, "<Atom id=33, name=OXT, <AtomType=o>>"),
    ],
)
def test_mol2_sorted(issorted, charge, msg):
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    molecule = ACTopol("AAA.mol2", chargeType="gas", debug=True, is_sorted=issorted)
    molecule.createACTopol()
    molecule.createMolTopol()
    assert molecule
    assert molecule.molTopol.atomTypes[0].__repr__() == "<AtomType=nz>"
    assert len(molecule.molTopol.atoms) == 33
    assert len(molecule.molTopol.properDihedrals) == 98
    assert len(molecule.molTopol.improperDihedrals) == 5
    assert molecule.molTopol.totalCharge == 0
    assert molecule.molTopol.atoms[0].charge == approx(charge)
    assert molecule.molTopol.atoms[-1].__repr__() == msg
    shutil.rmtree(molecule.absHomeDir)


def test_pdb(capsys):
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    molecule = ACTopol("FFF.pdb", chargeType="gas", debug=True)
    molecule.createACTopol()
    molecule.createMolTopol()
    assert len(molecule.molTopol.atoms) == 63
    assert len(molecule.molTopol.properDihedrals) == 188
    assert len(molecule.molTopol.improperDihedrals) == 23
    assert molecule.molTopol.atoms[0].__repr__() == "<Atom id=1, name=N, <AtomType=nz>>"
    # check gaff2 and force
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    molecule = ACTopol("FFF.pdb", chargeType="gas", debug=True, atomType="gaff2", force=True)
    molecule.createACTopol()
    molecule.createMolTopol()
    captured = capsys.readouterr()
    assert len(molecule.molTopol.atoms) == 63
    assert len(molecule.molTopol.properDihedrals) == 188
    assert len(molecule.molTopol.improperDihedrals) == 23
    assert molecule.molTopol.atoms[0].__repr__() == "<Atom id=1, name=N, <AtomType=nz>>"
    assert "==> Overwriting pickle file FFF.pkl" in captured.out
    # check for already present
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    molecule = ACTopol("FFF.mol2", chargeType="gas", debug=True)
    molecule.createACTopol()
    molecule.createMolTopol()
    captured = capsys.readouterr()
    assert molecule
    assert "==> Pickle file FFF.pkl already present... doing nothing" in captured.out
    shutil.rmtree(molecule.absHomeDir)


@pytest.mark.parametrize(
    ("force", "at", "ndih"),
    [(False, "amber", 189), (True, "amber2", 187)],
)
def test_amber(force, at, ndih):
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    molecule = ACTopol("FFF.mol2", chargeType="gas", debug=True, atomType=at, force=force)
    molecule.createACTopol()
    molecule.createMolTopol()
    assert molecule
    assert len(molecule.molTopol.atoms) == 63
    assert len(molecule.molTopol.properDihedrals) == ndih
    assert len(molecule.molTopol.improperDihedrals) == 23
    assert molecule.molTopol.atoms[0].__repr__() == "<Atom id=1, name=N, <AtomType=N3>>"
    shutil.rmtree(molecule.absHomeDir)


def test_charges_chiral():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    molecule = ACTopol("KKK.mol2", chargeType="gas", debug=True, chiral=True)
    molecule.createACTopol()
    molecule.createMolTopol()
    assert molecule
    assert len(molecule.molTopol.atoms) == 69
    assert len(molecule.molTopol.properDihedrals) == 218
    assert len(molecule.molTopol.improperDihedrals) == 5
    assert len(molecule.molTopol.chiralGroups) == 3
    assert molecule.chargeVal == "3"
    assert molecule.molTopol.totalCharge == 3
    assert molecule.molTopol.chiralGroups[-1][-1] == approx(66.713290)
    shutil.rmtree(molecule.absHomeDir)


@pytest.mark.parametrize(
    ("base", "msg"),
    [(None, "smiles_molecule.mol2"), ("thalidomide_smiles", "thalidomide_smiles.mol2")],
)
def test_smiles(base, msg):
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    smiles = "c1ccc2c(c1)C(=O)N(C2=O)C3CCC(=NC3=O)O"
    molecule = ACTopol(smiles, basename=base, chargeType="gas", debug=True)
    molecule.createACTopol()
    molecule.createMolTopol()
    assert molecule
    assert molecule.inputFile == msg
    assert len(molecule.molTopol.atoms) == 29
    shutil.rmtree(molecule.absHomeDir)
    os.remove(molecule.absInputFile)


@pytest.mark.parametrize(
    ("ct", "msg"),
    [
        (
            "bcc",
            "-dr no -i benzene.mol2 -fi mol2 -o benzene_bcc_gaff2.mol2 -fo mol2 -c bcc -nc 0 -m 1 -s 2 -df 2 -at gaff2",
        ),
        ("user", "cannot read charges from a PDB file"),
    ],
)
def test_sqm_tleap(capsys, ct, msg):
    # check chargeType user with PDB -> use bcc
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    molecule = ACTopol("benzene.pdb", chargeType=ct, debug=True)
    molecule.createACTopol()
    molecule.createMolTopol()
    captured = capsys.readouterr()
    assert molecule
    assert len(molecule.molTopol.atoms) == 12
    assert len(molecule.molTopol.properDihedrals) == 24
    assert len(molecule.molTopol.improperDihedrals) == 6
    assert molecule.molTopol.atoms[0].charge == approx(-0.13)
    assert molecule.molTopol.atoms[-1].charge == approx(0.13)
    assert msg in captured.out
    shutil.rmtree(molecule.absHomeDir)


def test_time_limit():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    molecule = ACTopol("KKK.pdb", chargeType="bcc", debug=True, timeTol=2)
    with pytest.raises(Exception) as e_info:
        molecule.createACTopol()
    assert e_info.value.args[0] == "Semi-QM taking too long to finish... aborting!"
    shutil.rmtree(molecule.absHomeDir)


def test_wrong_element(capsys):
    # Only elements are allowed: C, N, O, S, P, H, F, Cl, Br and I
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    with pytest.raises(Exception) as e_info:
        ACTopol("HEM.mol2", chargeType="user", debug=True)
    captured = capsys.readouterr()
    assert e_info.typename == "FileNotFoundError"
    assert "Unrecognized case-sensitive atomic symbol ( FE)." in captured.out


def test_charge_user():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    molecule = ACTopol("ADPMg.mol2", chargeType="user", debug=True)
    molecule.createACTopol()
    molecule.createMolTopol()
    assert molecule
    assert len(molecule.molTopol.atoms) == 39
    assert len(molecule.molTopol.properDihedrals) == 128
    assert len(molecule.molTopol.improperDihedrals) == 7
    assert molecule.molTopol.atoms[0].charge == 0.1667
    assert molecule.molTopol.atoms[15].charge == -0.517
    shutil.rmtree(molecule.absHomeDir)


@pytest.mark.parametrize(
    ("argv"),
    [["-di", "cccc"], ["-x", "Base.inpcrd", "-p", "Base.prmtop"]],
)
def test_inputs(capsys, argv):
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    temp_base = "vir_temp"
    init_main(argv=argv + ["-b", temp_base])
    captured = capsys.readouterr()
    assert "Total time of execution:" in captured.out
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    _getoutput(f"rm -fr {temp_base}*")


@pytest.mark.parametrize(
    ("argv", "code", "msg"),
    [
        (None, 2, " error: "),  # NOTE: None -> sys.argv from pytest
        (["-v"], 0, version),
        ([], 2, "error: missing input files"),
        (["-d", "-w"], 2, "error: argument -w/--verboseless: not allowed with argument -d/--debug"),
        (["-di", "AAAx.mol2"], 19, "ACPYPE FAILED: Input file AAAx.mol2 DOES NOT EXIST"),
        (["-di", " 123"], 19, "ACPYPE FAILED: [Errno 2] No such file or directory"),
        (["-di", " 123", "-x", "abc"], 2, "either '-i' or ('-p', '-x'), but not both"),
        (["-di", " 123", "-u"], 2, "option -u is only meaningful in 'amb2gmx' mode (args '-p' and '-x')"),
    ],
)
def test_args_wrong_inputs(capsys, argv, code, msg):
    with pytest.raises(SystemExit) as e_info:
        init_main(argv=argv)
    captured = capsys.readouterr()
    assert msg in captured.err + captured.out
    assert e_info.typename == "SystemExit"
    assert e_info.value.code == code
