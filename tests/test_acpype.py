import os

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
def test_mol2_sorted(janitor, issorted, charge, msg):
    molecule = ACTopol("AAA.mol2", chargeType="gas", debug=True, is_sorted=issorted)
    molecule.createACTopol()
    molecule.createMolTopol()
    assert molecule
    assert molecule.molTopol.atomTypes[0].__repr__() == "<AtomType=nz>"
    assert len(molecule.molTopol.atoms) == 33
    assert len(molecule.molTopol.properDihedrals) == 91
    assert len(molecule.molTopol.improperDihedrals) == 5
    assert molecule.molTopol.totalCharge == 0
    assert molecule.molTopol.atoms[0].charge == approx(charge)
    assert molecule.molTopol.atoms[-1].__repr__() == msg
    janitor.append(molecule.absHomeDir)
    janitor.append(molecule.tmpDir)


def test_pdb(janitor, capsys):
    molecule = ACTopol("FFF.pdb", chargeType="gas", debug=True, force=False)
    molecule.createACTopol()
    molecule.createMolTopol()
    assert len(molecule.molTopol.atoms) == 63
    assert len(molecule.molTopol.properDihedrals) == 181
    assert len(molecule.molTopol.improperDihedrals) == 23
    assert molecule.molTopol.atoms[0].__repr__() == "<Atom id=1, name=N, <AtomType=nz>>"
    # check gaff2 and force
    molecule = ACTopol("FFF.pdb", chargeType="gas", debug=True, atomType="gaff2", force=True)
    molecule.createACTopol()
    molecule.createMolTopol()
    captured = capsys.readouterr()
    assert len(molecule.molTopol.atoms) == 63
    assert len(molecule.molTopol.properDihedrals) == 181
    assert len(molecule.molTopol.improperDihedrals) == 23
    assert molecule.molTopol.atoms[0].__repr__() == "<Atom id=1, name=N, <AtomType=nz>>"
    assert "==> Overwriting pickle file FFF.pkl" in captured.out
    # check for already present
    molecule = ACTopol("FFF.mol2", chargeType="gas", debug=True)
    molecule.createACTopol()
    molecule.createMolTopol()
    captured = capsys.readouterr()
    assert molecule
    assert "==> Pickle file FFF.pkl already present... doing nothing" in captured.out
    janitor.append(molecule.absHomeDir)
    janitor.append(molecule.tmpDir)


@pytest.mark.parametrize(
    ("force", "at", "ndih"),
    [(False, "amber", 189), (True, "amber2", 197)],
)
def test_amber(janitor, force, at, ndih):
    molecule = ACTopol("FFF.mol2", chargeType="gas", debug=True, atomType=at, force=force)
    molecule.createACTopol()
    molecule.createMolTopol()
    assert molecule
    assert len(molecule.molTopol.atoms) == 63
    assert len(molecule.molTopol.properDihedrals) == ndih
    assert len(molecule.molTopol.improperDihedrals) == 23
    assert molecule.molTopol.atoms[0].__repr__() == "<Atom id=1, name=N, <AtomType=N3>>"
    janitor.append(molecule.absHomeDir)
    janitor.append(molecule.tmpDir)


def test_charges_chiral(janitor):
    molecule = ACTopol("KKK.mol2", chargeType="gas", debug=True, chiral=True)
    molecule.createACTopol()
    molecule.createMolTopol()
    assert molecule
    assert len(molecule.molTopol.atoms) == 69
    assert len(molecule.molTopol.properDihedrals) == 205
    assert len(molecule.molTopol.improperDihedrals) == 5
    assert len(molecule.molTopol.chiralGroups) == 3
    assert molecule.chargeVal == "3"
    assert molecule.molTopol.totalCharge == 3
    assert molecule.molTopol.chiralGroups[-1][-1] == approx(66.713290)
    janitor.append(molecule.absHomeDir)
    janitor.append(molecule.tmpDir)


@pytest.mark.parametrize(
    ("base", "msg"),
    [(None, "smiles_molecule.mol2"), ("thalidomide_smiles", "thalidomide_smiles.mol2")],
)
def test_smiles(janitor, base, msg):
    smiles = "c1ccc2c(c1)C(=O)N(C2=O)C3CCC(=NC3=O)O"
    molecule = ACTopol(smiles, basename=base, chargeType="gas", debug=True)
    molecule.createACTopol()
    molecule.createMolTopol()
    assert molecule
    assert molecule.inputFile == msg
    assert len(molecule.molTopol.atoms) == 29
    janitor.append(molecule.absHomeDir)
    janitor.append(molecule.tmpDir)
    os.remove(molecule.absInputFile)


@pytest.mark.parametrize(
    ("ct", "ft", "msg"),
    [
        (
            "bcc",
            "pdb",
            "-dr no -i 'benzene.mol2' -fi mol2 -o 'benzene_bcc_gaff2.mol2' -fo mol2 -c bcc -nc 0 -m 1 -s 2 -df 2 -at gaff2",
        ),
        (
            "bcc",
            "mol",
            "-dr no -i 'benzene.mol' -fi mdl -o 'benzene_bcc_gaff2.mol2' -fo mol2 -c bcc -nc 0 -m 1 -s 2 -df 2 -at gaff2",
        ),
        (
            "bcc",
            "mdl",
            "-dr no -i 'benzene.mdl' -fi mdl -o 'benzene_bcc_gaff2.mol2' -fo mol2 -c bcc -nc 0 -m 1 -s 2 -df 2 -at gaff2",
        ),
        ("user", "pdb", "cannot read charges from a PDB file"),
    ],
)
def test_sqm_tleap(janitor, capsys, ct, ft, msg):
    # check chargeType user with PDB -> use bcc
    # .mol and .mdl are the same file type
    molecule = ACTopol(f"benzene.{ft}", chargeType=ct, debug=True)
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
    janitor.append(molecule.absHomeDir)
    janitor.append(molecule.tmpDir)


def test_ekFlag(janitor):
    molecule = ACTopol(
        "benzene.pdb",
        ekFlag='''"qm_theory='AM1', grms_tol=0.0005, scfconv=1.d-10, ndiis_attempts=700, qmcharge=0"''',
        gmx4=True,
    )
    molecule.createACTopol()
    molecule.createMolTopol()
    assert molecule
    assert len(molecule.molTopol.atoms) == 12
    janitor.append(molecule.absHomeDir)
    janitor.append(molecule.tmpDir)


def test_time_limit(janitor):
    molecule = ACTopol("KKK.pdb", chargeType="bcc", debug=True, timeTol=2)
    with pytest.raises(Exception) as e_info:
        molecule.createACTopol()
    assert e_info.value.args[0] == "Semi-QM taking too long to finish... aborting!"
    janitor.append(molecule.absHomeDir)
    janitor.append(molecule.tmpDir)


def test_charge_user(janitor):
    molecule = ACTopol("ADPMg.mol2", chargeType="user", debug=True)
    molecule.createACTopol()
    molecule.createMolTopol()
    assert molecule
    assert len(molecule.molTopol.atoms) == 39
    assert len(molecule.molTopol.properDihedrals) == 125
    assert len(molecule.molTopol.improperDihedrals) == 7
    assert molecule.molTopol.atoms[0].charge == 0.1667
    assert molecule.molTopol.atoms[15].charge == -0.517
    janitor.append(molecule.absHomeDir)
    janitor.append(molecule.tmpDir)


@pytest.mark.parametrize(
    ("argv", "msg"),
    [
        (["-i", "cccc"], "Total time of execution:"),
        (["-x", "Base.inpcrd", "-p", "Base.prmtop", "-b", "vir_temp"], "Total time of execution:"),
        (["-wi", "cccc", "-b", "vir_temp"], ""),
        (
            ["-i", "wrong_res_set.pdb", "-b", "vir_temp"],
            "In vir_temp_AC.lib, residue name will be 'RSET' instead of 'SET' elsewhere",
        ),
        (
            ["-i", "wrong_res_num.pdb", "-b", "vir_temp"],
            "In vir_temp_AC.lib, residue name will be 'R100' instead of '100' elsewhere",
        ),
        (
            ["-i", "wrong_res_sym.pdb", "-b", "vir_temp"],
            "In vir_temp_AC.lib, residue name will be 'MOL' instead of '+++' elsewhere",
        ),
        (
            ["-i", "lower_res.pdb", "-b", "vir_temp"],
            "WARNING: this may raise problem with some applications like CNS",
        ),
        (
            ["-fi", "too_close.pdb", "-b", "vir_temp"],
            "You chose to proceed anyway with '-f' option. GOOD LUCK!",
        ),
        (
            ["-i", "no_res.pdb", "-b", "vir_temp"],
            "No residue name identified, using default resname: 'LIG'",
        ),
        (
            ["-i", "drift.mol2", "-c", "user", "-b", "vir_temp"],
            "Net charge drift '0.02020' bigger than tolerance '0.01000'",
        ),
    ],
)
def test_inputs(janitor, capsys, argv, msg):
    temp_base = "vir_temp"
    init_main(argv=argv)
    captured = capsys.readouterr()
    assert msg in captured.out
    _getoutput(f"rm -vfr {temp_base}* .*{temp_base}* smiles_molecule.acpype")


@pytest.mark.parametrize(
    ("argv", "code", "msg"),
    [
        (None, 2, " error: "),  # NOTE: None -> sys.argv from pytest
        (["-v"], 0, version),
        ([], 2, "error: missing input files"),
        (["-d", "-w"], 2, "error: argument -w/--verboseless: not allowed with argument -d/--debug"),
        (["-di", "AAAx.mol2"], 19, "ACPYPE FAILED: Input file AAAx.mol2 DOES NOT EXIST"),
        (["-zx", "ILDN.inpcrd", "-p", "ILDN.prmtop"], 19, "Likely trying to convert ILDN to RB"),
        (["-x", "glycam_exe.inpcrd", "-p", "glycam_corrupt.prmtop"], 19, "Skipping non-existent attributes dihedral_p"),
        (["-x", "glycam_exe.inpcrd", "-p", "glycam_empty.prmtop"], 19, "ERROR: ACPYPE FAILED: PRMTOP file empty?"),
        (["-x", "glycam_empty.inpcrd", "-p", "glycam_exe.prmtop"], 19, "ERROR: ACPYPE FAILED: INPCRD file empty?"),
        (["-di", "cccccc", "-n", "-1", "-b", "vir_temp"], 19, "Fatal Error!"),
        (["-di", " 123", "-b", "vir_temp"], 19, "ACPYPE FAILED: [Errno 2] No such file or directory"),
        (["-i", "double_res.pdb", "-b", "vir_temp"], 19, "Only ONE Residue is allowed for ACPYPE to work"),
        (["-i", "same_coord.pdb", "-b", "vir_temp"], 19, "Atoms with same coordinates in"),
        (["-i", "too_close.pdb", "-b", "vir_temp"], 19, "Atoms TOO close (<"),
        (["-i", "too_far.pdb", "-b", "vir_temp"], 19, "Atoms TOO scattered (>"),
        (["-di", " 123", "-x", "abc"], 2, "either '-i' or ('-p', '-x'), but not both"),
        (["-di", " 123", "-u"], 2, "option -u is only meaningful in 'amb2gmx' mode (args '-p' and '-x')"),
        (
            ["-di", "HEM.pdb", "-b", "vir_temp"],
            19,
            "No Gasteiger parameter for atom (ID: 42, Name: FE, Type: DU)",
            # Only Allowed: C, N, O, S, P, H, F, Cl, Br and I
        ),
    ],
)
def test_args_wrong_inputs(janitor, capsys, argv, code, msg):
    with pytest.raises(SystemExit) as e_info:
        init_main(argv=argv)
    captured = capsys.readouterr()
    assert msg in captured.err + captured.out
    assert e_info.typename == "SystemExit"
    assert e_info.value.code == code
    _getoutput("rm -vfr vir_temp* .*vir_temp*")
