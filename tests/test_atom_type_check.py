"""Tests for the early check that antechamber's atom types exist in the force field.

A nitro group through ``-a amber`` used to reach tleap with ``NO`` and ``DU`` atom
types nothing defines, after parmchk2 had silently written an empty frcmod that passed
as "Parmchk OK". The run then ended in a traceback. It now stops right after
antechamber, naming the atoms.
"""

from pathlib import Path

import pytest

from acpype.topol import ACTopol
from acpype.utils import readParmAtomTypes, unknownMol2Types

PARM_DAT = """PARM99 + frcmod.ff99SB, a synthetic excerpt
C  12.01         0.616  !            sp2 C carbonyl group
CA 12.01         0.360               sp2 C pure aromatic
N2 14.01         0.530               sp2 N in amino groups

C   H   HO  N   NA  NB  NC  N2  NT  N2  N3  N*  O   OH  OS  P   O2
CT-CT  310.0    1.526
"""

FRCMOD = """Remark line goes here
MASS
CO 12.01         0.616
2C 12.01         0.878

BOND
CO-2C  317.0    1.522
"""

MOL2 = """@<TRIPOS>MOLECULE
nitro
    4     0     1     0     0
SMALL
gas

@<TRIPOS>ATOM
      1 C1          0.0000     0.0000     0.0000 CA        1 UNL      0.000000
      2 N1          1.0000     0.0000     0.0000 NO        1 UNL      0.000000
      3 O1          2.0000     0.0000     0.0000 DU        1 UNL      0.000000
      4 C2          3.0000     0.0000     0.0000 CO        1 UNL      0.000000
@<TRIPOS>BOND
     1     1     2 1
"""

NITRO_SMILES = "CC(=O)Nc1ccc(cc1)[N+](=O)[O-]"


@pytest.fixture
def parm_files(tmp_path: Path) -> list[Path]:
    """Write one ``.dat`` and one frcmod and return their paths."""
    dat = tmp_path / "parm.dat"
    dat.write_text(PARM_DAT)
    frc = tmp_path / "frcmod.x"
    frc.write_text(FRCMOD)
    return [dat, frc]


def test_read_parm_atom_types_handles_both_layouts(parm_files: list[Path]) -> None:
    """The MASS block is found after the title of a .dat and after the MASS header of a frcmod."""
    types = readParmAtomTypes(parm_files)

    assert types == {"C", "CA", "N2", "CO", "2C"}
    assert "CT-CT" not in types, "reading must stop at the first blank line"


def test_unknown_mol2_types_flags_du_and_undefined_types() -> None:
    """DU and any type outside the known set are reported; defined types are not."""
    offending = unknownMol2Types(MOL2.splitlines(keepends=True), {"CA", "CO"})

    assert offending == [("N1", "NO"), ("O1", "DU")]


def test_unknown_mol2_types_always_flags_du() -> None:
    """A dummy type is an error even if a parameter file happened to list it."""
    assert unknownMol2Types(MOL2.splitlines(keepends=True), {"CA", "CO", "NO", "DU"}) == [("O1", "DU")]


def test_nitro_through_amber_fails_early_with_a_clear_message(
    janitor: list[str], capsys: pytest.CaptureFixture[str]
) -> None:
    """The nitro compound stops after antechamber, names the untypable atoms and points at gaff2."""
    molecule = ACTopol(NITRO_SMILES, basename="nitro", chargeType="gas", atomType="amber", debug=True, verbose=True)
    janitor.append(molecule.absHomeDir)

    assert molecule.createACTopol() is True

    out = capsys.readouterr().out
    assert "No parameters for atom(s)" in out
    assert "(NO)" in out
    assert "(DU)" in out
    assert "-a gaff2" in out
    assert "Parmchk OK" not in out, "parmchk2 must not run, let alone report success"
    assert "Could not find" not in out, "tleap must not be reached"


def test_nitro_through_gaff2_still_works(janitor: list[str]) -> None:
    """The same molecule is fine with GAFF2, which is what the message recommends."""
    molecule = ACTopol(NITRO_SMILES, basename="nitro2", chargeType="gas", atomType="gaff2", debug=True, verbose=False)
    janitor.append(molecule.absHomeDir)

    assert not molecule.createACTopol()
    assert (Path(molecule.absHomeDir) / molecule.acTopFileName).is_file()
