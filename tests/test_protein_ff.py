"""Tests for the protein force field handling behind ``-a amber``.

ff14SB keys its backbone torsions on residue-specific atom types (CX for the C-alpha,
CO, 2C, 3C, C8 in side chains) that antechamber never emits, and parmchk2 replaces
AMBER terms it already has with GAFF analogues when GAFF is merged into its parameter
file. Together those produced peptide topologies whose proper dihedrals disagreed with
GROMACS' own AMBER ports by 26-148%. These tests pin both halves of the fix.
"""

from pathlib import Path

import pytest

from acpype.params import FF14SB_RETYPES, PROTEIN_FF, PROTEIN_FF_LIBS
from acpype.topol import ACTopol
from acpype.utils import mergeFrcmodGaps, readAmberLibTypes, retypeMol2Atoms

LIB = """!!index array str
 "ALA"
!entry.ALA.unit.atoms table  str name  str type  int typex  int resx  int flags  int seq  int elmnt  dbl chg
 "N" "N" 0 1 131072 1 7 -0.415700
 "H" "H" 0 1 131072 2 1 0.271900
 "CA" "CX" 0 1 131072 3 6 0.033700
 "HA" "H1" 0 1 131072 4 1 0.082300
 "CB" "CT" 0 1 131072 5 6 -0.182500
 "C" "C" 0 1 131072 9 6 0.597300
 "O" "O" 0 1 131072 10 8 -0.567900
!entry.ALA.unit.atomspertinfo table  str pname  str ptype  int ptypex  int pelmnt  dbl pchg
 "N" "N" 0 -1 0.0
!entry.ASP.unit.atoms table  str name  str type  int typex  int resx  int flags  int seq  int elmnt  dbl chg
 "CA" "CX" 0 1 131072 3 6 0.0
 "CB" "2C" 0 1 131072 5 6 0.0
 "CG" "CO" 0 1 131072 8 6 0.0
!entry.NALA.unit.atoms table  str name  str type  int typex  int resx  int flags  int seq  int elmnt  dbl chg
 "N" "N3" 0 1 131072 1 7 0.0
 "CA" "CX" 0 1 131072 5 6 0.0
!entry.CALA.unit.atoms table  str name  str type  int typex  int resx  int flags  int seq  int elmnt  dbl chg
 "CA" "CX" 0 1 131072 3 6 0.0
 "OXT" "O2" 0 1 131072 11 8 0.0
"""

# An antechamber-typed tripeptide: CA is CT everywhere, the C-terminal oxygens carry
# the GROMACS names OC1/OC2, and residue numbering starts at 4 on purpose.
MOL2 = """@<TRIPOS>MOLECULE
AAA
    9     0     1     0     0
SMALL
gas

@<TRIPOS>ATOM
      1 N          -1.8940     0.4980     0.0640 N3        4 ALA     -0.415700
      2 CA         -0.5470     0.0730     0.4520 CT        4 ALA      0.033700
      3 CB         -0.1110    -1.1560    -0.3340 CT        4 ALA     -0.182500
      4 C           0.4630     1.2280     0.2580 C         4 ALA      0.597300
      5 CA          1.8990     2.2440    -1.0230 CT        5 ALA      0.033700
      6 CA          4.5430     3.5040    -0.6840 CT        6 ALA      0.033700
      7 HA          4.9040     4.1940    -1.4460 H1        6 ALA      0.082300
      8 OC1         5.7670     2.2560     0.7240 O2        6 ALA     -0.805500
      9 OC2         5.2300     4.2350     1.4840 O2        6 ALA     -0.805500
@<TRIPOS>BOND
"""

FRCMOD_AMBER_ONLY = """Remark line goes here
MASS

BOND

ANGLE

DIHE
CT-CT-N -XX   1    0.000         0.000           2.000      ATTN, need revision

IMPROPER

NONBON

"""

FRCMOD_WITH_GAFF = """Remark line goes here
MASS

BOND

ANGLE

DIHE
CT-C -N -CT   1    0.000         0.000          -2.000      same as c3-c -n -c3
CT-C -N -CT   1    1.500       180.000           1.000      same as c3-c -n -c3, penalty score=  0.0
CT-CT-N -XX   1    0.000         0.000          -2.000      same as c3-c3-n -xx
CT-CT-N -XX   1    0.700       180.000           3.000      same as c3-c3-n -xx, penalty score=  0.0

IMPROPER

NONBON

"""


@pytest.fixture
def templates(tmp_path: Path) -> dict[tuple[str, str], str]:
    """Parse the synthetic library above."""
    lib = tmp_path / "amino.lib"
    lib.write_text(LIB)
    return readAmberLibTypes([lib])


def test_read_amber_lib_types_reads_atoms_tables_only(templates: dict[tuple[str, str], str]) -> None:
    """Only the ``unit.atoms`` tables contribute, one type per residue and atom name."""
    assert templates[("ALA", "CA")] == "CX"
    assert templates[("ASP", "CG")] == "CO"
    assert templates[("NALA", "N")] == "N3"
    assert ("ALA", "pname") not in templates
    assert len(templates) == 14


def test_retype_applies_only_ff14sb_types(templates: dict[tuple[str, str], str]) -> None:
    """Every CA becomes CX; atoms whose template type is not an ff14SB type are untouched."""
    lines, changes = retypeMol2Atoms(MOL2.splitlines(keepends=True), templates, FF14SB_RETYPES)

    assert [(res, atom, old, new) for res, atom, old, new in changes] == [
        ("ALA", "CA", "CT", "CX"),
        ("ALA", "CA", "CT", "CX"),
        ("ALA", "CA", "CT", "CX"),
    ]
    types = {line.split()[1]: line.split()[5] for line in lines if line.startswith("      ")}
    assert types["CA"] == "CX"
    assert types["CB"] == "CT", "CB is CT in ff14SB ALA and must not change"
    assert types["N"] == "N3", "N3 matches the template but is not an ff14SB type"
    assert types["OC1"] == "O2", "unknown atom names are left alone"


def test_retype_preserves_column_layout(templates: dict[tuple[str, str], str]) -> None:
    """The retyped line differs from the original only in the type token."""
    original = MOL2.splitlines(keepends=True)
    lines, _ = retypeMol2Atoms(original, templates, FF14SB_RETYPES)

    before = next(line for line in original if len(line.split()) > 5 and line.split()[1] == "CA")
    after = next(line for line in lines if len(line.split()) > 5 and line.split()[1] == "CA")
    assert len(before) == len(after)
    assert before.replace(" CT ", " CX ", 1) == after


def test_retype_uses_terminal_templates_by_position(tmp_path: Path) -> None:
    """The first and last residues are matched against the N- and C-terminal templates."""
    lib = tmp_path / "term.lib"
    lib.write_text(
        "!entry.NALA.unit.atoms table  str name  str type  int typex  int resx  int flags  int seq  int elmnt  dbl chg\n"
        ' "CA" "3C" 0 1 131072 1 6 0.0\n'
        "!entry.CALA.unit.atoms table  str name  str type  int typex  int resx  int flags  int seq  int elmnt  dbl chg\n"
        ' "CA" "C8" 0 1 131072 1 6 0.0\n'
        "!entry.ALA.unit.atoms table  str name  str type  int typex  int resx  int flags  int seq  int elmnt  dbl chg\n"
        ' "CA" "CX" 0 1 131072 1 6 0.0\n'
    )
    _, changes = retypeMol2Atoms(MOL2.splitlines(keepends=True), readAmberLibTypes([lib]), FF14SB_RETYPES)

    assert [new for _, _, _, new in changes] == ["3C", "CX", "C8"]


def test_retype_without_atom_section_is_a_noop() -> None:
    """A file with no ``@<TRIPOS>ATOM`` block comes back unchanged."""
    lines = ["@<TRIPOS>MOLECULE\n", "x\n"]
    assert retypeMol2Atoms(lines, {("ALA", "CA"): "CX"}, FF14SB_RETYPES) == (lines, [])


def test_merge_frcmod_gaps_keeps_gaff_only_for_genuine_gaps() -> None:
    """GAFF values survive for keys AMBER lacks; GAFF overrides of AMBER terms are dropped."""
    merged = "".join(
        mergeFrcmodGaps(FRCMOD_AMBER_ONLY.splitlines(keepends=True), FRCMOD_WITH_GAFF.splitlines(keepends=True))
    )

    assert "CT-C -N -CT" not in merged, "the amide torsion is covered by AMBER's X-C-N-X"
    assert "CT-CT-N -XX   1    0.700       180.000           3.000" in merged, "the gap takes the GAFF value"
    assert "ATTN" not in merged
    for header in ("MASS", "BOND", "ANGLE", "DIHE", "IMPROPER", "NONBON"):
        assert f"\n{header}\n" in merged


def test_merge_frcmod_gaps_falls_back_to_amber_run_values() -> None:
    """A gap the GAFF run did not fill keeps the AMBER run's line, ATTN flag included."""
    merged = "".join(
        mergeFrcmodGaps(FRCMOD_AMBER_ONLY.splitlines(keepends=True), FRCMOD_AMBER_ONLY.splitlines(keepends=True))
    )

    assert "CT-CT-N -XX   1    0.000         0.000           2.000      ATTN, need revision" in merged


def test_protein_ff_table_is_consistent() -> None:
    """Every configured force field names a leaprc, a parm file and a frcmod."""
    for name, ff in PROTEIN_FF.items():
        assert set(ff) == {"leaprc", "parm", "frcmod"}, name
        assert ff["leaprc"].startswith(("leaprc", "oldff/leaprc")), name
        assert name in PROTEIN_FF_LIBS, name
    assert PROTEIN_FF_LIBS["ff14SB"], "ff14SB needs templates to retype from"
    assert not PROTEIN_FF_LIBS["ff99SB"], "ff99SB keys torsions on antechamber's own types"


def test_actopol_rejects_unknown_protein_ff(janitor: list[str]) -> None:
    """An unsupported protein force field name is refused up front."""
    with pytest.raises(ValueError, match="protein_ff must be one of"):
        ACTopol("AAA.pdb", chargeType="gas", atomType="amber", protein_ff="ff03", verbose=False)


@pytest.mark.parametrize(
    ("protein_ff", "ca_type"),
    [("ff14SB", "CX"), ("ff99SB", "CT")],
)
def test_amber_tripeptide_end_to_end(janitor: list[str], protein_ff: str, ca_type: str) -> None:
    """A real tripeptide through -a amber gets the right C-alpha type and no GAFF backbone terms."""
    molecule = ACTopol("AAA.pdb", chargeType="gas", atomType="amber", protein_ff=protein_ff, debug=True, verbose=False)
    janitor.append(molecule.absHomeDir)
    assert not molecule.createACTopol()

    home = Path(molecule.absHomeDir)
    mol2 = (home / molecule.acMol2FileName).read_text()
    atoms = [
        line.split() for line in mol2.split("@<TRIPOS>ATOM")[1].split("@<TRIPOS>BOND")[0].splitlines() if line.strip()
    ]
    assert {t[5] for t in atoms if t[1] == "CA"} == {ca_type}

    frcmod = (home / molecule.acFrcmodFileName).read_text()
    assert "same as c3-c -n -c3" not in frcmod, "omega must come from AMBER's X-C-N-X, not GAFF"
    assert "same as h1-c3-c -o" not in frcmod, "the carboxylate torsion is zero in AMBER, not GAFF's 0.8"
    assert "ATTN" not in frcmod, "a standard peptide needs no invented parameters"
