"""Tests for reading typed records out of an AMBER prmtop.

Each block is introduced by a ``%FLAG`` and a ``%FORMAT`` record, and the format states
the type: ``E``, ``F`` and ``D`` are Fortran reals, ``I`` an integer, ``a`` text. The
type used to be guessed from the data instead, by looking for a ``.``, which mistook any
text field holding one for a line of floats.
"""

from pathlib import Path

import pytest

from acpype.topol import AbstractTopol

# An atom named N.4 is the case the old heuristic could not survive: a text record
# carrying a dot. Sybyl writes such names, and nothing stops them reaching a prmtop.
PRMTOP = """%VERSION  VERSION_STAMP = V0001.000  DATE = 01/01/26
%FLAG TITLE
%FORMAT(20a4)
made.up
%FLAG ATOM_NAME
%FORMAT(20a4)
C1  H1  N.4 O2
%FLAG AMBER_ATOM_TYPE
%FORMAT(20a4)
c3  hc  n4  o
%FLAG RESIDUE_LABEL
%FORMAT(20a4)
MOL
%FLAG RESIDUE_POINTER
%FORMAT(10I8)
       1
%FLAG CHARGE
%FORMAT(5E16.8)
  1.00000000E+00 -2.50000000E-01
%FLAG BOND_EQUIL_VALUE
%FORMAT(5E16.8)
  1.09690000E+00
%FLAG CMAP_PARAMETER_01
%FORMAT(8F9.5)
 -0.40490 -0.91563
"""


class Reader(AbstractTopol):
    """A bare topology that only exposes prmtop record reading."""

    def __init__(self, path: Path) -> None:
        """Load the prmtop, dropping %COMMENT records as MolTopol does."""
        self.topFileData = [
            line for line in path.read_text().splitlines(keepends=True) if not line.startswith("%COMMENT")
        ]
        self.level = 20


@pytest.fixture
def reader(tmp_path: Path) -> Reader:
    """A reader over the synthetic prmtop above."""
    path = tmp_path / "tiny.prmtop"
    path.write_text(PRMTOP)
    return Reader(path)


def test_text_record_holding_a_dot_stays_text(reader: Reader) -> None:
    """An atom named N.4 is read as text, where guessing by dot raised ValueError."""
    assert reader.getFlagData("ATOM_NAME") == ["C1", "H1", "N.4", "O2"]


@pytest.mark.parametrize(
    ("flag", "expected"),
    [
        ("CHARGE", [1.0, -0.25]),
        ("BOND_EQUIL_VALUE", [1.0969]),
        ("CMAP_PARAMETER_01", [-0.4049, -0.91563]),
    ],
)
def test_real_records_are_floats(reader: Reader, flag: str, expected: list[float]) -> None:
    """E and F formats are read as floats, whatever the data happens to look like."""
    assert reader.getFlagData(flag) == pytest.approx(expected)


def test_integer_record_is_integers(reader: Reader) -> None:
    """An I format is read as integers, not floats or text."""
    values = reader.getFlagData("RESIDUE_POINTER")

    assert values == [1]
    assert all(isinstance(v, int) for v in values)


def test_text_records_are_text(reader: Reader) -> None:
    """An a format is read as text, residue labels included."""
    assert reader.getFlagData("AMBER_ATOM_TYPE") == ["c3", "hc", "n4", "o"]
    assert reader.getFlagData("RESIDUE_LABEL") == ["MOL"]


def test_absent_flag_is_empty(reader: Reader) -> None:
    """A flag the file does not carry yields nothing rather than raising."""
    assert reader.getFlagData("CMAP_COUNT") == []


@pytest.mark.parametrize(
    ("flag", "kind", "head"),
    [
        ("CHARGE", float, [-3.6809046, 5.6853576]),
        ("MASS", float, [14.01, 1.008]),
        ("ATOM_NAME", str, ["N", "H2"]),
        ("AMBER_ATOM_TYPE", str, ["N3", "H"]),
        ("RESIDUE_LABEL", str, ["PRO", "GLN"]),
        ("RESIDUE_POINTER", int, [1, 17]),
        ("ATOMS_PER_MOLECULE", int, [1560, 1560]),
    ],
)
def test_a_real_prmtop_reads_unchanged(flag: str, kind: type, head: list[float | int | str]) -> None:
    """Every record type of a real prmtop keeps the type and values it always had."""
    values = Reader(Path("tests") / "ILDN.prmtop").getFlagData(flag)

    assert all(isinstance(v, kind) for v in values[:2])
    if kind is float:
        assert values[:2] == pytest.approx(head)
    else:
        assert values[:2] == head
