"""Tests for the parameter file helpers in :mod:`acpype.utils`."""

from pathlib import Path

import pytest

from acpype.utils import parseFrcmod

# A frcmod whose MASS and IMPR sections carry no entries. frcmod.ff14SB, the only one
# acpype loads by default, has none of these, which is why the crash below stayed
# hidden; frcmod.ff99SB has several.
FRCMOD_WITH_EMPTY_SECTIONS = """remark goes here
MASS

BOND
CT-C   317.00   1.522

ANGL

DIHE
CT-C -N -CT   4   10.00   180.0   2.

IMPR

NONB
"""


@pytest.fixture
def frcmod_file(tmp_path: Path) -> Path:
    """Write a frcmod containing empty sections and return its path."""
    path = tmp_path / "frcmod.empty_sections"
    path.write_text(FRCMOD_WITH_EMPTY_SECTIONS)
    return path


def test_parse_frcmod_drops_empty_sections(frcmod_file: Path) -> None:
    """Sections with no entries are dropped rather than raising during iteration."""
    parsed = parseFrcmod(frcmod_file.read_text().splitlines(keepends=True))

    assert sorted(parsed) == ["BOND", "DIHE"]


def test_parse_frcmod_keeps_entries(frcmod_file: Path) -> None:
    """The entries of a populated section survive parsing, keyed by their atom types."""
    parsed = parseFrcmod(frcmod_file.read_text().splitlines(keepends=True))

    assert list(parsed["DIHE"]) == ["CT-C-N-CT"]
    assert parsed["DIHE"]["CT-C-N-CT"] == ["CT-C -N -CT   4   10.00   180.0   2."]


def test_parse_frcmod_handles_only_empty_sections(tmp_path: Path) -> None:
    """A frcmod with nothing but empty sections parses to an empty mapping."""
    path = tmp_path / "frcmod.blank"
    path.write_text("remark\nMASS\n\nBOND\n\nDIHE\n\n")

    assert parseFrcmod(path.read_text().splitlines(keepends=True)) == {}


def test_parse_frcmod_on_the_bundled_ff99SB() -> None:
    """The bundled frcmod.ff99SB parses; it is what exposed the iteration crash."""
    for platform in ("amber_macos", "amber_linux"):
        path = Path("src/acpype") / platform / "dat/leap/parm/frcmod.ff99SB"
        if not path.is_file():
            continue
        parsed = parseFrcmod(path.read_text().splitlines(keepends=True))
        assert "DIHE" in parsed
        assert all(parsed.values()), "an empty section survived parsing"
        return
    pytest.skip("no bundled AmberTools to read frcmod.ff99SB from")
