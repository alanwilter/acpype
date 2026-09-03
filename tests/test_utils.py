"""Tests for the parameter file helpers in :mod:`acpype.utils`."""

import re
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


def test_parse_frcmod_ignores_a_trailing_cmap_block(tmp_path: Path) -> None:
    """ff19SB's CMAP block after NONBON is skipped instead of polluting the NONBON section."""
    path = tmp_path / "frcmod.cmap"
    path.write_text(
        "remark\nMASS\nXC 12.01 0.360\n\nNONBON\n  XC  1.9080  0.1094\n\nCMAP\n"
        "%FLAG CMAP_COUNT   1\n%FLAG CMAP_RESLIST  1\nGLY\n%FLAG CMAP_PARAMETER\n"
        "  0.82366   1.09817   1.13106\n  1.38658   2.11377   3.63194\n"
    )

    parsed = parseFrcmod(path.read_text().splitlines(keepends=True))

    assert sorted(parsed) == ["MASS", "NONB"]
    assert list(parsed["NONB"]) == ["XC"]


def test_parse_frcmod_on_the_bundled_ff19SB() -> None:
    """The bundled frcmod.ff19SB parses with only real parameter keys."""
    for platform in ("amber_macos", "amber_linux"):
        path = Path("src/acpype") / platform / "dat/leap/parm/frcmod.ff19SB"
        if not path.is_file():
            continue
        parsed = parseFrcmod(path.read_text().splitlines(keepends=True))
        # 2C and 3C are real types; junk is a %FLAG record or a bare grid number.
        junk = [
            key
            for key in parsed.get("NONB", {})
            if key.startswith(("%", "CMAP")) or re.fullmatch(r"-?\d+(\.\d+)?", key)
        ]
        assert not junk, junk[:5]
        return
    pytest.skip("no bundled AmberTools to read frcmod.ff19SB from")
