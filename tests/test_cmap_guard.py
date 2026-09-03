"""Tests for refusing prmtops whose CMAP terms the GROMACS writer cannot express yet.

ff19SB keeps its backbone correction in CMAP grids. ACPYPE's writer knows nothing of
them, so converting such a prmtop used to produce a topology that grompp accepted and
that silently lacked the correction. The fixtures are built with the bundled tleap so
the test exercises a real ff19SB prmtop rather than a hand-edited one.
"""

import shutil
import subprocess
from pathlib import Path

import pytest

from acpype.cli import init_main
from acpype.errors import AcpypeError, UnsupportedTopologyError
from acpype.topol import MolTopol

LEAP = "source {leaprc}\nm = sequence {{ NALA ALA CALA }}\nsaveamberparm m {prefix}.prmtop {prefix}.inpcrd\nquit\n"


def build_tripeptide(prefix: str, leaprc: str) -> tuple[str, str]:
    """Build an ALA tripeptide with the given leaprc and return its prmtop/inpcrd names."""
    tleap = shutil.which("tleap")
    if tleap is None:
        pytest.skip("no tleap available to build the fixture")
    Path(f"{prefix}.leap.in").write_text(LEAP.format(leaprc=leaprc, prefix=prefix))
    subprocess.run([tleap, "-f", f"{prefix}.leap.in"], capture_output=True, text=True, check=False)
    assert Path(f"{prefix}.prmtop").is_file(), f"tleap did not write {prefix}.prmtop"
    return f"{prefix}.prmtop", f"{prefix}.inpcrd"


@pytest.fixture
def ff19sb(janitor: list[str]) -> tuple[str, str]:
    """A tripeptide prmtop that carries CMAP terms."""
    top, crd = build_tripeptide("ala19", "leaprc.protein.ff19SB")
    assert "%FLAG CMAP_COUNT" in Path(top).read_text()
    return top, crd


@pytest.fixture
def ff14sb(janitor: list[str]) -> tuple[str, str]:
    """The same tripeptide without CMAP terms, as a control."""
    top, crd = build_tripeptide("ala14", "leaprc.protein.ff14SB")
    assert "CMAP" not in Path(top).read_text()
    return top, crd


def test_unsupported_topology_is_an_acpype_error() -> None:
    """The refusal is a deliberate one, so the CLI prints it without a traceback."""
    assert issubclass(UnsupportedTopologyError, AcpypeError)


def test_moltopol_refuses_a_cmap_prmtop(ff19sb: tuple[str, str]) -> None:
    """A prmtop with CMAP terms is refused before anything is written."""
    top, crd = ff19sb

    with pytest.raises(UnsupportedTopologyError, match=r"1 CMAP term\(s\) of 1 type\(s\), as ff19SB"):
        MolTopol(acFileXyz=crd, acFileTop=top, amb2gmx=True, verbose=False)

    assert not list(Path().glob("*/*_GMX.top")), "no GROMACS topology may be written"


def test_moltopol_converts_the_same_peptide_without_cmap(ff14sb: tuple[str, str], janitor: list[str]) -> None:
    """The guard keys on the CMAP flags, not on the force field: ff14SB still converts."""
    top, crd = ff14sb

    molecule = MolTopol(acFileXyz=crd, acFileTop=top, amb2gmx=True, verbose=False)
    janitor.append(molecule.absHomeDir)
    molecule.writeGromacsTopolFiles()

    assert (Path(molecule.absHomeDir) / "ala14_GMX.top").is_file()


def test_cli_prints_the_refusal_without_a_traceback(
    ff19sb: tuple[str, str], capsys: pytest.CaptureFixture[str]
) -> None:
    """The command line reports the refusal in one line and exits with ACPYPE's failure code."""
    top, crd = ff19sb

    with pytest.raises(SystemExit) as exc:
        init_main(argv=["-p", top, "-x", crd])

    assert exc.value.code == 19
    out = capsys.readouterr()
    text = out.out + out.err
    assert "ACPYPE FAILED" in text
    assert "CMAP" in text
    assert "ff14SB or ff99SB" in text
    assert "Traceback" not in text
