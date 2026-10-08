"""Tests for converting an AMBER periodic cell into GROMACS box vectors.

AMBER stores a box as three lengths and three angles; GROMACS wants three cell vectors
reduced so v1 lies along x and v2 in the xy plane. For a truncated octahedron the
reduction gives a negative v2x and v3y, which is the orientation AMBER writes its
coordinates in. ACPYPE used to emit GROMACS' own canonical octahedron instead, which
mirrors those two components: grompp accepted the result and every atom sat in the
wrong periodic image, with a Lennard-Jones energy ten orders of magnitude too large.
"""

import math
import shutil
import subprocess
from pathlib import Path

import pytest

from acpype.topol import MolTopol
from acpype.utils import cellToBoxVectors

GMX = shutil.which("gmx") or ""
OCTAHEDRON = 109.4712190

SPE_MDP = """integrator = md
nsteps = 0
cutoff-scheme = Verlet
nstlist = 10
coulombtype = PME
rcoulomb = 0.9
rvdw = 0.9
rlist = 0.9
pbc = xyz
continuation = yes
nstcalcenergy = 1
nstenergy = 1
"""


def determinant(v1: tuple[float, ...], v2: tuple[float, ...], v3: tuple[float, ...]) -> float:
    """Return the signed volume the three cell vectors enclose."""
    return (
        v1[0] * (v2[1] * v3[2] - v2[2] * v3[1])
        - v1[1] * (v2[0] * v3[2] - v2[2] * v3[0])
        + v1[2] * (v2[0] * v3[1] - v2[1] * v3[0])
    )


def cell_volume(a: float, b: float, c: float, alpha: float, beta: float, gamma: float) -> float:
    """Return the cell volume straight from the lengths and angles, as a cross-check."""
    ca, cb, cg = (math.cos(math.radians(x)) for x in (alpha, beta, gamma))
    return a * b * c * math.sqrt(1 - ca * ca - cb * cb - cg * cg + 2 * ca * cb * cg)


def tleap(script: str, prefix: str) -> tuple[str, str]:
    """Run the bundled tleap and return the prmtop and inpcrd it wrote."""
    exe = shutil.which("tleap")
    if exe is None:
        pytest.skip("no tleap available to build the fixture")
    Path(f"{prefix}.leap.in").write_text(script)
    subprocess.run([exe, "-f", f"{prefix}.leap.in"], capture_output=True, text=True, check=False)
    assert Path(f"{prefix}.prmtop").is_file(), f"tleap did not write {prefix}.prmtop"
    return f"{prefix}.prmtop", f"{prefix}.inpcrd"


@pytest.fixture
def octahedral_system(janitor: list[str]) -> tuple[str, str]:
    """A peptide solvated in a truncated octahedron, as `solvateoct` builds it."""
    script = (
        "source leaprc.protein.ff14SB\nsource leaprc.water.tip3p\n"
        "m = sequence { NALA ALA ALA ALA CALA }\nsolvateoct m TIP3PBOX 8.0\n"
        "saveamberparm m oct.prmtop oct.inpcrd\nquit\n"
    )
    return tleap(script, "oct")


def test_rectangular_cell_is_diagonal() -> None:
    """Right angles give a diagonal box, with no off-diagonal components to write."""
    v1, v2, v3 = cellToBoxVectors([3.79678, 3.33140, 2.61916], [90.0, 90.0, 90.0])

    assert (v1[0], v2[1], v3[2]) == pytest.approx((3.79678, 3.33140, 2.61916))
    assert max(abs(v2[0]), abs(v3[0]), abs(v3[1])) < 1e-6


def test_truncated_octahedron_matches_ambers_orientation() -> None:
    """v2x and v3y come out negative: the orientation AMBER's coordinates assume."""
    d = 3.345979
    v1, v2, v3 = cellToBoxVectors([d, d, d], [OCTAHEDRON] * 3)

    assert v1 == pytest.approx((d, 0.0, 0.0))
    assert v2 == pytest.approx((-d / 3, 2 * math.sqrt(2) * d / 3, 0.0))
    assert v3 == pytest.approx((-d / 3, -math.sqrt(2) * d / 3, math.sqrt(6) * d / 3))
    assert v2[0] < 0 and v3[1] < 0, "mirroring these is what broke the periodic images"


@pytest.mark.parametrize(
    ("lengths", "angles"),
    [
        ([3.79678, 3.33140, 2.61916], [90.0, 90.0, 90.0]),
        ([3.345979] * 3, [OCTAHEDRON] * 3),
        ([3.0] * 3, [60.0, 60.0, 90.0]),
        ([2.5, 3.0, 3.5], [104.0, 98.0, 110.0]),
    ],
)
def test_volume_is_preserved(lengths: list[float], angles: list[float]) -> None:
    """The reduced vectors enclose the cell's own volume, and are right-handed."""
    volume = determinant(*cellToBoxVectors(lengths, angles))

    assert volume > 0, "GROMACS requires a right-handed box"
    assert volume == pytest.approx(cell_volume(*lengths, *angles))


@pytest.mark.parametrize("alpha", [60.0, 104.0, 120.0])
def test_angles_the_old_code_could_not_handle(alpha: float) -> None:
    """Any angle reduces; only 90 and 109.47 used to, the rest raised UnboundLocalError."""
    v1, v2, v3 = cellToBoxVectors([3.0, 3.0, 3.0], [alpha, alpha, alpha])

    assert determinant(v1, v2, v3) > 0
    assert v1[1] == v1[2] == v2[2] == 0.0, "GROMACS needs v1 along x and v2 in the xy plane"


def test_gro_box_line_for_an_octahedron(octahedral_system: tuple[str, str], janitor: list[str]) -> None:
    """The written box carries nine fields with AMBER's signs on v2x, v3x and v3y."""
    top, crd = octahedral_system
    molecule = MolTopol(acFileXyz=crd, acFileTop=top, amb2gmx=True, verbose=False)
    janitor.append(molecule.absHomeDir)
    molecule.writeGromacsTopolFiles()

    fields = [float(x) for x in (Path(molecule.absHomeDir) / "oct_GMX.gro").read_text().splitlines()[-1].split()]
    assert len(fields) == 9
    v1x, v2y, v3z, v1y, v1z, v2x, v2z, v3x, v3y = fields
    assert (v1y, v1z, v2z) == (0.0, 0.0, 0.0)
    assert v2x < 0 and v3x < 0 and v3y < 0
    assert v2y == pytest.approx(2 * math.sqrt(2) * v1x / 3, rel=1e-4)
    assert v3z == pytest.approx(math.sqrt(6) * v1x / 3, rel=1e-4)


@pytest.mark.skipif(not GMX, reason="needs a GROMACS install")
def test_octahedral_system_is_physical(octahedral_system: tuple[str, str], janitor: list[str]) -> None:
    """GROMACS reads the converted octahedron back with a sane Lennard-Jones energy.

    The mirrored box passed grompp just as happily; it was the energy that gave it
    away, at 3.2e10 kJ/mol against the few thousand a solvated peptide should show.
    """
    top, crd = octahedral_system
    molecule = MolTopol(acFileXyz=crd, acFileTop=top, amb2gmx=True, verbose=False)
    janitor.append(molecule.absHomeDir)
    molecule.writeGromacsTopolFiles()
    home = Path(molecule.absHomeDir)
    (home / "spe.mdp").write_text(SPE_MDP)

    grompp = subprocess.run(
        [GMX, "grompp", "-c", "oct_GMX.gro", "-p", "oct_GMX.top", "-f", "spe.mdp", "-o", "spe.tpr", "-maxwarn", "0"],
        cwd=home,
        capture_output=True,
        text=True,
        check=False,
    )
    assert grompp.returncode == 0, grompp.stderr[-2000:]
    subprocess.run(
        [GMX, "mdrun", "-deffnm", "spe", "-s", "spe.tpr", "-nt", "2"], cwd=home, capture_output=True, check=False
    )
    result = subprocess.run(
        [GMX, "energy", "-f", "spe.edr", "-o", "spe.xvg"],
        cwd=home,
        input="LJ-(SR)\n\n",
        capture_output=True,
        text=True,
        check=False,
    )

    energies = [
        float(line.split()[-5]) for line in (result.stdout + result.stderr).splitlines() if line.endswith("(kJ/mol)")
    ]
    assert energies, "gmx energy reported no Lennard-Jones term"
    assert abs(energies[0]) < 1.0e5, f"LJ (SR) is {energies[0]:.3e} kJ/mol; the periodic images overlap"
