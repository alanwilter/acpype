"""Tests for the ``[ molecules ]`` table when a system carries ions.

The table is positional: it names molecules in the order the coordinates give them. ACPYPE
built it by counting each species instead, which is only right while every species sits in
one contiguous block. ``addions m Na+ 8 Cl- 8`` interleaves them, so the table claimed
eight sodiums followed by eight chlorides where the coordinates alternate, and every ion
was handed the other species' charge and mass.

Two ``addions`` calls, one species each, do give contiguous blocks. That is how every
fixture in this directory was built, and why the fault went unseen for so long.
"""

import shutil
import subprocess
from pathlib import Path

import pytest

from acpype.topol import MolTopol
from acpype.utils import solventTailRuns

GMX = shutil.which("gmx") or ""
KNOWN = ["Na+", "Cl-", "K+", "WAT"]

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
"""

INTERLEAVED = """source leaprc.protein.ff14SB
source leaprc.water.tip3p
m = sequence {{ ACE ALA ALA ALA NME }}
solvatebox m TIP3PBOX 6.0
addions m Na+ 4 Cl- 4
saveamberparm m {0}.prmtop {0}.inpcrd
quit
"""

BLOCKED = """source leaprc.protein.ff14SB
source leaprc.water.tip3p
m = sequence {{ ACE ALA ALA ALA NME }}
solvatebox m TIP3PBOX 6.0
addions m Na+ 4
addions m Cl- 4
saveamberparm m {0}.prmtop {0}.inpcrd
quit
"""


def tleap(script: str, prefix: str) -> tuple[str, str]:
    """Run the bundled tleap and return the prmtop and inpcrd it wrote."""
    exe = shutil.which("tleap")
    if exe is None:
        pytest.skip("no tleap available to build the fixture")
    Path(f"{prefix}.leap.in").write_text(script.format(prefix))
    subprocess.run([exe, "-f", f"{prefix}.leap.in"], capture_output=True, text=True, check=False)
    assert Path(f"{prefix}.prmtop").is_file(), f"tleap did not write {prefix}.prmtop"
    return f"{prefix}.prmtop", f"{prefix}.inpcrd"


def moleculesTable(path: Path) -> list[tuple[str, int]]:
    """Return the [ molecules ] entries of a written topology, in order."""
    lines = path.read_text().splitlines()
    start = lines.index("[ molecules ]")
    entries = []
    for line in lines[start + 1 :]:
        if line.startswith(";") or not line.strip():
            continue
        name, count = line.split()
        entries.append((name, int(count)))
    return entries


def convert(top: str, crd: str, janitor: list[str]) -> Path:
    """Convert a prmtop/inpcrd pair and return the directory ACPYPE wrote into."""
    molecule = MolTopol(acFileXyz=crd, acFileTop=top, amb2gmx=True, verbose=False)
    janitor.append(molecule.absHomeDir)
    molecule.writeGromacsTopolFiles()
    return Path(molecule.absHomeDir)


def test_interleaved_ions_become_one_entry_each() -> None:
    """Alternating ions yield a run each, in file order, rather than a tally per species."""
    labels = ["ALA", "ALA"] + ["Na+", "Cl-"] * 3 + ["WAT"] * 9

    assert solventTailRuns(labels, KNOWN) == [("Na+", 1), ("Cl-", 1)] * 3 + [("WAT", 9)]


def test_contiguous_ions_stay_grouped() -> None:
    """Blocked ions collapse back to one entry per species, as the old table always wrote."""
    labels = ["ALA", "ALA"] + ["Na+"] * 3 + ["Cl-"] * 4 + ["WAT"] * 9

    assert solventTailRuns(labels, KNOWN) == [("Na+", 3), ("Cl-", 4), ("WAT", 9)]


def test_ions_scattered_after_the_water() -> None:
    """addIonsRand drops ions in among the waters; the runs follow them there."""
    labels = ["ALA"] + ["WAT"] * 4 + ["Na+"] + ["WAT"] * 2 + ["Na+"] + ["WAT"] * 5

    assert solventTailRuns(labels, KNOWN) == [("WAT", 4), ("Na+", 1), ("WAT", 2), ("Na+", 1), ("WAT", 5)]


def test_solute_is_skipped_however_it_is_named() -> None:
    """Counting starts at the first solvent residue, whatever precedes it."""
    labels = ["MOL", "LIG", "ZN2", "WAT", "WAT"]

    assert solventTailRuns(labels, KNOWN) == [("WAT", 2)]


def test_a_system_without_solvent_has_no_entries() -> None:
    """A bare solute contributes nothing to the table rather than raising."""
    assert solventTailRuns(["MOL", "LIG"], KNOWN) == []
    assert solventTailRuns([], KNOWN) == []


@pytest.mark.parametrize(
    ("name", "script", "expected"),
    [
        ("inter", INTERLEAVED, [("Na+", 1), ("Cl-", 1)] * 4),
        ("block", BLOCKED, [("Na+", 4), ("Cl-", 4)]),
    ],
)
def test_written_table_follows_the_prmtop(
    name: str, script: str, expected: list[tuple[str, int]], janitor: list[str]
) -> None:
    """The emitted table reproduces the residue order tleap stored, either way round."""
    home = convert(*tleap(script, name), janitor)
    entries = moleculesTable(home / f"{name}_GMX.top")

    assert entries[0] == (name, 1), "the solute comes first"
    assert [(n.upper(), c) for n, c in expected] == entries[1 : 1 + len(expected)]
    assert entries[-1][0] == "WAT"


@pytest.mark.skipif(not GMX, reason="needs a GROMACS install")
@pytest.mark.parametrize(("name", "script"), [("inter", INTERLEAVED), ("block", BLOCKED)])
def test_grompp_accepts_the_result(name: str, script: str, janitor: list[str]) -> None:
    """grompp reads topology and coordinates back without a single name mismatch.

    This is the one of the two ion faults that GROMACS does catch. It catches it as a
    warning, though, and offers to carry on using the topology's names over the
    coordinates' -- which is precisely the wrong half to keep.
    """
    home = convert(*tleap(script, name), janitor)
    (home / "spe.mdp").write_text(SPE_MDP)

    grompp = subprocess.run(
        [GMX, "grompp", "-c", f"{name}_GMX.gro", "-p", f"{name}_GMX.top", "-f", "spe.mdp", "-o", "spe.tpr"],
        cwd=home,
        capture_output=True,
        text=True,
        check=False,
    )

    assert "does not match" not in grompp.stdout + grompp.stderr, "topology and coordinates disagree"
    assert grompp.returncode == 0, (grompp.stdout + grompp.stderr)[-2000:]
