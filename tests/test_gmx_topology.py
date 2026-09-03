"""Characterisation tests for the emitted GROMACS topology.

The rest of the suite asserts on in-memory objects -- atom counts, dihedral counts,
charges -- and never reads the ``.top`` file that is actually written. That leaves the
topology text itself unpinned, so a change to ``writeGromacsTop`` can rearrange or
corrupt it without turning a single test red.

These tests parse the generated file and pin its structure: which ``[ moleculetype ]``
blocks appear, how many atoms each holds, and what the ``[ molecules ]`` table says.
"""

import re
from dataclasses import dataclass, field
from pathlib import Path

import pytest

from acpype.topol import MolTopol

# A data line is anything that is not blank, not a ';' comment and not a preprocessor
# directive -- '#ifdef FLEXIBLE' and friends appear inside the water block.
_SKIP = re.compile(r"^\s*(?:[;#].*)?$")


@dataclass
class Moleculetype:
    """One ``[ moleculetype ]`` block of a GROMACS topology."""

    name: str
    sections: dict[str, list[str]] = field(default_factory=dict)

    def count(self, section: str) -> int:
        """Return how many data lines the named section holds."""
        return len(self.sections.get(section, []))


def parse_top(path: Path) -> tuple[list[Moleculetype], list[tuple[str, int]]]:
    """Split a GROMACS topology into its moleculetypes and its ``[ molecules ]`` table.

    Args:
        path: the ``.top`` file to read.

    Returns:
        The moleculetype blocks in file order, and the compound/count pairs of the
        ``[ molecules ]`` table.
    """
    blocks: list[Moleculetype] = []
    molecules: list[tuple[str, int]] = []
    current: Moleculetype | None = None
    section = ""
    in_molecules = False

    for raw in path.read_text().splitlines():
        header = re.match(r"^\s*\[\s*(.+?)\s*\]", raw)
        if header:
            section = header.group(1)
            in_molecules = section == "molecules"
            if section == "moleculetype":
                current = None  # the name is on the next data line
            elif current is not None:
                current.sections.setdefault(section, [])
            continue

        if _SKIP.match(raw):
            continue

        if in_molecules:
            name, count = raw.split()[:2]
            molecules.append((name, int(count)))
        elif section == "moleculetype" and current is None:
            current = Moleculetype(name=raw.split()[0])
            blocks.append(current)
        elif current is not None:
            current.sections.setdefault(section, []).append(raw)

    return blocks, molecules


def write_topology(base: str) -> tuple[MolTopol, Path]:
    """Convert an AMBER prmtop/inpcrd pair and return the molecule and its ``.top``."""
    molecule = MolTopol(acFileTop=f"{base}.prmtop", acFileXyz=f"{base}.inpcrd")
    molecule.writeGromacsTopolFiles()
    return molecule, Path(molecule.absHomeDir) / f"{base}_GMX.top"


# Current output, captured before any moleculetype splitting exists. Each entry is
# (basename, [(moleculetype name, atom count)], [(compound, nmols)]).
REFERENCE: list[tuple[str, list[tuple[str, int]], list[tuple[str, int]]]] = [
    (
        # Protein homodimer (1560 + 1560) plus an 80-atom DMP ligand, all merged.
        "ILDN",
        [("ILDN", 3200), ("NA+", 1), ("CL-", 1), ("WAT", 3)],
        [("ILDN", 1), ("NA+", 23), ("CL-", 27), ("WAT", 7227)],
    ),
    (
        # Protein (3313) + two ADP (39 each) + a Zn2+ that ionsDict does not know
        # about, so it lands in the solute block rather than getting its own.
        "ComplexG1",
        [("ComplexG1", 3392), ("NA+", 1), ("WAT", 3)],
        [("ComplexG1", 1), ("NA+", 15), ("WAT", 8855)],
    ),
    (
        # A single-molecule solute: nothing for a split to do.
        "RAMP1_ion",
        [("RAMP1_ion", 1289), ("NA+", 1), ("CL-", 1), ("WAT", 3)],
        [("RAMP1_ion", 1), ("NA+", 21), ("CL-", 19), ("WAT", 5763)],
    ),
]


@pytest.mark.parametrize(("base", "expected_blocks", "expected_molecules"), REFERENCE, ids=[r[0] for r in REFERENCE])
def test_moleculetype_blocks(
    janitor: list[str],
    base: str,
    expected_blocks: list[tuple[str, int]],
    expected_molecules: list[tuple[str, int]],
) -> None:
    """The emitted topology has the expected moleculetypes, sizes and molecules table."""
    molecule, top = write_topology(base)
    janitor.append(molecule.absHomeDir)

    blocks, molecules = parse_top(top)

    assert [(b.name, b.count("atoms")) for b in blocks] == expected_blocks
    assert molecules == expected_molecules


@pytest.mark.parametrize("base", [r[0] for r in REFERENCE])
def test_solute_block_merges_every_amber_molecule(janitor: list[str], base: str) -> None:
    """Every non-solvent AMBER molecule is currently merged into one moleculetype."""
    molecule, top = write_topology(base)
    janitor.append(molecule.absHomeDir)

    blocks, molecules = parse_top(top)
    solute = blocks[0]

    # tleap records the molecule decomposition it derived from bonded connectivity.
    # Atoms are contiguous per molecule, which is what a split would key off.
    per_molecule = molecule.getFlagData("ATOMS_PER_MOLECULE")
    solvent_atoms = sum(
        count * next(b.count("atoms") for b in blocks if b.name == name) for name, count in molecules[1:]
    )

    assert solute.count("atoms") == sum(per_molecule) - solvent_atoms


@pytest.mark.parametrize("base", [r[0] for r in REFERENCE])
def test_atom_indices_are_local_to_their_block(janitor: list[str], base: str) -> None:
    """Atom indices restart at 1 in every moleculetype and run without gaps."""
    molecule, top = write_topology(base)
    janitor.append(molecule.absHomeDir)

    blocks, _ = parse_top(top)

    for block in blocks:
        ids = [int(line.split()[0]) for line in block.sections["atoms"]]
        assert ids == list(range(1, len(ids) + 1)), f"{block.name} atom ids are not 1..N"


@pytest.mark.parametrize("base", [r[0] for r in REFERENCE])
def test_bonded_terms_stay_inside_their_block(janitor: list[str], base: str) -> None:
    """No bonded term references an atom outside the moleculetype that declares it."""
    molecule, top = write_topology(base)
    janitor.append(molecule.absHomeDir)

    blocks, _ = parse_top(top)

    # Guard against a parser that silently finds nothing, which would make the
    # containment assertion below pass for every topology ever written.
    solute = blocks[0]
    for section in ("bonds", "pairs", "angles", "dihedrals"):
        assert solute.count(section) > 0, f"parsed no {section} for {solute.name}"

    for block in blocks:
        natoms = block.count("atoms")
        for section, width in (("bonds", 2), ("pairs", 2), ("angles", 3), ("dihedrals", 4)):
            for line in block.sections.get(section, []):
                refs = [int(tok) for tok in line.split()[:width]]
                assert all(1 <= ref <= natoms for ref in refs), f"{block.name} {section}: {line.strip()!r}"
