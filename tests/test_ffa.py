"""Unit tests for the GROMACS force field addon writer (``--ffa``, acpype.ffa).

The addon extends a base force field with one ligand residue for
``gmx pdb2gmx -ff <base> -ffa <addon>`` (GROMACS 2027+). These tests convert the
single-residue K11 ligand (GAFF parameters in ``K11.prmtop``/``K11.inpcrd``) and
pin what the addon directory must contain:

* the file layout pdb2gmx expects, including the ``#include "ffbonded.itp"``
  in ``forcefield-addon.itp`` through which the bonded type tables reach the
  generated topology (and hence grompp);
* atom type definitions, bonded type tables and the residue topology that are
  consistent with the parsed molecule: same types, same Fourier terms, same
  charges.

The molecule-level numbers are taken from the same in-memory ``MolTopol`` that
produced the files, so the tests pin the writer, not a hand-copied expectation.
"""

import re
from pathlib import Path

import pytest

from acpype.errors import UnsupportedTopologyError
from acpype.params import cal, radPi
from acpype.topol import MolTopol

# Data lines only: skip blanks, ';' comments and '#' preprocessor directives.
_SKIP = re.compile(r"^\s*(?:[;#].*)?$")


def parse_sections(path: Path) -> dict[str, list[str]]:
    """Split an itp-like file into ``{section name: [data line, ...]}``."""
    sections: dict[str, list[str]] = {}
    section = ""
    for raw in path.read_text().splitlines():
        header = re.match(r"^\s*\[\s*(.+?)\s*\]", raw)
        if header:
            section = header.group(1)
            sections.setdefault(section, [])
            continue
        if section and not _SKIP.match(raw):
            sections[section].append(raw)
    return sections


def canon(types: tuple[str, ...]) -> tuple[str, ...]:
    return min(types, types[::-1])


def write_addon(base: str) -> tuple[MolTopol, Path]:
    """Convert a prmtop/inpcrd pair with ``ffa=True`` and return molecule + addon dir."""
    molecule = MolTopol(acFileTop=f"{base}.prmtop", acFileXyz=f"{base}.inpcrd", ffa=True)
    molecule.writeGromacsTopolFiles()
    return molecule, Path(molecule.absHomeDir) / f"{molecule.baseName}.ffa"


def solute_atoms(mol: MolTopol) -> list:
    solvent = (mol.freeIonResidues() if mol.amb2gmx else set()) | {"WAT"}
    return [a for a in mol.atomsGromacs if mol.residueLabel[a.resid] not in solvent]


def test_addon_files_are_written(janitor: list[str]) -> None:
    """All files of the addon directory exist under <baseName>.ffa."""
    molecule, addon = write_addon("K11")
    janitor.append(molecule.absHomeDir)
    assert addon.is_dir()
    for name in (
        "compatibility.dat",
        "residuetypes.dat",
        "atomtypes.atp",
        "forcefield-addon.itp",
        "ffbonded.itp",
        "K11.rtp",
        "README.txt",
    ):
        assert (addon / name).is_file(), f"missing {name}"


def test_residue_and_compatibility(janitor: list[str]) -> None:
    """The residue becomes a molecule of its own and any base ff is allowed."""
    molecule, addon = write_addon("K11")
    janitor.append(molecule.absHomeDir)
    assert addon.joinpath("compatibility.dat").read_text().strip() == "*"
    assert addon.joinpath("residuetypes.dat").read_text().strip() == "K11 K11"


def test_bonded_types_reach_the_topology(janitor: list[str]) -> None:
    """forcefield-addon.itp includes ffbonded.itp -- the only path to grompp."""
    molecule, addon = write_addon("K11")
    janitor.append(molecule.absHomeDir)
    text = (addon / "forcefield-addon.itp").read_text()
    assert '#include "ffbonded.itp"' in text
    atomtypes = parse_sections(addon / "forcefield-addon.itp")["atomtypes"]
    # Full LJ definitions: name bondA mass charge ptype sigma epsilon.
    assert all(len(line.split()) == 7 for line in atomtypes)
    assert all(float(line.split()[5]) > 0.0 and float(line.split()[6]) > 0.0 for line in atomtypes)


def test_atom_types_match_the_molecule(janitor: list[str]) -> None:
    """atomtypes.atp lists exactly the solute atom types."""
    molecule, addon = write_addon("K11")
    janitor.append(molecule.absHomeDir)
    atp = {line.split()[0] for line in addon.joinpath("atomtypes.atp").read_text().splitlines() if line.strip()}
    assert atp == {a.atomType.atomTypeName for a in solute_atoms(molecule)}


def test_bond_and_angle_type_tables(janitor: list[str]) -> None:
    """Every canonical bond/angle type appears in both orientations with the model parameters."""
    molecule, addon = write_addon("K11")
    janitor.append(molecule.absHomeDir)
    bonded = parse_sections(addon / "ffbonded.itp")

    expected = {}
    for bond in molecule.bonds:
        types = tuple(a.atomType.atomTypeName for a in bond.atoms)
        expected.setdefault(canon(types), []).append((bond.rEq * 0.1, bond.kBond * 200 * cal))
    written = {}
    for line in bonded["bondtypes"]:
        t = line.split()
        assert int(t[2]) == 1
        written.setdefault(canon(tuple(t[:2])), []).append((float(t[3]), float(t[4])))
    assert set(written) == set(expected)

    expected = {}
    for angle in molecule.angles:
        types = tuple(a.atomType.atomTypeName for a in angle.atoms)
        expected.setdefault(canon(types), []).append((angle.thetaEq * radPi, 2 * cal * angle.kTheta))
    written = {}
    for line in bonded["angletypes"]:
        t = line.split()
        assert int(t[3]) == 1
        written.setdefault(canon(tuple(t[:3])), []).append((float(t[4]), float(t[5])))
    assert set(written) == set(expected)


def test_dihedral_fourier_terms(janitor: list[str]) -> None:
    """The [ dihedraltypes ] Fourier terms are exactly the AMBER torsion terms."""
    molecule, addon = write_addon("K11")
    janitor.append(molecule.absHomeDir)
    lines = parse_sections(addon / "ffbonded.itp")["dihedraltypes"]

    def terms_of(func: int) -> dict[tuple, set]:
        out: dict[tuple, set] = {}
        for line in lines:
            t = line.split()
            if int(t[4]) != func:
                continue
            phase, kd, pn = float(t[5]), float(t[6]), int(t[7])
            out.setdefault(canon(tuple(t[:4])), set()).add((phase, kd, pn))
        return out

    proper = terms_of(9)
    expected: dict[tuple, set] = {}
    for dih in molecule.properDihedralsGmx45:
        types = tuple(a.atomType.atomTypeName for a in dih[0])
        expected.setdefault(canon(types), set()).add((round(dih[1], 2), round(dih[2] * cal, 5), int(dih[3])))
    assert set(proper) == set(expected)
    for quartet, want in expected.items():
        got = proper[quartet]
        assert len(want) == len(got), f"{quartet}: {want} != {got}"
        for w, g in zip(sorted(want), sorted(got), strict=True):
            assert abs(w[0] - g[0]) < 6e-3 and abs(w[1] - g[1]) < 2e-5 and w[2] == g[2], f"{quartet}: {want} != {got}"

    improper = terms_of(4)
    expected = {}
    for dih in molecule.improperDihedrals:
        types = tuple(a.atomType.atomTypeName for a in dih.atoms)
        expected.setdefault(canon(types), set()).add(
            (round(dih.phase * radPi, 2), round(dih.kPhi * cal, 5), int(dih.period))
        )
    assert set(improper) == set(expected)


def test_rtp_atoms_charges_and_connectivity(janitor: list[str]) -> None:
    """The .rtp residue block carries the model atoms, charges and connectivity."""
    molecule, addon = write_addon("K11")
    janitor.append(molecule.absHomeDir)
    rtp = parse_sections(addon / "K11.rtp")

    atoms = rtp["atoms"]
    model = {a.atomName: a.charge for a in solute_atoms(molecule)}
    assert len(atoms) == len(model)
    for line in atoms:
        name, _type, charge = line.split()[:3]
        assert float(charge) == pytest.approx(model[name], abs=1e-6)
    assert sum(float(t.split()[2]) for t in atoms) == pytest.approx(sum(model.values()), abs=1e-4)

    assert len(rtp["bonds"]) == len({tuple(a.atomName for a in b.atoms) for b in molecule.bonds})
    assert len(rtp["angles"]) == len({tuple(a.atomName for a in ang.atoms) for ang in molecule.angles})
    quartets = {tuple(a.atomName for a in dih[0]) for dih in molecule.properDihedralsGmx45}
    assert len(rtp["dihedrals"]) == len(quartets)


def test_addon_rejects_multi_residue_solute(janitor: list[str]) -> None:
    """A multi-residue solute cannot be represented by a single-residue addon."""
    molecule = MolTopol(acFileTop="glycam_exe.prmtop", acFileXyz="glycam_exe.inpcrd", ffa=True)
    with pytest.raises(UnsupportedTopologyError, match="exactly one solute residue"):
        molecule.writeGromacsTopolFiles()
