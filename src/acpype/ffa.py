"""GROMACS force field addon (``.ffa``) writer for ACPYPE.

Writes a ``<basename>.ffa`` directory that extends a base force field with
the parameters of a single ligand residue, for use with
``gmx pdb2gmx -ff <base> -ffa <basename>`` (GROMACS 2027 and later, see
the ``extendable-force-field`` merge requests).

The addon carries:

``compatibility.dat``
    Base force fields the addon may be used with (``*`` = any).
``residuetypes.dat``
    Maps the ligand residue name to a residue type of its own, so pdb2gmx
    makes it a molecule on its own.
``atomtypes.atp``
    The new atom types with their masses.
``forcefield-addon.itp``
    The full Lennard-Jones definitions plus ``#include "ffbonded.itp"``;
    pdb2gmx includes this file into the generated topology directly after
    the base force field, so it is the only path through which the addon
    bonded type tables reach grompp.  It starts with a ``[ addon_priority ]``
    declaration so that the addon definitions win conflicts with the base
    force field; grompp propagates it into the included files.
``ffbonded.itp``
    Bonded *type* definitions (``[ bondtypes ]``, ``[ angletypes ]``,
    ``[ dihedraltypes ]``) merged into pdb2gmx's type tables. Proper
    torsions are written as multi-term dihedral type 9 blocks, one
    consecutive block per atom-type quartet, matching the AMBER Fourier
    terms of that torsion.
``<residue>.rtp``
    The residue topology: atom names, types and charges plus the bonded
    interaction connectivity (atom names only; parameters are resolved
    from the type tables above).
"""

from __future__ import annotations

import math
import os

from acpype.errors import UnsupportedTopologyError
from acpype.params import cal, radPi


def _type_name(atom) -> str:
    return atom.atomType.atomTypeName


def _sigma_epsilon(a_type):
    """Sigma (nm) and epsilon (kJ/mol) from AMBER A/B coefficients."""
    a_coef = a_type.ACOEF
    b_coef = a_type.BCOEF
    if b_coef == 0.0:
        return 0.0, 0.0
    sigma = 0.1 * math.pow(a_coef / b_coef, 1.0 / 6)
    epsilon = 0.25 * b_coef * b_coef / a_coef * cal
    return sigma, epsilon


def _canonical(types: tuple[str, ...]):
    """Return the lexicographically smaller of a type tuple and its reverse."""
    return min(types, types[::-1])


def write_ffa(mol):
    """Write the ``<mol.baseName>.ffa`` force field addon for one ligand."""
    if mol.cmaps:
        raise UnsupportedTopologyError(
            "the force field addon writer does not support CMAP terms; "
            "use the regular GROMACS topology output for this molecule"
        )

    solvent_names = (mol.freeIonResidues() if mol.amb2gmx else set()) | {"WAT"}
    solute_atoms = [a for a in mol.atomsGromacs if mol.residueLabel[a.resid] not in solvent_names]
    solute_residues = {mol.residueLabel[a.resid] for a in solute_atoms}
    if len(solute_residues) != 1:
        raise UnsupportedTopologyError(
            f"the force field addon writer supports exactly one solute residue, found {sorted(solute_residues)}"
        )
    residue = solute_residues.pop()

    if mol.topo14Data.hasNondefault14():
        mol.printWarn(
            "Non-default 1-4 scale factors detected. The force field addon applies "
            "the uniform scaling of the base force field; per-torsion factors are "
            "not carried over.\n"
        )

    addon_dir = os.path.join(os.path.abspath("."), mol.baseName + ".ffa")
    os.makedirs(addon_dir, exist_ok=True)
    stamp = f"; Written by ACPYPE {mol.baseName} for gmx pdb2gmx -ffa\n"

    def write(name: str, text):
        path = os.path.join(addon_dir, name)
        with open(path, "w") as fhd:
            fhd.writelines(text)
        mol.printMess(f"  wrote {os.path.relpath(path)}\n")

    # compatibility.dat: allow any base force field; tighten by hand if needed.
    write("compatibility.dat", "*\n")

    # residuetypes.dat: the residue becomes a molecule of its own.
    write("residuetypes.dat", f"{residue} {residue}\n")

    # atomtypes.atp: new atom types with masses.
    atp = [f"{t.atomTypeName:<8s} {t.mass:9.5f}\n" for t in mol.atomTypesGromacs]
    write("atomtypes.atp", atp)

    # forcefield-addon.itp: full atom type definitions, included into the
    # generated topology right after the base force field.
    # [ addon_priority ] makes the addon definitions win type conflicts with the
    # base force field, and grompp propagates it into the included files below.
    itp = [
        stamp,
        "\n[ addon_priority ]\n10\n",
        "\n[ atomtypes ]\n",
        ";name   bond_type     mass     charge   ptype   sigma         epsilon\n",
    ]
    for a_type in mol.atomTypesGromacs:
        sigma, epsilon = _sigma_epsilon(a_type)
        if a_type.BCOEF == 0.0:
            mol.printWarn(
                f"Atom type {a_type.atomTypeName} has zero van der Waals "
                "epsilon; it will be non-interacting in the addon.\n"
            )
        itp.append(
            f" {a_type.atomTypeName:<8s} {a_type.atomTypeName:<11s} {a_type.mass:9.5f}  0.00000   A"
            f"   {sigma:13.5e} {epsilon:13.5e}\n"
        )
    # The bonded type tables must reach the topology through this file as
    # well: pdb2gmx includes only forcefield-addon.itp into the generated
    # top, while its own read of the addon's ffbonded.itp only fills the
    # internal lookup tables.
    itp.append('\n#include "ffbonded.itp"\n')
    write("forcefield-addon.itp", itp)

    # ffbonded.itp: bonded types used by the residue.
    bonded = [stamp]

    bond_types = {}
    for bond in mol.bonds:
        types = tuple(_type_name(a) for a in bond.atoms)
        key = _canonical(types)
        params = (bond.rEq * 0.1, bond.kBond * 200 * cal)
        params = (round(params[0], 8), round(params[1], 4))
        seen = bond_types.setdefault(key, [])
        if params not in seen:
            if seen:
                mol.printWarn(
                    f"Bond types {key[0]}-{key[1]} have more than one parameter set; pdb2gmx will use the first one.\n"
                )
            seen.append(params)
    bonded.append("\n[ bondtypes ]\n; i    j   funct   length    force.c.\n")
    for types, param_list in sorted(bond_types.items()):
        orders = [types, types[::-1]] if types[::-1] != types else [types]
        for order in orders:
            for r_eq, k_bond in param_list:
                bonded.append(f"{order[0]:<6s} {order[1]:<6s} 1   {r_eq:.4e}   {k_bond:.4e}\n")

    angle_types = {}
    for angle in mol.angles:
        types = tuple(_type_name(a) for a in angle.atoms)
        key = _canonical(types)
        params = (round(angle.thetaEq * radPi, 4), round(2 * cal * angle.kTheta, 4))
        seen = angle_types.setdefault(key, [])
        if params not in seen:
            if seen:
                mol.printWarn(
                    f"Angle types {key[0]}-{key[1]}-{key[2]} have more than one parameter set; "
                    "pdb2gmx will use the first one.\n"
                )
            seen.append(params)
    bonded.append("\n[ angletypes ]\n; i    j    k  funct   angle     force.c.\n")
    for types, param_list in sorted(angle_types.items()):
        orders = [types, types[::-1]] if types[::-1] != types else [types]
        for order in orders:
            for theta, k_theta in param_list:
                bonded.append(f"{order[0]:<6s} {order[1]:<6s} {order[2]:<6s} 1   {theta:.3f}   {k_theta:.4e}\n")

    # Proper torsions: one multi-term (type 9) block per atom-type quartet,
    # with the Fourier terms of the AMBER torsion on consecutive lines.
    dihedral_blocks = {}
    for dih in mol.properDihedralsGmx45:
        types = tuple(_type_name(a) for a in dih[0])
        key = _canonical(types)
        term = (round(dih[1], 4), round(dih[2] * cal, 6), int(dih[3]))
        seen = dihedral_blocks.setdefault(key, [])
        if term not in seen:
            seen.append(term)
    bonded.append("\n[ dihedraltypes ]\n;i    j    k    l   funct   phase      kd        pn\n")
    for types, terms in sorted(dihedral_blocks.items()):
        for phase, kd, pn in terms:
            bonded.append(
                f"{types[0]:<6s} {types[1]:<6s} {types[2]:<6s} {types[3]:<6s} 9   {phase:8.2f}   {kd:9.5f}   {pn:2d}\n"
            )

    # Impropers (periodic improper function 4) and the AMBER alpha/gamma
    # special propers, which the regular output also writes as function 4.
    improper_types = {}
    for dih in mol.improperDihedrals:
        types = tuple(_type_name(a) for a in dih.atoms)
        key = _canonical(types)
        term = (round(dih.phase * radPi, 4), round(dih.kPhi * cal, 6), int(dih.period))
        seen = improper_types.setdefault(key, [])
        if term not in seen:
            if seen:
                mol.printWarn(
                    f"Improper types {key} have more than one term; GROMACS takes "
                    "the first one. Merge the terms by hand if that is not intended.\n"
                )
            seen.append(term)
    for dih in mol.properDihedralsAlphaGamma:
        types = tuple(_type_name(a) for a in dih[0])
        key = _canonical(types)
        term = (round(dih[1], 4), round(dih[2] * cal, 6), int(dih[3]))
        seen = improper_types.setdefault(key, [])
        if term not in seen:
            seen.append(term)
    for types, terms in sorted(improper_types.items()):
        for phase, kd, pn in terms:
            bonded.append(
                f"{types[0]:<6s} {types[1]:<6s} {types[2]:<6s} {types[3]:<6s} 4   {phase:8.2f}   {kd:9.5f}   {pn:2d}\n"
            )
    write("ffbonded.itp", bonded)

    # <residue>.rtp: residue topology referencing the types above by name.
    rtp = [
        stamp,
        "\n[ bondedtypes ]\n",
        "; bonds  angles  dihedrals  impropers  hydrogens  nrexcl\n",
        "     1       1          9          4          1       3\n",
        f"\n[ {residue} ]\n",
        " [ atoms ]\n",
        ";name   atomtype   charge    cgenr\n",
    ]
    for atom in solute_atoms:
        rtp.append(f" {atom.atomName:<6s} {atom.atomType.atomTypeName:<6s} {atom.charge:10.6f}   {atom.cgnr}\n")

    def emit_section(header: str, items):
        if not items:
            return
        rtp.append(f" [ {header} ]\n")
        for combo in items:
            rtp.append("  " + "  ".join(f"{name:<6s}" for name in combo) + "\n")

    solute_ids = {a.id for a in solute_atoms}
    bonds = sorted({tuple(_atom_name(a) for a in b.atoms) for b in mol.bonds if {a.id for a in b.atoms} <= solute_ids})
    angles = sorted(
        {tuple(_atom_name(a) for a in ang.atoms) for ang in mol.angles if {a.id for a in ang.atoms} <= solute_ids}
    )
    emit_section("bonds", bonds)
    emit_section("angles", angles)

    seen_quartets = set()
    dihedral_combos = []
    for dih in mol.properDihedralsGmx45:
        combo = tuple(_atom_name(a) for a in dih[0])
        if combo in seen_quartets:
            continue
        seen_quartets.add(combo)
        dihedral_combos.append(combo)
    emit_section("dihedrals", dihedral_combos)

    seen_quartets = set()
    improper_combos = []
    for dih_atoms in [dih.atoms for dih in mol.improperDihedrals] + [dih[0] for dih in mol.properDihedralsAlphaGamma]:
        combo = tuple(_atom_name(a) for a in dih_atoms)
        if combo in seen_quartets:
            continue
        seen_quartets.add(combo)
        improper_combos.append(combo)
    emit_section("impropers", improper_combos)

    write(f"{residue}.rtp", rtp)

    write(
        "README.txt",
        [
            f"Force field addon {mol.baseName!r}: residue {residue} with GAFF/AMBER parameters,\n"
            "generated by ACPYPE for gmx pdb2gmx's -ffa option (GROMACS 2027+).\n"
            "\n"
            "Usage (run next to this directory, or give its path after -ffa):\n"
            "\n"
            f"  gmx pdb2gmx -ff <base-force-field> -ffa {mol.baseName} -f <ligand>.pdb -o out.gro -p topol.top\n"
            "\n"
            "The ligand residue in the input structure must be named "
            f"{residue!r}.\n",
        ],
    )

    mol.printMess(f"Force field addon written to {addon_dir}\n")


def _atom_name(atom) -> str:
    return atom.atomName
