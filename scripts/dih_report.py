#!/usr/bin/env python3
"""List every proper dihedral where an ACPYPE topology disagrees with its GROMACS reference.

Companion to ``check_acpype.py``, which reports *that* the bonded energies of a
tripeptide disagree; this says *which* dihedrals, on which atoms, with which parameters
on each side. It reads the ``.tpr`` pairs the harness leaves behind, so the parameters
compared are the ones GROMACS actually resolved after every wildcard and include, not
anything re-parsed from topology text.

For each peptide it prints the proper-dihedral energy on both sides and the number of
offending quartets, then the offending quartets grouped by atom-type pattern across all
peptides, then the per-peptide detail with atom names and types. Where the reference
typed an atom differently from ACPYPE, both patterns are shown: that difference is
usually the whole explanation, as with tyrosine's CZ.

Usage:
    uv run python scripts/dih_report.py WORKDIR [--codes AAA,GGG] [--out report.txt]

``WORKDIR`` is the directory ``check_acpype.py --dir`` wrote: ``aXXX.tpr``/``aXXX.edr``
and ``aXXXp.top`` for the GROMACS reference, ``agXXX.tpr``/``agXXX.edr`` plus
``agXXX.acpype/`` for ACPYPE. Needs ``gmx`` on ``PATH`` or in ``$GMX``.
"""

import argparse
import collections
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

AA = {
    "A": "ALA",
    "C": "CYS",
    "D": "ASP",
    "E": "GLU",
    "F": "PHE",
    "G": "GLY",
    "H": "HIE",
    "J": "HIP",
    "O": "HID",
    "I": "ILE",
    "K": "LYS",
    "L": "LEU",
    "M": "MET",
    "N": "ASN",
    "P": "PRO",
    "Q": "GLN",
    "R": "ARG",
    "S": "SER",
    "T": "THR",
    "V": "VAL",
    "W": "TRP",
    "Y": "TYR",
}

# (phase in degrees, force constant in kJ/mol, multiplicity)
Term = tuple[int, float, int]
Quartet = tuple[int, int, int, int]
AtomInfo = tuple[str, str, str]  # atom type, residue name, atom name
UNKNOWN: AtomInfo = ("?", "?", "?")


def gmx() -> str:
    """Return the GROMACS executable to use, or exit with a clear message."""
    exe = os.environ.get("GMX") or shutil.which("gmx") or shutil.which("gmx_mpi")
    if not exe:
        sys.exit("ERROR: no 'gmx' found. Install GROMACS or point $GMX at it.")
    return exe


def atoms_from_top(path: Path) -> dict[int, AtomInfo]:
    """Read ``1-based id -> (type, residue, atom name)`` from a topology's ``[ atoms ]``.

    Args:
        path: a ``.top`` or ``.itp`` holding one moleculetype.

    Returns:
        The atoms keyed by their index within the moleculetype; empty if the file is
        missing.
    """
    atoms: dict[int, AtomInfo] = {}
    if not path.is_file():
        return atoms
    section = None
    for line in path.read_text().splitlines():
        header = re.match(r"^\s*\[\s*(.+?)\s*\]", line)
        if header:
            section = header.group(1)
            continue
        if section != "atoms" or not line.strip() or line.lstrip()[0] in ";#":
            continue
        tokens = line.split()
        try:
            atoms[int(tokens[0])] = (tokens[1], tokens[3], tokens[4])
        except (ValueError, IndexError):
            continue
    return atoms


def pdihs_from_tpr(path: Path) -> dict[Quartet, list[Term]]:
    """Read the resolved proper dihedrals of a run input file.

    Args:
        path: the ``.tpr``.

    Returns:
        For each atom quartet, 1-based and in a canonical direction, the sorted list of
        Fourier terms GROMACS assigned to it.
    """
    dump = subprocess.run([gmx(), "dump", "-s", str(path)], capture_output=True, text=True, check=False).stdout
    types: dict[int, Term] = {}
    grouped: dict[Quartet, list[Term]] = collections.defaultdict(list)
    for line in dump.splitlines():
        functype = re.search(
            r"functype\[(\d+)\]=PDIHS, phiA=\s*([-\d.e+]+), cpA=\s*([-\d.e+]+).*mult=(\d+)",
            line,
        )
        if functype:
            index, phase, k, mult = functype.groups()
            types[int(index)] = (round(float(phase)), round(float(k), 3), int(mult))
            continue
        instance = re.match(r"\s*\d+ type=(\d+) \(PDIHS\)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)", line)
        if instance:
            a, b, c, d = (int(instance.group(i)) + 1 for i in range(2, 6))
            quartet: Quartet = (a, b, c, d) if a <= d else (d, c, b, a)
            grouped[quartet].append(types[int(instance.group(1))])
    return {quartet: sorted(terms) for quartet, terms in grouped.items()}


def proper_dihedral_energy(edr: Path) -> float:
    """Return the ``Proper Dih.`` energy of an energy file, in kJ/mol."""
    result = subprocess.run(
        [gmx(), "energy", "-f", str(edr), "-o", str(edr.with_suffix(".dih.xvg"))],
        input="Proper-Dih.\n\n",
        capture_output=True,
        text=True,
        check=False,
    )
    for line in (result.stdout + result.stderr).splitlines():
        if line.startswith("Proper Dih."):
            return float(line.split()[-5])
    return float("nan")


def fmt(terms: list[Term]) -> str:
    """Render a Fourier series as ``k@phase°xn + ...`` or a dash when empty."""
    return " + ".join(f"{k:.3f}@{phase}°x{mult}" for phase, k, mult in terms) or "-"


def types_of(quartet: Quartet, atoms: dict[int, AtomInfo]) -> str:
    """Join the atom types of a quartet as ``A-B-C-D``."""
    return "-".join(atoms.get(i, UNKNOWN)[0] for i in quartet)


def type_pattern(quartet: Quartet, acp: dict[int, AtomInfo], ref: dict[int, AtomInfo]) -> str:
    """Name a quartet by ACPYPE's atom types, adding the reference's when they differ."""
    mine = types_of(quartet, acp)
    if not ref:
        return mine
    theirs = types_of(quartet, ref)
    return mine if theirs == mine else f"{mine} (ref {theirs})"


def compare(workdir: Path, codes: list[str] | None) -> str:
    """Build the report for every peptide found in ``workdir``.

    Args:
        workdir: the ``check_acpype.py`` working directory.
        codes: peptide codes to restrict to, e.g. ``["AAA", "YYY"]``; all when ``None``.

    Returns:
        The report text.
    """
    found = sorted(p.name[1:4] for p in workdir.glob("a???.tpr"))
    selected = [c for c in found if not codes or c in codes]
    if not selected:
        sys.exit(f"ERROR: no aXXX.tpr files in {workdir}; run check_acpype.py --dir {workdir} --keep first")

    rows = []
    by_pattern: collections.Counter[str] = collections.Counter()
    example: dict[str, tuple[list[Term], list[Term]]] = {}
    for code in selected:
        ref = pdihs_from_tpr(workdir / f"a{code}.tpr")
        acp = pdihs_from_tpr(workdir / f"ag{code}.tpr")
        names = atoms_from_top(workdir / f"ag{code}.acpype" / f"ag{code}_GMX.itp")
        # grompp -pp output of the reference. Its atom types can differ from ACPYPE's,
        # and when they do that is usually the whole story.
        ref_names = atoms_from_top(workdir / f"a{code}p.top")
        e_ref = proper_dihedral_energy(workdir / f"a{code}.edr")
        e_acp = proper_dihedral_energy(workdir / f"ag{code}.edr")
        diffs = [(q, ref.get(q, []), acp.get(q, [])) for q in sorted(set(ref) | set(acp)) if ref.get(q) != acp.get(q)]
        rows.append((code, e_ref, e_acp, diffs, names, ref_names))
        for quartet, a, b in diffs:
            pattern = type_pattern(quartet, names, ref_names)
            by_pattern[pattern] += 1
            example.setdefault(pattern, (a, b))

    out = []
    out.append("=" * 96)
    out.append(f"{'pep':4} {'res':4} {'Proper Dih. kJ/mol':>22}   {'gap':>7}  offending quartets")
    out.append(f"{'':4} {'':4} {'reference':>10} {'acpype':>10}")
    out.append("-" * 96)
    for code, e_ref, e_acp, diffs, _, _ in rows:
        gap = 100 * abs(e_acp - e_ref) / abs(e_ref) if e_ref else float("nan")
        out.append(f"{code:4} {AA.get(code[0], '?'):4} {e_ref:10.2f} {e_acp:10.2f}   {gap:6.1f}%  {len(diffs)}")

    out.append("")
    out.append("=" * 96)
    out.append("OFFENDING QUARTETS BY ATOM-TYPE PATTERN  (acpype types; reference's in brackets when different)")
    out.append("-" * 96)
    for pattern, count in by_pattern.most_common():
        a, b = example[pattern]
        out.append(f"{pattern:28} x{count:<3}  reference: {fmt(a)}")
        out.append(f"{'':28}       acpype   : {fmt(b)}")

    out.append("")
    out.append("=" * 96)
    out.append("PER-PEPTIDE DETAIL  (atom name + index (acpype type) ... : reference -> acpype)")
    out.append("-" * 96)
    for code, _, _, diffs, names, ref_names in rows:
        if not diffs:
            continue
        out.append(f"\n{code} ({AA.get(code[0], '?')}):")
        for quartet, a, b in diffs:
            label = "-".join(f"{names.get(i, UNKNOWN)[2]}{i}({names.get(i, UNKNOWN)[0]})" for i in quartet)
            if ref_names and types_of(quartet, ref_names) != types_of(quartet, names):
                label += f"   [reference types: {types_of(quartet, ref_names)}]"
            out.append(f"   {label}")
            out.append(f"        reference: {fmt(a)}")
            out.append(f"        acpype   : {fmt(b)}")
    return "\n".join(out) + "\n"


def main() -> None:
    """Parse the command line and write the report."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("workdir", type=Path, help="directory written by check_acpype.py --dir ... --keep")
    parser.add_argument("--codes", default=None, help="comma-separated peptide codes to restrict to, e.g. AAA,YYY")
    parser.add_argument("--out", type=Path, default=None, help="write the report here instead of stdout")
    args = parser.parse_args()

    codes = [c.strip().upper() for c in args.codes.split(",")] if args.codes else None
    report = compare(args.workdir, codes)
    if args.out:
        args.out.write_text(report)
        print(f"wrote {args.out}")
    else:
        sys.stdout.write(report)


if __name__ == "__main__":
    main()
