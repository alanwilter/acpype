"""Tests for converting CMAP terms and four-site water from AMBER to GROMACS.

ff19SB keeps its backbone correction in per-residue CMAP grids and is meant to be run
with OPC, a four-site water. Neither existed in ACPYPE's writer, which silently dropped
the CMAP and wrote a three-site water against four-site coordinates. The fixtures are
built with the bundled tleap so every test runs against real AmberTools output.

With GROMACS on ``PATH`` the energies are also compared against ``pdb2gmx`` on the very
same coordinates; those tests skip otherwise.
"""

import re
import shutil
import subprocess
from pathlib import Path

import pytest

from acpype.errors import UnsupportedTopologyError
from acpype.topol import MolTopol

KCAL = 4.184
# Empty rather than None so the subprocess argument lists stay lists of str.
GMX = shutil.which("gmx") or ""

SPE_MDP = """integrator = md
nsteps = 0
cutoff-scheme = Verlet
nstlist = 1
verlet-buffer-tolerance = -1
rlist = 1.0
rcoulomb = 1.0
rvdw = 1.0
coulombtype = {coulombtype}
pbc = xyz
nstcalcenergy = 1
nstenergy = 1
continuation = yes
"""


def tleap(script: str, prefix: str) -> tuple[str, str]:
    """Run the bundled tleap on ``script`` and return the prmtop and inpcrd it wrote."""
    exe = shutil.which("tleap")
    if exe is None:
        pytest.skip("no tleap available to build the fixture")
    Path(f"{prefix}.leap.in").write_text(script)
    subprocess.run([exe, "-f", f"{prefix}.leap.in"], capture_output=True, text=True, check=False)
    assert Path(f"{prefix}.prmtop").is_file(), f"tleap did not write {prefix}.prmtop"
    return f"{prefix}.prmtop", f"{prefix}.inpcrd"


def gmx_energies(tpr_prefix: str, coords: str, top: str, coulombtype: str) -> dict[str, float]:
    """Single-point energies of ``top`` on ``coords`` with GROMACS, by term name."""
    Path("spe.mdp").write_text(SPE_MDP.format(coulombtype=coulombtype))
    grompp = subprocess.run(
        [GMX, "grompp", "-c", coords, "-p", top, "-f", "spe.mdp", "-o", f"{tpr_prefix}.tpr", "-maxwarn", "0"],
        capture_output=True,
        text=True,
        check=False,
    )
    assert grompp.returncode == 0, grompp.stderr[-2500:]
    subprocess.run(
        [GMX, "mdrun", "-deffnm", tpr_prefix, "-s", f"{tpr_prefix}.tpr", "-nt", "2"], capture_output=True, check=False
    )
    result = subprocess.run(
        [GMX, "energy", "-f", f"{tpr_prefix}.edr", "-o", f"{tpr_prefix}.xvg"],
        input="Bond\nProper-Dih.\nCMAP-Dih.\nLJ-(SR)\nCoulomb-(SR)\nPotential\n\n",
        capture_output=True,
        text=True,
        check=False,
    )
    energies = {}
    for line in (result.stdout + result.stderr).splitlines():
        if line.endswith("(kJ/mol)"):
            tokens = line.split()
            energies[" ".join(tokens[:-5])] = float(tokens[-5])
    return energies


@pytest.fixture
def ff19sb_tripeptide(janitor: list[str]) -> tuple[str, str]:
    """An ALA tripeptide in vacuum with one CMAP term on the middle residue."""
    script = "source leaprc.protein.ff19SB\nm = sequence { NALA ALA CALA }\nsaveamberparm m ala19.prmtop ala19.inpcrd\nsavepdb m ala19.pdb\nquit\n"
    top, crd = tleap(script, "ala19")
    assert "%FLAG CMAP_COUNT" in Path(top).read_text()
    return top, crd


@pytest.fixture
def ff19sb_in_opc(janitor: list[str]) -> tuple[str, str]:
    """The same tripeptide in a box of OPC water with a pair of ions."""
    script = (
        "source leaprc.protein.ff19SB\nsource leaprc.water.opc\nm = sequence { NALA ALA CALA }\n"
        "solvatebox m OPCBOX 8.0\naddions m Na+ 1\naddions m Cl- 1\n"
        "saveamberparm m sol19.prmtop sol19.inpcrd\nsavepdb m sol19.pdb\nquit\n"
    )
    return tleap(script, "sol19")


def test_cmap_grids_are_read(ff19sb_tripeptide: tuple[str, str], janitor: list[str]) -> None:
    """One 24x24 grid in kcal/mol and one term on the middle residue's C-N-CA-C-N."""
    top, crd = ff19sb_tripeptide
    molecule = MolTopol(acFileXyz=crd, acFileTop=top, amb2gmx=True, verbose=False)
    janitor.append(molecule.absHomeDir)

    assert molecule.cmaps == [((11, 13, 15, 21, 23), 1)]
    side, grid = molecule.cmapGrids[0]
    assert side == 24
    assert len(grid) == 576
    assert grid[0] == pytest.approx(-0.4049)


def test_gromacs_topology_carries_the_cmap(ff19sb_tripeptide: tuple[str, str], janitor: list[str]) -> None:
    """The topology gets a residue-qualified [ cmaptypes ] entry and a [ cmap ] term, values in kJ/mol."""
    top, crd = ff19sb_tripeptide
    molecule = MolTopol(acFileXyz=crd, acFileTop=top, amb2gmx=True, verbose=False)
    janitor.append(molecule.absHomeDir)
    molecule.writeGromacsTopolFiles()
    text = (Path(molecule.absHomeDir) / "ala19_GMX.top").read_text()

    assert "\n[ cmaptypes ]\n" in text
    header = re.search(r"^C-\* N-ALA XC-ALA C-ALA N-\* 1 24 24\\\n([-\d.]+) ", text, re.M)
    assert header, "cmaptypes key must be the one GROMACS' amber19sb port uses"
    assert float(header.group(1)) == pytest.approx(-0.4049 * KCAL, abs=1e-6)
    assert re.search(r"\[ cmap \]\n;[^\n]*\n\s+11\s+13\s+15\s+21\s+23\s+1\n", text)
    assert text.index("[ cmaptypes ]") < text.index("[ moleculetype ]"), "cmaptypes are file-level parameters"


def test_gromacs4_flavour_refuses_cmap(ff19sb_tripeptide: tuple[str, str], janitor: list[str]) -> None:
    """GROMACS 4 output cannot carry CMAP, so asking for it is refused rather than degraded."""
    top, crd = ff19sb_tripeptide
    molecule = MolTopol(acFileXyz=crd, acFileTop=top, amb2gmx=True, gmx4=True, verbose=False)
    janitor.append(molecule.absHomeDir)

    with pytest.raises(UnsupportedTopologyError, match="GROMACS 4"):
        molecule.writeGromacsTopolFiles()


def test_four_site_water_is_built_from_the_prmtop(ff19sb_in_opc: tuple[str, str], janitor: list[str]) -> None:
    """OPC comes out as a four-atom moleculetype with settles and a virtual site matching GROMACS' opc.itp."""
    top, crd = ff19sb_in_opc
    molecule = MolTopol(acFileXyz=crd, acFileTop=top, amb2gmx=True, verbose=False)
    janitor.append(molecule.absHomeDir)
    molecule.writeGromacsTopolFiles()
    text = (Path(molecule.absHomeDir) / "sol19_GMX.top").read_text()

    water = text.split("four-site water built from the prmtop")[1].split("[ moleculetype ]")[0]
    atoms = [
        line.split()
        for line in water.split("[ atoms ]")[1].split("#ifdef")[0].splitlines()
        if line.strip() and not line.startswith(";")
    ]
    assert [a[4] for a in atoms] == ["O", "H1", "H2", "EPW"]
    assert [float(a[6]) for a in atoms] == pytest.approx([0.0, 0.67914, 0.67914, -1.35828], abs=1e-5)
    settles = re.search(r"\[ settles \]\n;[^\n]*\n1\s+1\s+([\d.]+)\s+([\d.]+)", water)
    assert settles and float(settles.group(1)) == pytest.approx(0.0872433, abs=1e-6)
    assert float(settles.group(2)) == pytest.approx(0.1371205, abs=1e-6)
    vsite = re.search(r"\[ virtual_sites3 \]\n;[^\n]*\n4\s+1\s+2\s+3\s+1\s+([\d.]+)\s+([\d.]+)", water)
    assert vsite and float(vsite.group(1)) == pytest.approx(0.147722, abs=2e-6)
    assert re.search(r"\[ molecules \][\s\S]*WAT\s+\d+", text)


@pytest.mark.skipif(not GMX, reason="needs a GROMACS install")
def test_cmap_energy_matches_gromacs_amber19sb(ff19sb_tripeptide: tuple[str, str], janitor: list[str]) -> None:
    """On identical coordinates the CMAP energy equals GROMACS' own amber19sb port to 5 decimals."""
    top, crd = ff19sb_tripeptide
    molecule = MolTopol(acFileXyz=crd, acFileTop=top, amb2gmx=True, verbose=False)
    janitor.append(molecule.absHomeDir)
    molecule.writeGromacsTopolFiles()

    subprocess.run(
        [GMX, "pdb2gmx", "-ff", "amber19sb", "-water", "none", "-f", "ala19.pdb", "-o", "ref.pdb", "-p", "ref.top"],
        capture_output=True,
        check=False,
    )
    for src, dst in (("ala19.pdb", "a_box.pdb"), ("ref.pdb", "b_box.pdb")):
        subprocess.run(
            [GMX, "editconf", "-f", src, "-o", dst, "-bt", "cubic", "-d", "3.0"], capture_output=True, check=False
        )
    ours = gmx_energies("a", "a_box.pdb", str(Path(molecule.absHomeDir) / "ala19_GMX.top"), "Cut-off")
    theirs = gmx_energies("b", "b_box.pdb", "ref.top", "Cut-off")

    assert ours["CMAP Dih."] == pytest.approx(theirs["CMAP Dih."], abs=1e-5)
    assert ours["CMAP Dih."] != 0.0
    # pdb2gmx rewrites the coordinates it was given, so the two sides differ by PDB
    # rounding; on a small extended peptide that is a few 1e-4 kJ/mol.
    for term in ("Bond", "Proper Dih.", "Potential"):
        assert ours[term] == pytest.approx(theirs[term], abs=2e-3), term


@pytest.mark.skipif(not GMX, reason="needs a GROMACS install")
def test_solvated_ff19sb_matches_gromacs(ff19sb_in_opc: tuple[str, str], janitor: list[str]) -> None:
    """Peptide, OPC water and ions together reproduce pdb2gmx's energies on the same coordinates."""
    top, crd = ff19sb_in_opc
    molecule = MolTopol(acFileXyz=crd, acFileTop=top, amb2gmx=True, verbose=False)
    janitor.append(molecule.absHomeDir)
    molecule.writeGromacsTopolFiles()

    # tleap's PDB carries the box as a CRYST1 record, so both sides can use it directly.
    # pdb2gmx wants GROMACS names (SOL/OW/HW1/HW2/MW, NA, CL); acpype's ion templates
    # write NA+/CL- atom names, so its copy only uppercases the ions.
    ions = {"Na+": "NA", "Cl-": "CL", "K+": "K"}
    for_gmx, for_acpype = [], []
    for line in Path("sol19.pdb").read_text().splitlines(keepends=True):
        acp = line
        if line.startswith(("ATOM", "HETATM")):
            residue, name = line[17:20].strip(), line[12:16].strip()
            if residue == "WAT":
                name = {"O": "OW", "H1": "HW1", "H2": "HW2", "EPW": "MW"}.get(name, name)
                line = line[:12] + f"{name:^4}" + line[16:17] + "SOL" + line[20:]
            elif residue in ions:
                line = line[:12] + f"{ions[residue]:^4}" + line[16:17] + f"{ions[residue]:<3}" + line[20:]
                acp = acp[:12] + f"{name.upper():^4}" + acp[16:]
        for_gmx.append(line)
        for_acpype.append(acp)
    Path("sol19_gmx.pdb").write_text("".join(for_gmx))
    Path("sol19_acp.pdb").write_text("".join(for_acpype))
    pdb2gmx = subprocess.run(
        [GMX, "pdb2gmx", "-ff", "amber19sb", "-water", "opc", "-f", "sol19_gmx.pdb", "-o", "ref.pdb", "-p", "ref.top"],
        capture_output=True,
        text=True,
        check=False,
    )
    assert pdb2gmx.returncode == 0, pdb2gmx.stderr[-2000:]
    ours = gmx_energies("a", "sol19_acp.pdb", str(Path(molecule.absHomeDir) / "sol19_GMX.top"), "PME")
    theirs = gmx_energies("b", "ref.pdb", "ref.top", "PME")

    for term in ("Bond", "Proper Dih.", "CMAP Dih.", "LJ (SR)", "Coulomb (SR)", "Potential"):
        assert ours[term] == pytest.approx(theirs[term], rel=1e-4), term
