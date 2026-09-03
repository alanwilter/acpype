# Changelog

All notable changes to ACPYPE. Versions are dates, `YYYY.M.D`, stamped at release time
by `scripts/ver_today.sh`. This file starts at the first release after the long gap
that followed 2023.10.27.

## [2026.9.3] - 2026-09-03

Force field correctness. `-a amber` produced peptide backbone torsions that had been
wrong since before 2023.10.27, and `amb2gmx` silently dropped ff19SB's CMAP correction;
both are fixed and validated against GROMACS' own AMBER ports.

### Added

- `amb2gmx` converts ff19SB systems, CMAP backbone correction included. The `prmtop`
  CMAP grids become `[ cmaptypes ]` keyed the way GROMACS' own `amber19sb.ff` port keys
  them (`C-* N-ALA XC-ALA C-ALA N-*`, matched on atom type and residue name), plus a
  `[ cmap ]` section per moleculetype. Verified with GROMACS 2026.3: an ALA tripeptide
  reproduces `pdb2gmx -ff amber19sb` exactly (CMAP `-2.31952` kJ/mol on both sides).
- `amb2gmx` writes four-site waters (OPC, TIP4P-Ew) from the `prmtop`'s own geometry
  and charges, settles plus a `virtual_sites3`, instead of a three-site template that
  no longer matched the coordinates. A peptide in 789 OPC waters with ions reproduces
  `pdb2gmx -water opc`'s total potential to 0.001%.
- `-S/--split_molecules` (`amb2gmx`): one `[ moleculetype ]` per molecule AMBER
  identified in `ATOMS_PER_MOLECULE`, so a `combine { target ligand }` complex comes out
  as separate target and ligand blocks (#136). Charges are rounded per molecule, since
  AMBER's stored precision leaves an individual molecule off by ~1e-3. Bonded and
  Lennard-Jones energies are bit-identical to the merged topology; only the intended
  charge correction moves the electrostatics.
- `-F/--protein_ff ff14SB|ff99SB` selects the protein force field behind `-a amber` and
  `-a amber2`; ff14SB is the default.
- A clear, early error when antechamber assigns atom types the force field has no
  parameters for, such as `NO`/`DU` on a nitro group through `-a amber`. One line naming
  the atoms and suggesting `-a gaff2`, exit code 19, no traceback. tleap failures also
  get a one-line summary of the missing parameters before the raw log.
- `scripts/dih_report.py`: lists every proper dihedral on which an ACPYPE topology and
  its GROMACS reference disagree, with atom names, both sides' types and parameters.
- `scripts/check_acpype.py`, the 22-tripeptide validation harness behind `NOTE.txt`,
  runs again against a modern GROMACS and AmberTools, with `--protein_ff` and `--gmxff`.
- Characterisation tests that parse the emitted GROMACS topology, and tests that run
  `grompp -maxwarn 0` and compare energies with `pdb2gmx` when GROMACS is installed.

### Fixed

- `-a amber` produced peptide backbone torsions that disagreed with GROMACS' AMBER
  ports by 26-148%, on every one of the 22 standard tripeptides. Two causes, one
  origin: `-a amber` had moved to ff14SB, whose torsions are keyed on atom types
  (`CX`, `CO`, `2C`, `3C`, `C8`) antechamber never emits, so tleap fell back to zero
  wildcards; and parmchk2, given a parameter file with GAFF merged in, replaced AMBER
  terms it already had with GAFF analogues. ff14SB's residue types are now applied
  from AMBER's own templates after antechamber, and parmchk2 runs twice so GAFF fills
  only genuine gaps. Bonded energies now agree with `amber14sb.ff` to within 0.05% for
  21 residues and 1% for tyrosine, the exception on record since 2010. The fault
  predates 2023.10.27.
- `parseFrcmod` crashed on any frcmod with an empty section (`frcmod.ff99SB`) and swept
  `frcmod.ff19SB`'s trailing CMAP block into the non-bonded section as junk.
- `check_acpype.py` had been dead for years: hard-coded home-directory paths, an mdp
  still carrying its shell heredoc wrapper, `cutoff-scheme = group` (removed in GROMACS
  2020), an energy parser broken by the three-word `Per. Imp. Dih.`, and a failing
  residue that stranded the loop in its output directory.
- Sphinx could no longer import the package for its version after the `src/` move.
- The CLI help and README now state that `-q mopac` and `-q divcon` need an external
  AmberTools; the bundled one, like conda-forge's, only carries `sqm` (#139).

### Changed

- `-a amber` means ff14SB plus GAFF: AMBER types where antechamber assigns them,
  ff14SB's own residue types applied afterwards, GAFF only where AMBER has no
  parameter. `-F ff99SB` restores the original 2010 parameter set.
- The GROMACS 4 flavour (`-z`) refuses a topology with CMAP terms rather than dropping
  them; the CNS and CHARMM writers warn that they drop CMAP.
- The release workflow publishes from the same run that creates the tag and Release,
  reads the tag from `pyproject.toml` instead of the clock, and attaches PEP 740
  attestations that verify.

## [2026.9.2] - 2026-09-02

The first release in almost three years. The tooling was rebuilt end to end; the
scientific core is the same code with the fixes listed below.

### Added

- `-c abcg2`, the ABCG2 charge method new in AmberTools 26 (#129, PR #149).
- `-r/--predindex`, antechamber's atom/bond type prediction index, for molecules whose
  default perception fails (#26, PR #148).
- Python 3.12, 3.13 and 3.14 (#146). CI tests all three.
- AmberTools 26 vendored, replacing AmberTools 22, on both platforms. `charmmgen`,
  dropped from modern AmberTools, is rebuilt from the last release that shipped it, as
  a universal binary on macOS and natively on Linux.
- Wheels per platform (macOS arm64, Linux x86_64) plus a slim source distribution, to
  stay under PyPI's 100 MB per-file limit; the bundled AmberTools kept 24 shared-library
  symlinks the old build dropped.
- Publishing to PyPI through trusted publishing (OIDC) and to Docker Hub
  (`acpype/acpype:<version>` and `:latest`) from CI on release.
- `ty` type checking, expanded `ruff` rules, `uv audit`, a `justfile`, coverage
  reported through OIDC with a 90% floor, Dependabot.
- README section on what one `acpype` command replaces (#144).

### Changed

- **Breaking:** the build moved from poetry to uv (`uv_build` backend, `uv.lock`);
  Python 3.12 or newer is required.
- **Breaking:** the package lives under `src/acpype/`. Imports are unchanged; source
  checkouts that added the repository root to `sys.path` are not.
- **Breaking:** the command line is rebuilt on typer with rich help and error output.
  Every flag and short form is preserved, option values are now validated, and usage
  errors exit with code 2.
- The default branch is `main`.
- `mypy` and CodeQL are gone; `ty` and the audit step replace them.

### Fixed

- AmberTools would not start on Apple Silicon: duplicated `LC_RPATH` entries in the
  vendored binaries made dyld refuse them. The vendoring script now de-duplicates and
  re-signs.
- The GROMACS `[ atomtypes ]` block reported a mass of zero for every type (#132;
  fix by @davide-grheco).
- `parmMerge` cached merged parameter files by name only, so an AmberTools upgrade
  silently reused merges built from the previous release's parameters. The cache is
  now keyed on the inputs' content.
- `scripts/ver_today.sh` had stopped stamping the version at all (Homebrew dropped
  `pcregrep`, and the glob missed `src/`), and never failed; it now uses `grep`, names
  its files, refreshes `uv.lock`, and exits non-zero on trouble.
- Paths with spaces in the AmberTools location broke the binaries' wrappers.
- Tests ran in the shared `tests/` directory and collided; each now runs in its own
  temporary directory, and coverage no longer varies with a warm `/tmp` cache.
- CLI error assertions no longer depend on whether the terminal renders colour.

### Closed

- #24 and #128, reviewed and closed as resolved by the rebuilt toolchain.

## [2023.10.27]

Last release before this changelog.
