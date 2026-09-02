#!/usr/bin/env python3
"""Copy the AmberTools files ACPYPE ships out of a conda environment.

ACPYPE bundles a trimmed-down AmberTools so that ``antechamber`` and friends work
without a separate install. Historically the file list was hand-maintained in
``update_linux_bins.sh`` / ``update_macos_bins.sh``, which silently rotted: an
AmberTools upgrade renames shared libraries (``libnetcdf.19`` became ``libnetcdf.22``
in AmberTools 26) and the copy of a now-missing file was ignored, producing a
half-broken bundle.

This script instead derives the shared libraries from the executables themselves,
following ``NEEDED``/``LC_LOAD_DYLIB`` entries transitively and keeping any library
that the source environment actually provides. Only the executables and data
directories stay declared by hand, because those are a genuine product decision.
"""

import argparse
import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path

#: Executables ACPYPE invokes, relative to the environment root.
PROGRAMS = (
    "bin/am1bcc",
    "bin/antechamber",
    "bin/atomtype",
    "bin/bondtype",
    "bin/parmchk2",
    "bin/sqm",
    "bin/teLeap",
    "bin/tleap",
    "bin/wrapped_progs/am1bcc",
    "bin/wrapped_progs/antechamber",
    "bin/wrapped_progs/atomtype",
    "bin/wrapped_progs/bondtype",
    "bin/wrapped_progs/parmchk2",
)

# Note: charmmgen is deliberately absent. Modern AmberTools dropped it, so ACPYPE
# ships an old build from charmmgen_{linux,macos}.tgz for legacy compatibility; the
# update scripts untar it after this script runs. Because it is outside the closure
# below, scripts/check_amber_bundle.py exists to confirm it still loads.

#: Data directories and stray files copied verbatim.
DATA = (
    "dat/antechamber",
    "dat/chamber",
    "dat/leap",
    "LICENSE",
    "GNU_LGPL_v3",
    "amber.sh",
)

#: Excluded from the copied data directories.
EXCLUDE = shutil.ignore_patterns("__pycache__", "*.pyc", "*.pyo", "*.log", "pixmaps")

#: Libraries left to the host on Linux.
#:
#: ELF binaries resolve their ``NEEDED`` entries through the system loader, so anything
#: a mainstream distribution ships can stay unbundled -- these are exactly the packages
#: the Dockerfile and CI already install. Bundling them instead would add well over
#: 100 MB (OpenBLAS alone is 40 MB, ICU 32 MB, libstdc++ 23 MB).
#:
#: macOS gets no such list: conda links its dylibs through ``@rpath``, so a dependency
#: that is not bundled simply fails to load.
#: Only sonames a distribution genuinely provides may appear here. conda-forge builds
#: some libraries with sonames no distro ships -- libxml2 is ``libxml2.so.16`` under
#: conda but ``libxml2.so.2`` on Debian -- and excluding those breaks teLeap at load
#: time. Bundling libxml2 in turn drags in ICU, which is where 35 MB of the Linux
#: bundle goes.
SYSTEM_LIBS_LINUX = (
    "libblas",
    "libc++",
    "libcurl",
    "libgcc_s",
    "libgfortran",
    "libgomp",
    "liblapack",
    "libopenblas",
    "libquadmath",
    "libstdc++",
)


def is_system_lib(name: str, is_macos: bool) -> bool:
    """Report whether a library should be left to the host rather than bundled.

    Args:
        name: Library file name, e.g. ``libstdc++.so.6.0.36``.
        is_macos: Whether the target platform is macOS.

    Returns:
        ``True`` when the host is expected to provide the library.
    """
    if is_macos:
        return False
    stem = name.split(".so")[0]
    return stem in SYSTEM_LIBS_LINUX


def run(cmd: list[str]) -> str:
    """Run a command and return stdout, tolerating non-UTF-8 bytes.

    Args:
        cmd: Command and arguments to execute.

    Returns:
        Standard output, decoded leniently. Empty when the command fails.
    """
    result = subprocess.run(cmd, capture_output=True, check=False)
    return result.stdout.decode("utf-8", errors="replace")


def needed_libraries(path: Path, is_macos: bool) -> list[str]:
    """List the shared-library names a binary links against.

    Args:
        path: Executable or shared library to inspect.
        is_macos: Whether to read Mach-O rather than ELF metadata.

    Returns:
        Bare library file names, without any directory part.
    """
    names: list[str] = []
    if is_macos:
        # Skip the first line, which is the file's own install name.
        for line in run(["otool", "-L", str(path)]).splitlines()[1:]:
            entry = line.strip().split(" (", 1)[0]
            if entry:
                names.append(os.path.basename(entry))
    else:
        for line in run(["objdump", "-p", str(path)]).splitlines():
            parts = line.split()
            if len(parts) == 2 and parts[0] == "NEEDED":
                names.append(parts[1])
    return names


def resolve_closure(source: Path, programs: list[Path], is_macos: bool) -> set[Path]:
    """Walk link dependencies transitively, keeping those the environment provides.

    A dependency is bundled only when it exists under ``<source>/lib``; anything else
    is assumed to come from the host system, which is how the previous hand-written
    lists behaved.

    Args:
        source: Root of the conda environment.
        programs: Executables to start the walk from.
        is_macos: Whether the environment holds Mach-O binaries.

    Returns:
        Absolute paths of the libraries to bundle, including symlinks encountered.
    """
    libdir = source / "lib"
    found: set[Path] = set()
    queue = list(programs)
    seen: set[Path] = set(programs)

    while queue:
        current = queue.pop()
        for name in needed_libraries(current, is_macos):
            candidate = libdir / name
            if not candidate.exists() or candidate in seen:
                continue
            if is_system_lib(name, is_macos):
                seen.add(candidate)
                continue
            seen.add(candidate)
            found.add(candidate)
            # Keep the symlink *and* its target, so the SONAME the loader asks for
            # and the real file both end up in the bundle.
            if candidate.is_symlink():
                target = candidate.resolve()
                if target.exists():
                    found.add(target)
                    if target not in seen:
                        seen.add(target)
                        queue.append(target)
            queue.append(candidate)
    return found


def copy(src: Path, dest: Path) -> None:
    """Copy a file, directory or symlink, creating parent directories as needed.

    Args:
        src: Path to copy from.
        dest: Path to copy to.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() or dest.is_symlink():
        if dest.is_dir() and not dest.is_symlink():
            shutil.rmtree(dest)
        else:
            dest.unlink()
    if src.is_symlink():
        os.symlink(os.readlink(src), dest)
    elif src.is_dir():
        shutil.copytree(src, dest, symlinks=True, ignore=EXCLUDE)
    else:
        shutil.copy2(src, dest)


def main(argv: list[str] | None = None) -> int:
    """Vendor AmberTools from a conda environment into the ACPYPE tree.

    Args:
        argv: Command line arguments; defaults to ``sys.argv[1:]``.

    Returns:
        Process exit status.
    """
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--source",
        type=Path,
        required=True,
        help="conda environment holding AmberTools",
    )
    parser.add_argument(
        "--dest", type=Path, required=True, help="target, e.g. acpype/amber_linux"
    )
    parser.add_argument(
        "--macos", action="store_true", default=platform.system() == "Darwin"
    )
    args = parser.parse_args(argv)

    if not (args.source / "bin").is_dir():
        parser.error(f"{args.source} does not look like a conda environment")

    programs = []
    missing = []
    for rel in PROGRAMS:
        path = args.source / rel
        (programs if path.exists() else missing).append(path if path.exists() else rel)
    if missing:
        parser.error("missing expected programs: " + ", ".join(map(str, missing)))

    libraries = resolve_closure(args.source, programs, args.macos)

    if args.dest.exists():
        shutil.rmtree(args.dest)

    for path in programs:
        copy(path, args.dest / path.relative_to(args.source))
    for path in sorted(libraries):
        copy(path, args.dest / "lib" / path.name)
    for rel in DATA:
        path = args.source / rel
        if path.exists():
            copy(path, args.dest / rel)
        else:
            print(f"  note: {rel} absent from {args.source}", file=sys.stderr)

    print(
        f"vendored {len(programs)} programs and {len(libraries)} libraries into {args.dest}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
