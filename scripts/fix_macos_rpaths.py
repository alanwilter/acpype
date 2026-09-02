#!/usr/bin/env python3
"""De-duplicate ``LC_RPATH`` load commands in the vendored macOS AmberTools binaries.

conda-forge's macOS builds of AmberTools -- and of the GCC runtime libraries they
link against, notably ``libgfortran`` and ``libquadmath`` -- ship Mach-O files that
declare the same rpath more than once, for example both ``@loader_path`` and
``@loader_path/``. Since macOS 13, dyld rejects such a binary outright with::

    Library not loaded: @rpath/libgfortran.5.dylib
    Reason: tried: '.../lib/libgfortran.5.dylib' (duplicate LC_RPATH '@loader_path')

which makes ``sqm``, ``teLeap`` and ``antechamber`` unusable on Apple Silicon even
though the libraries are present and correctly signed. The defect is independent of
the AmberTools version, so this fixup must be re-applied every time the binaries are
re-vendored by ``update_macos_bins.sh``.

Removing the redundant entries with ``install_name_tool`` invalidates the ad-hoc code
signature that arm64 requires, so each modified file is re-signed afterwards.
"""

import argparse
import subprocess
import sys
from pathlib import Path

RPATH_PREFIX = "path "
OFFSET_MARKER = " (offset "


def _run(cmd: list[str]) -> str:
    """Run a command and return its stdout, tolerating non-UTF-8 bytes.

    ``otool`` echoes strings lifted straight out of the binary, which are not always
    valid UTF-8, so decoding errors are replaced rather than raised.

    Args:
        cmd: Command and arguments to execute.

    Returns:
        Standard output, decoded leniently.
    """
    result = subprocess.run(cmd, capture_output=True, check=False)
    return result.stdout.decode("utf-8", errors="replace")


def is_macho(path: Path) -> bool:
    """Report whether ``path`` is a regular (non-symlink) Mach-O file.

    Args:
        path: Candidate file.

    Returns:
        ``True`` if the file is Mach-O and not a symlink.
    """
    if not path.is_file() or path.is_symlink():
        return False
    return _run(["file", "-b", str(path)]).startswith("Mach-O")


def read_rpaths(path: Path) -> list[str]:
    """Return the ``LC_RPATH`` values declared by a Mach-O file, in order.

    Args:
        path: Mach-O file to inspect.

    Returns:
        The declared rpaths, including any duplicates.
    """
    lines = _run(["otool", "-l", str(path)]).splitlines()
    rpaths: list[str] = []
    for index, line in enumerate(lines):
        if "LC_RPATH" not in line:
            continue
        # otool prints `cmd LC_RPATH`, `cmdsize N`, then `path <value> (offset N)`.
        for candidate in lines[index + 1 : index + 4]:
            stripped = candidate.strip()
            if stripped.startswith(RPATH_PREFIX):
                value = stripped[len(RPATH_PREFIX) :]
                marker = value.rfind(OFFSET_MARKER)
                rpaths.append(value[:marker] if marker != -1 else value)
                break
    return rpaths


def redundant_rpaths(rpaths: list[str]) -> list[str]:
    """Pick the rpath entries that repeat a directory already declared earlier.

    Two entries are considered equivalent when they differ only by trailing slashes,
    which is how dyld normalises them before reporting a duplicate.

    Args:
        rpaths: Declared rpaths, in declaration order.

    Returns:
        The entries to delete, preserving the first occurrence of each directory.
    """
    seen: set[str] = set()
    redundant: list[str] = []
    for rpath in rpaths:
        key = rpath.rstrip("/")
        if key in seen:
            redundant.append(rpath)
        else:
            seen.add(key)
    return redundant


def fix_file(path: Path, dry_run: bool = False) -> list[str]:
    """Delete duplicate rpaths from one Mach-O file and re-sign it.

    Args:
        path: Mach-O file to repair.
        dry_run: When ``True``, report what would change without modifying the file.

    Returns:
        The rpath values removed; empty when the file was already correct.
    """
    to_delete = redundant_rpaths(read_rpaths(path))
    if not to_delete or dry_run:
        return to_delete

    for rpath in to_delete:
        subprocess.run(
            ["install_name_tool", "-delete_rpath", rpath, str(path)],
            capture_output=True,
            text=True,
            check=False,
        )
    # install_name_tool invalidates the signature; arm64 requires at least an ad-hoc one.
    subprocess.run(
        ["codesign", "--force", "--sign", "-", str(path)],
        capture_output=True,
        check=False,
    )
    return to_delete


def main(argv: list[str] | None = None) -> int:
    """Repair every Mach-O file beneath the given directory.

    Args:
        argv: Command line arguments; defaults to ``sys.argv[1:]``.

    Returns:
        Process exit status.
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("root", type=Path, nargs="?", default=Path("src/acpype/amber_macos"))
    parser.add_argument(
        "-n",
        "--dry-run",
        action="store_true",
        help="report changes without applying them",
    )
    args = parser.parse_args(argv)

    if sys.platform != "darwin":
        parser.error("this fixup only applies to the macOS binaries and must run on macOS")
    if not args.root.is_dir():
        parser.error(f"no such directory: {args.root}")

    scanned = 0
    fixed = 0
    for path in sorted(args.root.rglob("*")):
        if not is_macho(path):
            continue
        scanned += 1
        removed = fix_file(path, dry_run=args.dry_run)
        if removed:
            fixed += 1
            verb = "would fix" if args.dry_run else "fixed"
            print(f"  {verb} {path.relative_to(args.root)}: removed {' '.join(removed)}")

    print(f"scanned {scanned} Mach-O files, de-duplicated rpaths in {fixed}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
