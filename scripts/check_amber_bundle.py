#!/usr/bin/env python3
"""Verify every executable in a vendored AmberTools bundle can actually be loaded.

``scripts/vendor_amber.py`` decides which shared libraries to bundle by walking the
link dependencies of the AmberTools executables. Two things sit outside that walk:

``charmmgen``
    Modern AmberTools dropped it, so ACPYPE keeps a binary from an old release in
    ``charmmgen_{linux,macos}.tgz`` purely for legacy compatibility. It is untarred
    *after* vendoring and is not part of the closure, so nothing guarantees that the
    libraries it needs are present.

libraries left to the host
    On Linux the bundle deliberately relies on the system for libstdc++, BLAS and
    friends (see ``SYSTEM_LIBS_LINUX``), so a bundle can look complete on the build
    machine and still fail elsewhere.

Both failure modes show up the same way: the dynamic loader refuses to start the
program. This script runs each executable and fails if the loader complains, which is
exactly the breakage that made AmberTools unusable on Apple Silicon.
"""

import argparse
import subprocess
import sys
from pathlib import Path

#: Substrings that mean the dynamic loader, not the program, rejected the binary.
LOADER_ERRORS = (
    "Library not loaded",
    "image not found",
    "error while loading shared libraries",
    "cannot open shared object file",
    "symbol not found",
    "Bad CPU type",  # an x86_64 binary on Apple Silicon with no Rosetta installed
)


def is_executable(path: Path) -> bool:
    """Report whether a path is a regular file with the executable bit set.

    Args:
        path: Candidate file.

    Returns:
        ``True`` for a non-symlink executable file.
    """
    return path.is_file() and not path.is_symlink() and path.stat().st_mode & 0o111 != 0


def loader_failure(path: Path, amberhome: Path) -> str | None:
    """Run one executable and report a dynamic-loader failure, if any.

    The programs are run with no arguments; most print usage and exit non-zero, which
    is fine. Only loader diagnostics count as failures.

    Args:
        path: Executable to run.
        amberhome: Value to expose as ``AMBERHOME``.

    Returns:
        The offending output, or ``None`` when the binary loaded.
    """
    try:
        result = subprocess.run(
            [str(path)],
            capture_output=True,
            timeout=60,
            env={"AMBERHOME": str(amberhome), "PATH": f"{amberhome}/bin"},
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return f"could not execute: {exc}"
    output = (result.stdout + result.stderr).decode("utf-8", errors="replace")
    for marker in LOADER_ERRORS:
        if marker in output:
            return next((ln for ln in output.splitlines() if marker in ln), marker).strip()
    return None


def main(argv: list[str] | None = None) -> int:
    """Check every executable in a vendored bundle.

    Args:
        argv: Command line arguments; defaults to ``sys.argv[1:]``.

    Returns:
        Process exit status: non-zero if any executable failed to load.
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("root", type=Path, help="bundle to check, e.g. src/acpype/amber_linux")
    args = parser.parse_args(argv)

    if not (args.root / "bin").is_dir():
        parser.error(f"{args.root} does not look like a vendored AmberTools bundle")

    root = args.root.resolve()
    failures: list[tuple[Path, str]] = []
    checked = 0
    for path in sorted((root / "bin").rglob("*")):
        if not is_executable(path):
            continue
        checked += 1
        problem = loader_failure(path, root)
        if problem:
            failures.append((path.relative_to(root), problem))

    for rel, problem in failures:
        print(f"  FAIL {rel}: {problem}", file=sys.stderr)
    print(f"checked {checked} executables, {len(failures)} failed to load")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
