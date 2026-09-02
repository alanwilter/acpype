#!/usr/bin/env python3
"""Build the distributions ACPYPE publishes: two platform wheels and a slim sdist.

ACPYPE vendors AmberTools for both Linux and macOS. A single ``py3-none-any`` wheel
therefore carries both trees, which since AmberTools 26 comes to about 129 MB --
over PyPI's 100 MB per-file limit, and twice what any individual user needs.

This script builds the universal wheel once, then derives one wheel per platform by
dropping the other platform's tree and retagging. Each result is roughly half the
size and installs only where its binaries can actually run.

Platform tags are deliberately narrow:

``manylinux_2_35_x86_64``
    The vendored binaries themselves only need ``GLIBC_2.17``, but ACPYPE leaves
    ``libstdc++`` to the host and one of the bundled libraries needs
    ``GLIBCXX_3.4.30``, which first ships with GCC 12 (Ubuntu 22.04, Debian 12).
    Claiming ``manylinux_2_17`` would install on systems where teLeap cannot load.

``macosx_11_0_arm64``
    conda-forge builds AmberTools for ``osx-arm64`` with a macOS 11 floor. The
    bundle contains no x86_64 slice, so Intel Macs must not resolve this wheel.

The sdist has the vendored AmberTools stripped out. Left in it would be ~121 MB, over
the same PyPI limit; more importantly, without it ``pip install acpype`` on a platform
with no wheel silently resolves to the last release that had a universal wheel. A slim
sdist means those users get the current ACPYPE and supply their own AmberTools, which
is the conda route the README already documents.
"""

import argparse
import shutil
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

#: Platform tag -> the vendored tree that wheel keeps.
TARGETS = {
    "manylinux_2_35_x86_64": "amber_linux",
    "macosx_11_0_arm64": "amber_macos",
}

#: Every vendored tree, so the ones not kept can be dropped.
ALL_TREES = set(TARGETS.values())


def run(cmd: list[str]) -> None:
    """Run a command, raising if it fails.

    Args:
        cmd: Command and arguments to execute.

    Raises:
        SystemExit: If the command exits non-zero.
    """
    result = subprocess.run(cmd, check=False)
    if result.returncode != 0:
        raise SystemExit(f"command failed: {' '.join(cmd)}")


def build_universal(outdir: Path) -> Path:
    """Build the ordinary ``py3-none-any`` wheel.

    Args:
        outdir: Directory to place the wheel in.

    Returns:
        Path to the built wheel.

    Raises:
        SystemExit: If exactly one wheel is not produced.
    """
    run(["uv", "build", "--wheel", "--out-dir", str(outdir)])
    wheels = sorted(outdir.glob("*-py3-none-any.whl"))
    if len(wheels) != 1:
        raise SystemExit(f"expected one universal wheel in {outdir}, found {len(wheels)}")
    return wheels[0]


def specialise(wheel: Path, platform_tag: str, keep: str, outdir: Path) -> Path:
    """Derive a platform-specific wheel from the universal one.

    Args:
        wheel: The ``py3-none-any`` wheel to start from.
        platform_tag: Wheel platform tag to apply, e.g. ``macosx_11_0_arm64``.
        keep: Name of the vendored tree to retain, e.g. ``amber_macos``.
        outdir: Directory to write the finished wheel to.

    Returns:
        Path to the platform-specific wheel.

    Raises:
        SystemExit: If the repacked wheel cannot be located.
    """
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        run(["wheel", "unpack", str(wheel), "-d", str(work)])
        unpacked = next(work.iterdir())

        for tree in ALL_TREES - {keep}:
            victim = unpacked / "acpype" / tree
            if victim.is_dir():
                shutil.rmtree(victim)

        staging = work / "packed"
        staging.mkdir()
        # `wheel pack` regenerates RECORD, so removing files above stays consistent.
        run(["wheel", "pack", str(unpacked), "-d", str(staging)])
        packed = next(staging.glob("*.whl"))

        run(["wheel", "tags", "--platform-tag", platform_tag, "--remove", str(packed)])
        tagged = next(staging.glob("*.whl"))

        outdir.mkdir(parents=True, exist_ok=True)
        final = outdir / tagged.name
        if final.exists():
            final.unlink()
        shutil.move(str(tagged), final)
        return final


def build_slim_sdist(outdir: Path) -> Path:
    """Build an sdist with the vendored AmberTools removed.

    Args:
        outdir: Directory to write the finished sdist to.

    Returns:
        Path to the slim sdist.

    Raises:
        SystemExit: If exactly one sdist is not produced.
    """
    with tempfile.TemporaryDirectory() as tmp:
        staging = Path(tmp)
        run(["uv", "build", "--sdist", "--out-dir", str(staging)])
        sdists = sorted(staging.glob("*.tar.gz"))
        if len(sdists) != 1:
            raise SystemExit(f"expected one sdist in {staging}, found {len(sdists)}")
        fat = sdists[0]

        outdir.mkdir(parents=True, exist_ok=True)
        slim = outdir / fat.name
        if slim.exists():
            slim.unlink()
        with tarfile.open(fat) as src, tarfile.open(slim, "w:gz") as dst:
            for member in src:
                # Members look like `acpype-<version>/acpype/amber_linux/...`; drop the
                # tree's own directory entry as well as everything beneath it.
                parts = member.name.split("/")
                if len(parts) > 2 and parts[1] == "acpype" and parts[2] in ALL_TREES:
                    continue
                dst.addfile(member, src.extractfile(member) if member.isfile() else None)
        return slim


def main(argv: list[str] | None = None) -> int:
    """Build every distribution ACPYPE publishes.

    Args:
        argv: Command line arguments; defaults to ``sys.argv[1:]``.

    Returns:
        Process exit status.
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out-dir", type=Path, default=Path("dist"))
    parser.add_argument(
        "--keep-universal",
        action="store_true",
        help="also leave the combined py3-none-any wheel in place (not for publishing)",
    )
    parser.add_argument("--no-sdist", action="store_true", help="build only the wheels")
    args = parser.parse_args(argv)

    with tempfile.TemporaryDirectory() as staging:
        universal = build_universal(Path(staging))
        built = [specialise(universal, tag, keep, args.out_dir) for tag, keep in TARGETS.items()]
        if args.keep_universal:
            shutil.copy2(universal, args.out_dir / universal.name)

    if not args.no_sdist:
        built.append(build_slim_sdist(args.out_dir))

    limit = 100 * 1024 * 1024
    over = False
    for path in built:
        size = path.stat().st_size
        flag = "  OVER PyPI LIMIT" if size > limit else ""
        over = over or size > limit
        print(f"  {path.name}  {size / 1e6:.1f} MB{flag}")
    if over:
        print("at least one wheel exceeds PyPI's 100 MB per-file limit", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
