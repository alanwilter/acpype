#!/usr/bin/env bash

# Build the `charmmgen` binaries ACPYPE ships in charmmgen_{macos,linux}.tgz.
#
# antechamber shells out to `charmmgen` to produce CHARMM topologies, but modern
# AmberTools dropped it and conda-forge's `ambertools` package has never contained it,
# so ACPYPE has to supply its own. The binaries shipped until now were built around
# 2016-2017: the macOS one was x86_64 against the 10.12 SDK (so Apple Silicon needed
# Rosetta 2) and the Linux one still declared a GLIBC 2.3 floor.
#
# The source is still maintained in Amber-MD/AmberClassic, where charmmgen is a unity
# build needing only charmmgen.c, eprintf.c and libm.
#
#   macOS  universal arm64 + x86_64, so both Apple Silicon and Intel run it natively
#   Linux  x86_64, built on Debian 12; needs at most GLIBC_2.34, inside the
#          manylinux_2_35 floor the wheels promise
#
# Both were verified to produce byte-identical .rtf/.prm/.inp output to the binaries
# they replace.
#
# Usage: scripts/build_charmmgen.sh [macos|linux|all] [output-dir]

set -euo pipefail

target="${1:-all}"
outdir="${2:-$PWD}"
repo="https://github.com/Amber-MD/AmberClassic.git"
min_macos="11.0"

# AmberClassic renamed AMBERHOME to AMBERCLASSICHOME. ACPYPE sets AMBERHOME and points
# it at the vendored bundle, so prefer that and keep the upstream name as a fallback.
read -r -d '' patch_py <<'PYEOF' || true
import sys, pathlib
p = pathlib.Path(sys.argv[1])
s = p.read_text()
old = '''    amberhome = (char *) getenv("AMBERCLASSICHOME");
    if (amberhome == NULL) {
        fprintf(stdout, "AMBERCLASSICHOME is not set!\\n");
        exit(1);
    }'''
new = '''    amberhome = (char *) getenv("AMBERHOME");
    if (amberhome == NULL)
        amberhome = (char *) getenv("AMBERCLASSICHOME");
    if (amberhome == NULL) {
        fprintf(stdout, "AMBERHOME is not set!\\n");
        exit(1);
    }'''
if old not in s:
    sys.exit("charmmgen.c does not match the expected AMBERCLASSICHOME block; upstream changed")
p.write_text(s.replace(old, new, 1))
PYEOF

build_macos() {
    if [[ "$(uname)" != "Darwin" ]]; then
        echo ">>> Skipping macOS build: not running on macOS" >&2
        return 0
    fi
    local workdir
    workdir="$(mktemp -d)"
    trap 'rm -rf "$workdir"' RETURN

    echo ">>> macOS: fetching AmberClassic antechamber source"
    git clone --depth 1 --filter=blob:none --sparse "$repo" "$workdir/AmberClassic" >/dev/null 2>&1
    git -C "$workdir/AmberClassic" sparse-checkout set src/antechamber >/dev/null 2>&1
    local src="$workdir/AmberClassic/src/antechamber"

    echo ">>> macOS: patching charmmgen to honour AMBERHOME"
    python3 -c "$patch_py" "$src/charmmgen.c"

    echo ">>> macOS: building universal (arm64 + x86_64)"
    local stage_bin="$workdir/stage/src/acpype/amber_macos/bin"
    mkdir -p "$stage_bin"
    (
        cd "$src"
        clang -arch arm64 -arch x86_64 -O2 -mmacosx-version-min="$min_macos" \
            -o "$stage_bin/charmmgen" charmmgen.c eprintf.c -lm
    )
    strip -S "$stage_bin/charmmgen"
    codesign -f -s - "$stage_bin/charmmgen"

    tar czf "$outdir/charmmgen_macos.tgz" -C "$workdir/stage" src/acpype/amber_macos/bin/charmmgen
    file "$stage_bin/charmmgen"
    echo ">>> Wrote $outdir/charmmgen_macos.tgz"
}

build_linux() {
    echo ">>> Linux: building in a linux/amd64 container"
    local image="acpype-charmmgen-builder"
    docker build --platform linux/amd64 -t "$image" - >/dev/null <<'EOF'
FROM debian:12
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc libc6-dev git ca-certificates python3 \
    && rm -rf /var/lib/apt/lists/*
EOF

    local workdir
    workdir="$(mktemp -d)"
    trap 'rm -rf "$workdir"' RETURN
    mkdir -p "$workdir/src/acpype/amber_linux/bin"
    printf '%s' "$patch_py" > "$workdir/patch.py"

    docker run --rm --platform linux/amd64 -v "$workdir":/out "$image" bash -c '
        set -e
        git clone --depth 1 --filter=blob:none --sparse '"$repo"' /src >/dev/null 2>&1
        git -C /src sparse-checkout set src/antechamber >/dev/null 2>&1
        cd /src/src/antechamber
        python3 /out/patch.py charmmgen.c
        gcc -O2 -o /out/src/acpype/amber_linux/bin/charmmgen charmmgen.c eprintf.c -lm
        strip /out/src/acpype/amber_linux/bin/charmmgen
    '

    rm -f "$workdir/patch.py"
    tar czf "$outdir/charmmgen_linux.tgz" -C "$workdir" src/acpype/amber_linux/bin/charmmgen
    echo ">>> Wrote $outdir/charmmgen_linux.tgz"
}

case "$target" in
macos) build_macos ;;
linux) build_linux ;;
all)
    build_macos
    build_linux
    ;;
*)
    echo "usage: $0 [macos|linux|all] [output-dir]" >&2
    exit 1
    ;;
esac
