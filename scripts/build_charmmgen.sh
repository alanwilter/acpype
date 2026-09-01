#!/usr/bin/env bash

# Build the `charmmgen` binary ACPYPE ships in charmmgen_macos.tgz.
#
# Modern AmberTools dropped charmmgen, and conda-forge's `ambertools` package has
# never contained it, so ACPYPE carries its own build for legacy CHARMM output. The
# binary shipped until now was x86_64 built against the macOS 10.12 SDK, which meant
# Rosetta 2 was required on Apple Silicon even though the rest of the bundle is arm64.
#
# The source is still maintained in Amber-MD/AmberClassic, where charmmgen is a unity
# build needing only charmmgen.c, eprintf.c and libm. This script builds a universal
# (arm64 + x86_64) binary so both Apple Silicon and Intel Macs run it natively.
#
# Verified byte-identical output (.rtf/.prm/.inp) against the previous x86_64 binary.

set -euo pipefail

repo="https://github.com/Amber-MD/AmberClassic.git"
workdir="$(mktemp -d)"
outdir="${1:-$PWD}"
min_macos="11.0"

trap 'rm -rf "$workdir"' EXIT

if [[ "$(uname)" != "Darwin" ]]; then
    echo "This script builds the macOS binary and must run on macOS." >&2
    exit 1
fi

echo ">>> Fetching AmberClassic antechamber source"
git clone --depth 1 --filter=blob:none --sparse "$repo" "$workdir/AmberClassic" >/dev/null 2>&1
git -C "$workdir/AmberClassic" sparse-checkout set src/antechamber >/dev/null 2>&1

src="$workdir/AmberClassic/src/antechamber"

# AmberClassic renamed AMBERHOME to AMBERCLASSICHOME. ACPYPE sets AMBERHOME and points
# it at the vendored bundle, so prefer that and keep the upstream name as a fallback.
echo ">>> Patching charmmgen to honour AMBERHOME"
python3 - "$src/charmmgen.c" <<'PY'
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
PY

echo ">>> Building universal (arm64 + x86_64) charmmgen"
# The tarball is unpacked from the repository root, so it must carry the full path.
stage_bin="$workdir/stage/acpype/amber_macos/bin"
mkdir -p "$stage_bin"
(
    cd "$src"
    clang -arch arm64 -arch x86_64 -O2 -mmacosx-version-min="$min_macos" \
        -o "$stage_bin/charmmgen" charmmgen.c eprintf.c -lm
)
strip -S "$stage_bin/charmmgen"
codesign -f -s - "$stage_bin/charmmgen"

echo ">>> Packaging"
tar czf "$outdir/charmmgen_macos.tgz" -C "$workdir/stage" acpype/amber_macos/bin/charmmgen

file "$stage_bin/charmmgen"
echo ">>> Wrote $outdir/charmmgen_macos.tgz"
