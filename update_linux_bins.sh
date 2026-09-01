#!/usr/bin/env bash

# Refresh the vendored Linux AmberTools bundle.
#
# Unlike the macOS counterpart this cannot run against a local conda environment on a
# developer's Mac, so it builds a throwaway linux/amd64 image holding AmberTools and
# vendors out of that. Run it from the repository root with Docker available.
#
# The list of shared libraries lives in scripts/vendor_amber.py, which derives it from
# the executables; see that file for what is deliberately left to the host system.

set -euo pipefail

amber_version="${AMBERTOOLS_VERSION:-26}"
destination="acpype/amber_linux"
image="acpype-ambertools-linux:${amber_version}"

docker build --platform linux/amd64 -t "$image" - <<EOF
FROM condaforge/miniforge3
RUN mamba create -y -p /opt/amber -c conda-forge ambertools=${amber_version} && mamba clean -a -y
RUN apt-get update && apt-get install -y --no-install-recommends binutils file \
    && rm -rf /var/lib/apt/lists/*
EOF

rm -rf "$destination"
docker run --rm --platform linux/amd64 \
    -v "$PWD":/repo \
    -w /repo \
    "$image" \
    python3 scripts/vendor_amber.py --source /opt/amber --dest "$destination"

# charmmgen is ACPYPE's own build, not part of AmberTools.
tar xvfz charmmgen_linux.tgz

pre-commit run -a

find "$destination" -type f | wc -l
du -sh "$destination"
