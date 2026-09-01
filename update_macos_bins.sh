#!/usr/bin/env bash

# Refresh the vendored macOS AmberTools bundle from a conda environment.
#
# The list of shared libraries is no longer maintained here: scripts/vendor_amber.py
# derives it from the executables, so an AmberTools upgrade that renames a library
# (libnetcdf.19 -> libnetcdf.22 in AmberTools 26) is picked up automatically instead
# of silently producing a half-broken bundle.

set -euo pipefail

amber_version="${AMBERTOOLS_VERSION:-26}"
env_name="ambertools"
source="$HOME/mambaforge/envs/${env_name}"
destination="acpype/amber_macos"

if mamba env list | grep -q " $env_name "; then
    if [[ "${1:-}" == "-f" ]]; then
        mamba env remove --name "$env_name" --yes
        mamba create -n "$env_name" "ambertools=${amber_version}" --yes
    else
        echo "The '$env_name' environment already exists. Use '-f' to force re-run."
        exit 1
    fi
else
    mamba create -n "$env_name" "ambertools=${amber_version}" --yes
fi

python3 scripts/vendor_amber.py --source "$source" --dest "$destination"

# charmmgen is ACPYPE's own build, not part of AmberTools.
tar xvfz charmmgen_macos.tgz

# conda-forge's arm64 Mach-O files declare some rpaths twice, which dyld on
# macOS 13+ refuses to load. Must run after every re-vendoring.
python3 scripts/fix_macos_rpaths.py "$destination"

pre-commit run -a

find "$destination" -type f | wc -l
du -sh "$destination"
