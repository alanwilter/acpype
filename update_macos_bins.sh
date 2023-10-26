#!/usr/bin/env bash

env_name="ambertools"
source="$HOME/mambaforge/envs/${env_name}"
destination="acpype/amber_macos"

if mamba env list | grep -q " $env_name "; then
    if [[ "$1" == "-f" ]]; then
        mamba env remove --name "$env_name" --yes
        mamba create -n "$env_name" ambertools --yes
    else
        echo "The '$env_name' environment already exists. Use '-f' to force re-run."
        exit 1
    fi
fi

files=(
    bin/bondtype
    bin/tleap
    bin/atomtype
    bin/wrapped_progs/bondtype
    bin/wrapped_progs/atomtype
    bin/wrapped_progs/antechamber
    bin/wrapped_progs/parmchk2
    bin/wrapped_progs/am1bcc
    bin/antechamber
    bin/parmchk2
    bin/sqm
    bin/am1bcc
    bin/teLeap
    LICENSE
    GNU_LGPL_v3
    amber.sh
    dat/antechamber
    dat/chamber
    dat/leap
    lib/libcifparse.dylib
    lib/libsff_fortran.dylib
    lib/libnetcdf.19.dylib
    lib/libmfhdf.0.dylib
    lib/libmfhdf.dylib
    lib/libdf.dylib
    lib/libdf.0.dylib
    lib/libhdf5_hl.310.dylib
    lib/libhdf5_hl.dylib
    lib/libhdf5.310.dylib
    lib/libhdf5.dylib
    lib/libcrypto.3.dylib
    lib/libzip.5.5.dylib
    lib/libzip.5.dylib
    lib/libsz.2.dylib
    lib/libsz.2.0.1.dylib
    lib/libblosc.1.dylib
    lib/libblosc.1.21.5.dylib
    lib/libzstd.1.dylib
    lib/libzstd.1.5.5.dylib
    lib/liblz4.1.dylib
    lib/liblz4.1.9.4.dylib
    lib/libsnappy.1.dylib
    lib/libsnappy.1.1.10.dylib
    lib/libjpeg.8.dylib
    lib/libjpeg.8.3.2.dylib
    lib/libgfortran.5.dylib
    lib/libquadmath.0.dylib
    lib/libgcc_s.1.1.dylib
    lib/libarpack.2.dylib
    lib/libarpack.2.1.0.dylib
    lib/libopenblas.0.dylib
    lib/libblas.3.dylib
    lib/liblapack.3.dylib
    lib/libopenblas.0.dylib
    lib/libomp.dylib
)

exclude=(
    --exclude=__pycache__
    --exclude=*.pyc
    --exclude=*.pyo
    --exclude=*.log
    --exclude=pixmaps/
)

rm -fr $destination

for item in "${files[@]}"; do
    source_path="${source}/${item}"
    destination_path="${destination}/${item}"
    pdir=$(dirname "$destination_path")

    mkdir -p "$pdir"
    if [ -d "$source_path" ]; then
        rsync -av "${exclude[@]}" "$source_path" "$pdir"
    else
        rsync -av "${exclude[@]}" "$source_path" "$destination_path"
    fi
done

tar xvfz charmmgen_macos.tgz

pre-commit run -a

tree -d $destination

find $destination | wc -l # 564 files, 16 dirs
# acpype/amber_linux
# ├── bin
# │   └── wrapped_progs
# ├── dat
# │   ├── antechamber
# │   ├── chamber
# │   └── leap
# │       ├── cmd
# │       │   └── oldff
# │       ├── lib
# │       │   └── oldff
# │       ├── parm
# │       └── prep
# │           ├── oldff
# │           └── protonated_nucleic
# └── lib
# amber.sh
# GNU_LGPL_v3
# LICENSE
