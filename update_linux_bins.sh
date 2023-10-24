#!/usr/bin/env bash

env_name="ambertools"
source="$HOME/mambaforge/envs/${env_name}"
destination="acpype/amber_linux"

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
    lib/libcifparse.so
    lib/libsff_fortran.so
    lib/libnetcdf.so.19
    lib/libmfhdf.so.0
    lib/libmfhdf.so.0.0.0
    lib/libdf.so.0
    lib/libdf.so.0.0.0*
    lib/libhdf5_hl.so.310
    lib/libhdf5_hl.so.310.0.2*
    lib/libhdf5.so.310
    lib/libhdf5.so.310.2.0
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

tar xvfz charmmgen.tgz

tree -d $destination

find $destination | wc -l # 551 files, 15 dirs
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
