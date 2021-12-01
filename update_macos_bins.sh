#!/usr/bin/env bash
# update macOS antechamber bins

source="/usr/local/Caskroom/miniconda/base"

files=(bin/charmmgen bin/bondtype bin/tleap bin/atomtype bin/wrapped_progs/bondtype bin/wrapped_progs/atomtype bin/wrapped_progs/antechamber bin/wrapped_progs/parmchk2 bin/wrapped_progs/am1bcc bin/antechamber bin/parmchk2 bin/sqm bin/am1bcc bin/teLeap LICENSE GNU_LGPL_v3 amber.sh lib/libsff_fortran.dylib lib/libcifparse.dylib)

for item in "${files[@]}"; do
    dd=$(diff -q "${source}"/"${item}" amber21-11_os/"${item}")
    if [ -n "${dd}" ]; then
        cp -fv "${source}"/"${item}" amber21-11_os/"${item}"
    fi
done
