#!/bin/bash
# use docker pull lpkagami/acpype
# command: sh acpype.sh -i molecule.mol2

while [ "$#" -gt 0 ]; do
    case $1 in
        -i|--input) input="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

echo "Running Acpype Docker for $input file.."
docker run --name rm_tmp -v $PWD:/home/test lpkagami/acpype:latest acpype -i /home/test/$input
docker cp rm_tmp:/"${input%.*}".acpype $PWD
docker rm rm_tmp
echo "Done"
