#!/usr/bin/env bash

# Create releases for pip or docker or both

set -euo pipefail

function usage() {
    echo "syntax: $0 < [-p, -d] | -a > to create a release for pip or docker or both"
    echo " -p : for pip, create wheel and upload to https://pypi.org/project/acpype/ (if you have permission)"
    echo " -d : for docker, create images and upload to https://hub.docker.com/u/acpype (if you have permission)"
    echo " -a : do both above"
    echo " -v : verbose mode, print all commands"
    echo " -h : prints this message"
    exit 1
}

function run_pip() {
    echo ">>> Creating pip package"
    python3 -m build
    python3 -m twine upload --repository testpypi dist/*"$today"* # TestPyPI
    python3 -m twine upload --repository acpype dist/*"$today"*   # official release
    rm -vfr dist/*"$today"*
}

function run_docker() {
    echo ">>> Creating docker images"
    docker build -t acpype/acpype:latest -t acpype/acpype:"$today" .

    echo ">>> Pushing docker images"
    docker push acpype/acpype --all-tags
    docker image rm acpype/acpype:"$today"
}

function run_both() {
    run_docker
    run_pip
}

do_pip=false
do_doc=false
do_all=false
verb=false
no_args=true

while getopts "adpvh" optionName; do
    case "$optionName" in
    a) do_all=true ;;
    d) do_doc=true ;;
    p) do_pip=true ;;
    v) verb=true ;;
    h) usage ;;
    ?) usage ;;
    *) usage ;;
    esac
    no_args=false
done

if "${no_args}"; then
    usage
elif $do_all && ($do_doc || $do_pip); then
    usage
fi

if ${verb}; then
    set -x
fi

today="$(date +%Y.%m.%d)"

if $do_pip; then
    run_pip
fi

if $do_doc; then
    run_docker
fi

if $do_all; then
    run_both
fi
