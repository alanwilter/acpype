#!/usr/bin/env bash

# Manual release fallback. Normally .github/workflows/release.yml does all of this
# when a GitHub Release is published: it re-runs the checks, publishes to PyPI over
# OIDC and pushes the Docker images. Keep this for the cases CI cannot cover.
#
# `pcregrep` used to read the version here; Homebrew no longer packages it (only
# `pcre2grep`), so plain grep is used instead.

set -euo pipefail

version="$(grep -woE "[0-9]{4}\.[12]?[0-9]\.[123]?[0-9]" src/acpype/__init__.py | head -1)"

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
    # One wheel per platform: a combined wheel carries both AmberTools trees and
    # exceeds PyPI's 100 MB per-file limit.
    # Two platform wheels plus a slim sdist; a combined wheel would be ~129 MB and a
    # full sdist ~121 MB, both over PyPI's 100 MB per-file limit.
    uv run python scripts/build_dists.py --out-dir dist
    uv publish dist/*"$version"*
    rm -vfr dist/*"$version"*
}

function run_docker() {
    echo ">>> Creating docker images"
    docker buildx build --platform linux/amd64 -t acpype/acpype:latest -t acpype/acpype:"$version" .
    echo ">>> Pushing docker images"
    docker push acpype/acpype --all-tags
    docker image rm acpype/acpype:"$version"
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

if $do_pip; then
    run_pip
fi

if $do_doc; then
    run_docker
fi

if $do_all; then
    run_both
fi
