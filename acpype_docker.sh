#!/usr/bin/env bash

docker run -i -t -v "${PWD}":/wdir -w /wdir acpype/acpype acpype "$@"
