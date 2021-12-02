#!/usr/bin/env bash

docker run -i -t -v "${PWD}":/wdir -w /wdir -u "$(id -u "${USER}"):$(id -g "${USER}")" acpype/acpype acpype "$@"
