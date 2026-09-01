# To build:
# docker build --tag acpype:$(date +%Y.%m.%d) /path/to/acpype/

FROM ubuntu:24.04

LABEL maintainer="alanwilter@gmail.com,lucianopkagami@hotmail.com"

# set environment variables, to avoid pyc files and flushing buffer
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

ARG DEBIAN_FRONTEND=noninteractive

# openbabel provides the `obabel` executable and python3-openbabel its bindings;
# the rest are the libraries the vendored AmberTools bundle deliberately takes
# from the host (see SYSTEM_LIBS_LINUX in scripts/vendor_amber.py).
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 openbabel python3-openbabel \
    libgfortran5 libstdc++6 libgomp1 libblas3 liblapack3 libcurl4 \
    && apt-get autoremove -y && apt-get clean -y \
    && rm -rf /var/lib/apt/lists/*

# run_acpype.py imports `acpype` from its own directory, so both must sit in /home.
COPY run_acpype.py /home/
COPY acpype /home/acpype
# Repository files carry restrictive owner-only permissions; the image runs as a
# non-root user, so make the launcher and the bundled binaries world-readable.
RUN chmod -R a+rX /home/acpype /home/run_acpype.py \
    && ln -s /home/run_acpype.py /usr/local/bin/acpype

RUN useradd --create-home --home-dir /work --shell /bin/bash acpype
USER acpype
WORKDIR /work
