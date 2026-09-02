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
    python3 python3-venv openbabel python3-openbabel \
    libgfortran5 libstdc++6 libgomp1 libblas3 liblapack3 libcurl4 \
    && apt-get autoremove -y && apt-get clean -y \
    && rm -rf /var/lib/apt/lists/*

# ACPYPE's CLI needs typer and rich, and Ubuntu only carries typer 0.9 while the
# project requires >= 0.20. A venv with system site packages keeps apt's
# python3-openbabel visible while pulling the two CLI dependencies from PyPI.
ENV VIRTUAL_ENV=/opt/venv
ENV PATH="${VIRTUAL_ENV}/bin:${PATH}"
RUN python3 -m venv --system-site-packages "${VIRTUAL_ENV}" \
    && "${VIRTUAL_ENV}"/bin/pip install --no-cache-dir --upgrade pip \
    && "${VIRTUAL_ENV}"/bin/pip install --no-cache-dir "typer>=0.20.0" "rich>=14.2.0"

# run_acpype.py imports `acpype` from its own directory, so both must sit in /home.
COPY run_acpype.py /home/
COPY src/acpype /home/acpype
# Repository files carry restrictive owner-only permissions; the image runs as a
# non-root user, so make the launcher and the bundled binaries world-readable.
RUN chmod -R a+rX /home/acpype /home/run_acpype.py \
    && ln -s /home/run_acpype.py /usr/local/bin/acpype

RUN useradd --create-home --home-dir /work --shell /bin/bash acpype
USER acpype
WORKDIR /work
