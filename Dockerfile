# To build:
# docker build --tag acpype:$(date +%Y.%m.%d) /path/to/acpype/

FROM ubuntu:20.04

LABEL maintainer="alanwilter@gmail.com,lucianopkagami@hotmail.com"

# set environment variables, to avoid pyc files and flushing buffer
ENV PYTHONDONTWRITEBYTECODE 1

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    python3 \
    openbabel \
    python3-openbabel \
    libgfortran5 \
    libarpack++2-dev \
    python3-ipython \
    && apt-get autoremove -y && apt-get autoclean -y && apt-get clean -y \
    && rm -rf /var/lib/apt/lists/*

COPY run_acpype.py /home
COPY acpype /home/acpype
RUN ln -s /home/run_acpype.py /usr/local/bin/acpype
