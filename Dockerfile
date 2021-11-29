# To build:
# docker build --tag acpype:$(date +%Y.%m.%d) /path/to/acpype/

FROM ubuntu:20.04

LABEL maintainer="alanwilter@gmail.com,lucianopkagami@hotmail.com"

# set environment variables, to avoid pyc files and flushing buffer
ENV PYTHONDONTWRITEBYTECODE 1

COPY amber21-11_linux /home/amber21

RUN bash /home/amber21/amber.sh

RUN touch /root/.bashrc \
    && echo "export AMBERHOME=/home/amber21" >> /root/.bashrc \
    && echo "export LD_LIBRARY_PATH=/home/amber21/lib" >> /root/.bashrc \
    && echo "export PATH=\$PATH:\$AMBERHOME/bin" >> /root/.bashrc

ENV LD_LIBRARY_PATH="/home/amber21/lib"

ENV AMBERHOME="/home/amber21"

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
COPY acpype_lib /home/acpype_lib
RUN ln -s /home/run_acpype.py /usr/local/bin/acpype
ENV PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/home/amber21/bin" 