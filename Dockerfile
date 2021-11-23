#To build:
#docker build /path/to/acpype/ --tag acpype:2021.11.14

FROM ubuntu:18.04

LABEL maintainer="alanwilter@gmail.com,lucianopkagami@hotmail.com"

# set environment variables, to avoid pyc files and flushing buffer
ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1

RUN apt-get update && apt-get install -y \
    python3 \
    openbabel \
    git \
    gcc \
    libgfortran3

COPY amber21-11_linux /home/amber21

COPY amber21-11_linux/dat /usr/local/dat

COPY acpype_lib/acpype.py /home/

RUN bash /home/amber21/amber.sh

RUN touch /root/.bashrc \
    && echo "export AMBERHOME='/home/amber21'\nexport ACHOME='/home/amber21/bin'\nexport LD_LIBRARY_PATH='/home/amber21/lib'\n" >> /root/.bashrc

RUN cd /home/ && ln -s $PWD/run_acpype.py /usr/local/bin/acpype

RUN cd /home/amber21/bin && ln -s $PWD/antechamber /usr/local/bin/antechamber

RUN cd /home/amber21/bin && ln -s $PWD/charmmgen /usr/local/bin/charmmgen

RUN cd /home/amber21/bin && ln -s $PWD/tleap /usr/local/bin/tleap

RUN cd /home/amber21/bin && ln -s $PWD/teLeap /usr/local/bin/teLeap

RUN cd /home/amber21/bin && ln -s $PWD/parmchk2 /usr/local/bin/parmchk2

ENV ACHOME="/home/amber21/bin"

ENV LD_LIBRARY_PATH="/home/amber21/lib"

ENV AMBERHOME="/home/amber21"
