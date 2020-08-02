#To build:
#docker build /path/to/acpype/ --tag acpype:2020.07.25.08.41
#use acpype.sh to run

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

COPY amber19-0_linux /home/amber19

COPY amber19-0_linux/dat /usr/local/dat

COPY acpype_lib/acpype.py /home/

RUN touch /root/.bashrc \
 && echo "export AMBERHOME='/home/amber19'\nexport ACHOME='/home/amber19/bin'\nexport LD_LIBRARY_PATH='/home/amber19/lib'\n" >> /root/.bashrc

RUN cd /home/ && ln -s $PWD/acpype.py /usr/local/bin/acpype

RUN cd /home/amber19/bin/to_be_dispatched && ln -s $PWD/antechamber /usr/local/bin/antechamber

RUN cd /home/amber19/bin/to_be_dispatched && ln -s $PWD/charmmgen /usr/local/bin/charmmgen

RUN cd /home/amber19/bin/ && ln -s $PWD/tleap /usr/local/bin/tleap

RUN cd /home/amber19/bin/ && ln -s $PWD/teLeap /usr/local/bin/teLeap

RUN cd /home/amber19/bin/to_be_dispatched && ln -s $PWD/parmchk2 /usr/local/bin/parmchk2

ENV ACHOME="/home/amber19/bin"

ENV LD_LIBRARY_PATH="/home/amber19/lib"

ENV AMBERHOME="/home/amber19"
