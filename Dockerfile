#To build:
#docker build /path/to/acpype/ --tag acpype:2020.07.25.08.41

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

RUN mkdir /home/amber19 && mkdir /home/test

COPY amber19-0_linux /home/amber19

COPY acpype_lib/acpype.py /home/

COPY test /home/test

RUN touch /root/.bashrc \
 && echo "export AMBERHOME='/home/amber19'\nexport ACHOME='/home/amber19/bin'\nexport LD_LIBRARY_PATH='/home/amber19/lib'\n" >> /root/.bashrc

RUN cd /home/ && ln -s $PWD/acpype.py /usr/local/bin/acpype

RUN cd /home/amber19/bin && ln -s $PWD/antechamber /usr/local/bin/antechamber

ENTRYPOINT ["acpype"]

CMD ["-i"]
