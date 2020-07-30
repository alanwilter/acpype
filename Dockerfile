#To build:
#docker build /path/to/acpype/ --tag acpype:2020.07.25.08.41

FROM debian:buster-slim

LABEL maintainer="alanwilter@gmail.com,lucianopkagami@hotmail.com"

# set environment variables, to avoid pyc files and flushing buffer
ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1

RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    git

RUN pip3 install git+https://github.com/alanwilter/acpype.git
