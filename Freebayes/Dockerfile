FROM ubuntu:16.04
MAINTAINER Lindokuhle Nkambule lindonkambule116@gmail.com
RUN apt-get update && \
    apt-get install -y --no-install-recommends\
    apt-utils \
    automake \
    autoconf \
    build-essential \
    bash \
    bzip2 \
    ca-certificates \
    cmake \
    c++17 \
    less \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libconfig-dev \
    libncurses-dev \
    make \
    git \
    g++ \
    gcc \
    perl \
    sudo \
    wget \
    zlib1g-dev
RUN cd /tmp && \
    git clone --recursive git://github.com/ekg/freebayes.git && \
    chmod -R 777 freebayes/ && \
    cd freebayes/ && \
    make && \
    sudo make install
RUN rm -rf /tmp/*
