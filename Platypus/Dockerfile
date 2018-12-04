FROM ubuntu:16.04
MAINTAINER Lindokuhle Nkambule lindonkambule116@gmail.com
RUN apt-get update && \
    apt-get install -y --no-install-recommends\
    apt-utils \
    autoconf \
    bash \
    bzip2 \
    ca-certificates \
    gcc \
    git \
    less \
    libbz2-dev \
    libcurl4-gnutls-dev \
    liblzma-dev \
    make \
    perl \
    python python-dev \
    sudo \
    wget \
    zlib1g-dev
RUN cd /tmp && \
    wget https://github.com/cython/cython/archive/0.28.5.tar.gz && \
    chmod 777 0.28.5.tar.gz && \
    tar -xvzf 0.28.5.tar.gz && \
    cd cython-0.28.5 && \
    python setup.py install
RUN cd /tmp && \
    wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && \
    tar -xjvf htslib-1.9.tar.bz2 && \
    chmod -R 777 htslib-1.9 && \
    cd htslib-1.9 && \
    autoconf && \
    ./configure && \
    make && \
    make install
ENV C_INCLUDE_PATH /usr/local/include
ENV LIBRARY_PATH /usr/local/lib
ENV LD_LIBRARY_PATH /usr/local/lib
RUN cd /tmp && \
    git clone https://github.com/andyrimmer/Platypus.git && \
    chmod -R 777 Platypus && \
    cd Platypus && \
    make && \
    chmod -R 777 ./* && \
    cp -vrf ./bin/* /usr/local/bin && \
    ln -s .bin/* /usr/local/bin
RUN rm -rf /tmp/*
