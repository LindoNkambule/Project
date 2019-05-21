FROM ubuntu:18.04
MAINTAINER Lindo Nkambule (lindonkambule116@gmail.com)
RUN echo 'debconf debconf/frontend select Noninteractive' | debconf-set-selections
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    apt-utils  \
    ant \
    automake \
    autoconf \
    bash \
    build-essential \
    bzip2 \
    c++11 \
    c++17 \
    ca-certificates \
    cmake \
    g++ \
    git \
    gcc \
    gzip \
    less \
    libbz2-dev \
    libcurl4-gnutls-dev \
    liblzma-dev \
    libconfig-dev \
    libgsl-dev \
    libncurses-dev \
    libperl-dev \
    make \
    python python-dev python-yaml ncurses-dev zlib1g-dev python-numpy python-pip \
    r-base \
    sudo \
    wget \
    perl \
    software-properties-common \
    pkg-config \
    unzip \
    ncurses-dev \
    libx11-dev libxpm-dev libxft-dev libxext-dev \
    vim \
    xvfb \
    zlib1g-dev
RUN apt-get update && \
    apt-get install -y openjdk-8-jdk openjdk-8-jre-headless
RUN mkdir -p /usr/local/pipeline/Data
RUN mkdir -p /usr/local/pipeline/Scripts
RUN mkdir -p /usr/local/pipeline/Tools
COPY whole_pipeline.sh /usr/local/pipeline/Scripts
RUN wget -O /tmp/bwa-0.7.17.tar.bz2 https://liquidtelecom.dl.sourceforge.net/project/bio-bwa/bwa-0.7.17.tar.bz2 && \
    tar -xjvf /tmp/bwa-0.7.17.tar.bz2 -C /usr/local/pipeline/Tools && \
    chmod -R 777 /usr/local/pipeline/Tools && \
    cd /usr/local/pipeline/Tools/bwa-0.7.17 && make && \
    cp -v /usr/local/pipeline/Tools/bwa-0.7.17/bwa /usr/local/bin
RUN wget -O /tmp/samtools-1.9.tar.bz2 https://liquidtelecom.dl.sourceforge.net/project/samtools/samtools/1.9/samtools-1.9.tar.bz2 && \
    tar -xvjf /tmp/samtools-1.9.tar.bz2 -C /usr/local/pipeline/Tools && \
    chmod -R 777 /usr/local/pipeline/Tools && \
    cd /usr/local/pipeline/Tools/samtools-1.9 && make && \
    cp -v /usr/local/pipeline/Tools/samtools-1.9/samtools /usr/local/bin
RUN cd /tmp && \
    git clone git://github.com/pezmaster31/bamtools.git && \
    chmod -R 777 bamtools && \
    cd bamtools && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    make install
RUN cd /usr/local/pipeline/Tools && \
    git clone git://github.com/samtools/htslib.git && \
    git clone git://github.com/samtools/bcftools.git && \
    cd bcftools && \
    autoheader && autoconf && ./configure --enable-libgsl --enable-perl-filters && \
    make && \
    ln -s /usr/local/pipeline/Tools/bcftools/bcftools /usr/local/bin/
RUN cd /tmp && \
    git clone --recursive git://github.com/ekg/freebayes.git && \
    chmod -R 777 freebayes/ && \
    cd freebayes/ && \
    make && \
    sudo make install
RUN cd /usr/local/pipeline/Tools && \
    wget https://github.com/broadinstitute/picard/releases/download/2.18.26/picard.jar
RUN wget -O /tmp/gatk-4.1.0.0.zip https://github.com/broadinstitute/gatk/releases/download/4.1.0.0/gatk-4.1.0.0.zip && \
    unzip /tmp/gatk-4.1.0.0.zip -d /usr/local/pipeline/Tools
RUN rm -rf /tmp/*
