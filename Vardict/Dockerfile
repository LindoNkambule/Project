FROM ubuntu:16.04
MAINTAINER Lindo Nkambule (lindonkambule116@gmail.com)
RUN echo 'debconf debconf/frontend select Noninteractive' | debconf-set-selections
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    apt-utils \
    bzip2 \
    gcc \
    git \
    less \
    libbz2-dev \
    liblzma-dev \
    make \
    ncurses-dev \
    perl \
    r-base \
    sudo \
    vim \
    wget \
    zlib1g-dev
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y  software-properties-common && \
    add-apt-repository ppa:webupd8team/java && \
    apt-get update && \
    echo "oracle-java8-installer shared/accepted-oracle-license-v1-1 select true" | sudo debconf-set-selections && \
    apt-get install -y oracle-java8-installer && \
    apt-get install oracle-java8-set-default && \
    apt-get clean
RUN mkdir -p /usr/local/pipeline/Tools
RUN wget -O /tmp/samtools-1.9.tar.bz2 https://sourceforge.net/projects/samtools/files/samtools/1.9/samtools-1.9.tar.bz2 && \
    tar -xvjf /tmp/samtools-1.9.tar.bz2 -C /usr/local/pipeline/Tools && \
    chmod -R 777 /usr/local/pipeline/Tools && \
    cd /usr/local/pipeline/Tools/samtools-1.9 && make && \
    cp -v /usr/local/pipeline/Tools/samtools-1.9/samtools /usr/local/bin
RUN cd /usr/local/pipeline/Tools && \
    git clone --recursive https://github.com/AstraZeneca-NGS/VarDictJava.git && \
    cd VarDictJava && \
    ./gradlew clean installDist && \
    ./gradlew distZip && \
    cp -r /usr/local/pipeline/Tools/VarDictJava/VarDict/* /usr/local/bin/
RUN rm -rf /tmp/*
