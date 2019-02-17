FROM ubuntu:16.04
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
    python-software-properties \
    pkg-config \
    unzip \
    ncurses-dev \
    libx11-dev libxpm-dev libxft-dev libxext-dev \
    vim \
    xvfb \
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
RUN mkdir -p /usr/local/pipeline/Data
RUN mkdir -p /usr/local/pipeline/Scripts
RUN mkdir -p /usr/local/pipeline/Tools
COPY whole_pipeline.sh /usr/local/pipeline/Scripts
RUN wget -O /tmp/fastqc_v0.11.7.zip http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip && \
    unzip /tmp/fastqc_v0.11.7.zip -d /usr/local/pipeline/Tools && \
    cd /usr/local/pipeline/Tools/FastQC && \
    chmod +x /usr/local/pipeline/Tools/FastQC/fastqc && \
    ln -s /usr/local/pipeline/Tools/FastQC/fastqc /usr/local/bin/fastqc
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
RUN cd /usr/local/pipeline/Tools && \
    git clone --recursive https://github.com/AstraZeneca-NGS/VarDictJava.git && \
    cd VarDictJava && \
    ./gradlew clean installDist && \
    ./gradlew distZip && \
    cp -r /usr/local/pipeline/Tools/VarDictJava/VarDict/* /usr/local/bin/
RUN cd /usr/local/pipeline/Tools && \
    wget https://github.com/broadinstitute/picard/releases/download/2.18.26/picard.jar
RUN wget -O /tmp/gatk-4.1.0.0.zip https://github.com/broadinstitute/gatk/releases/download/4.1.0.0/gatk-4.1.0.0.zip && \
    unzip /tmp/gatk-4.1.0.0.zip -d /usr/local/pipeline/Tools
RUN rm -rf /tmp/*
