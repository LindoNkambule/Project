FROM ubuntu:18.04
MAINTAINER Lindo Nkambule (lindonkambule116@gmail.com)
RUN echo 'debconf debconf/frontend select Noninteractive' | debconf-set-selections
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    python \
    r-base \
    wget \
    unzip
RUN apt-get update && \
    apt-get install -y openjdk-8-jdk openjdk-8-jre-headless
ENV java=/usr/lib/jvm/java-1.8.0-openjdk-amd64
RUN mkdir -p /usr/local/pipeline/Tools
RUN cd /usr/local/pipeline/Tools && \
    wget https://github.com/broadinstitute/picard/releases/download/2.18.26/picard.jar
RUN wget -O /tmp/gatk-4.1.0.0.zip https://github.com/broadinstitute/gatk/releases/download/4.1.0.0/gatk-4.1.0.0.zip && \
    unzip /tmp/gatk-4.1.0.0.zip -d /usr/local/pipeline/Tools
RUN rm -rf /tmp/*
