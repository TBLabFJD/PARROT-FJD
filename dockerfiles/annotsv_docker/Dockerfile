# Buid the image on top of another base image 

FROM r-base



# Commands to be executed 

# Install bcftools
RUN apt-get update && apt-get upgrade -y && apt-get install -y \
  bzip2 \
  g++ \
  libbz2-dev \
  libcurl4-openssl-dev \
  liblzma-dev \
  make \
  ncurses-dev \
  wget \
  zlib1g-dev \
  git \
  curl

ENV BCFTOOLS_INSTALL_DIR=/opt/bcftools
ENV BCFTOOLS_VERSION=1.15

WORKDIR /tmp
RUN wget https://github.com/samtools/bcftools/releases/download/$BCFTOOLS_VERSION/bcftools-$BCFTOOLS_VERSION.tar.bz2 && \
  tar --bzip2 -xf bcftools-$BCFTOOLS_VERSION.tar.bz2

WORKDIR /tmp/bcftools-$BCFTOOLS_VERSION
RUN make prefix=$BCFTOOLS_INSTALL_DIR && \
  make prefix=$BCFTOOLS_INSTALL_DIR install

WORKDIR /
RUN ln -s $BCFTOOLS_INSTALL_DIR/bin/bcftools /usr/bin/bcftools && \
  rm -rf /tmp/bcftools-$BCFTOOLS_VERSION


# Install bedtools
RUN apt-get install -y bedtools

ENV DEBIAN_FRONTEND=noninteractive 
RUN apt install tclsh -y 

# Nextflow -with-report
 RUN apt-get update && apt install -y procps g++ && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* 


CMD /bin/bash

#docker build -t annotsv -f Dockerfile .
