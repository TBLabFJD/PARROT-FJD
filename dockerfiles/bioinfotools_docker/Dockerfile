# Buid the image on top of another base image 

FROM ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive


# Commands to be executed 



# Install R, python....

# Install the dependencies necessary to add a new repository over HTTPS:
RUN apt-get update && apt-get upgrade -y && apt install -y dirmngr gnupg apt-transport-https ca-certificates software-properties-common

# Add the CRAN repository to your system sources’ list:
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'

RUN apt-get install -y --no-install-recommends \
  build-essential \
  r-base \
  python3.6 \
  python3-pip \
  python3-setuptools \
  python3-dev \
  bzip2 \
  g++ \
  libbz2-dev \
  libcurl4-openssl-dev \
  liblzma-dev \
  make \
  ncurses-dev \
  wget \
  zlib1g-dev \
  git

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

RUN rm /tmp/bcftools-1.15.tar.bz2 



# Install bedtools
RUN apt-get install -y bedtools

 

# Install bgzip and tabix
RUN apt-get install -y tabix 

RUN R -e "install.packages('optparse')"
RUN R -e "install.packages('openxlsx')"


# Nextflow -with-report
 RUN apt-get update && apt install -y procps g++ && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* 



# Install mosdepth
RUN wget https://github.com/brentp/mosdepth/releases/download/v0.3.3/mosdepth
RUN mv mosdepth /bin/
RUN chmod 777 /bin/mosdepth


# Install samtools
RUN apt-get update && apt-get install -y samtools

# Install Exomedepth
RUN R -e "install.packages('BiocManager')"
RUN R -e "install.packages('VGAM')"
RUN R -e "BiocManager::install('GenomicRanges')"
RUN R -e "BiocManager::install('ExomeDepth')"
RUN R -e "install.packages('stringr')"


# Install Panelcn.MOPS
RUN R -e "BiocManager::install('panelcn.mops')"


#Install bedr package (bedtools)
RUN R -e "install.packages('bedr')"
RUN R -e "install.packages('plyr')"



# Install CONVADING requirements
RUN apt-get install -y libstatistics-normality-perl

RUN pip install numpy
RUN pip install pandas

RUN ln -s /usr/bin/python3.8 /usr/bin/python



# Install Automap
WORKDIR /home/docker/
RUN git clone https://github.com/mquinodo/AutoMap

RUN chmod 777 -R /home/docker/

RUN apt-get install -y bc


RUN apt-get install -y gfortran
RUN R -e "BiocManager::install('ExomeDepth')"




# Add  the necessary files 

CMD /bin/bash

# docker build -t bioinfotools_tmp2 -f Dockerfile2 .
# docker tag 9ab562a534ef tblabfjd/bioinfotools:2.0.0
# docker tag 9ab562a534ef tblabfjd/bioinfotools:latest

