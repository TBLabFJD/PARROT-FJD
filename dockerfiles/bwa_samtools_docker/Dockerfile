# Buid the image on top of another base image 

FROM biocontainers/samtools:v1.9-4-deb_cv1

USER root

# Commands to be executed 

RUN apt update
RUN apt upgrade -y
#RUN apt install git -y
RUN apt install make -y
RUN apt install libz-dev -y

WORKDIR /data
#RUN git clone https://github.com/lh3/bwa.git && cd bwa && make
COPY bwa /root/bwa/
WORKDIR /root/bwa/
RUN make 

WORKDIR /

ENV PATH=$PATH:/root/bwa
RUN chmod 777 -R /root


# Nextflow -with-report
RUN apt-get update && apt install -y procps g++ && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* 



# Add  the necessary files 

CMD /bin/bash

#docker build -t bwa_samtools -f Dockerfile .


