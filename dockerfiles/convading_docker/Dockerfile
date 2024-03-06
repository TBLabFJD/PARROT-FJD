# Buid the image on top of another base image 

FROM python:2



# Commands to be executed 

RUN apt-get update

# Install samtools and bedtools
RUN apt-get install -y samtools
RUN apt-get install -y bedtools


RUN apt-get install -y libstatistics-normality-perl

RUN pip install numpy
RUN pip install pandas


# Add  the necessary files 

CMD /bin/bash

#docker build -t convading -f Dockerfile .


