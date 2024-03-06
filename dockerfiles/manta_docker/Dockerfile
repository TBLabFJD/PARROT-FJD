# Buid the image on top of another base image 

FROM python:2.7.18



# Commands to be executed 

RUN apt-get update && apt-get install -y \
	bzip2 \
	gcc \
	g++ \
	make \
	python \
	zlib1g-dev \
	wget



# Install samtools and bedtools
RUN wget https://github.com/Illumina/manta/releases/download/v1.6.0/manta-1.6.0.release_src.tar.bz2
RUN tar -xjf manta-1.6.0.release_src.tar.bz2
RUN mkdir build && cd build && ../manta-1.6.0.release_src/configure --jobs=4 --prefix=/manta && make -j4 install
ENV PATH="${PATH}:/manta/bin"


# Add  the necessary files 

CMD /bin/bash

#docker build -t manta -f Dockerfile .


