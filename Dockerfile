FROM ubuntu:latest
MAINTAINER Shaun Jackman <sjackman@gmail.com>

RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
		libopenmpi1.6 make openmpi-bin ssh
ADD . /tmp/abyss
RUN apt-get install -y --no-install-recommends \
		automake g++ libboost-dev libopenmpi-dev libsparsehash-dev \
	&& cd /tmp/abyss \
	&& ./autogen.sh \
	&& mkdir build && cd build \
	&& ../configure --with-mpi=/usr/lib/openmpi \
	&& make install-strip \
	&& rm -rf /tmp/abyss \
	&& apt-get autoremove -y binutils \
		automake g++ libboost-dev libopenmpi-dev libsparsehash-dev
ENV SHELL=/bin/bash
ENTRYPOINT ["abyss-pe"]
CMD ["help"]
