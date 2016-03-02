FROM ubuntu:latest
MAINTAINER Shaun Jackman <sjackman@gmail.com>

ADD . /root/abyss
RUN apt-get update && \
	apt-get install -y \
		automake \
		g++ \
		libboost-dev \
		libopenmpi-dev \
		libsparsehash-dev \
		libsqlite3-dev \
		make
RUN cd /root/abyss \
	&& ./autogen.sh \
	&& mkdir build && cd build \
	&& ../configure --with-mpi=/usr/lib/openmpi \
	&& make install-strip \
	&& rm -rf /root/abyss
ENTRYPOINT ["abyss-pe"]
CMD ["help"]
