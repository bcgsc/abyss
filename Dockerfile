FROM ubuntu:18.04
MAINTAINER Shaun Jackman <sjackman@gmail.com>

RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
		bsdmainutils libgomp1 make openmpi-bin ssh sudo \
	&& useradd -m -s /bin/bash abyss \
	&& echo 'abyss ALL=(ALL) NOPASSWD:ALL' >>/etc/sudoers
ADD . /tmp/abyss
RUN apt-get install -y --no-install-recommends \
		automake g++ libboost-dev libopenmpi-dev libsparsehash-dev \
	&& cd /tmp/abyss \
	&& ./autogen.sh \
	&& mkdir build && cd build \
	&& ../configure --with-mpi=/usr/lib/x86_64-linux-gnu/openmpi \
	&& make install-strip \
	&& rm -rf /tmp/abyss \
	&& apt-get autoremove -y binutils \
		automake g++ libboost-dev libopenmpi-dev libsparsehash-dev
USER abyss
WORKDIR /home/abyss
ENV SHELL=/bin/bash USER=abyss
ENTRYPOINT ["abyss-pe"]
CMD ["help"]
