#!/usr/bin/env bash

# Clean up old stuff
rm -f DAssembler

# Make new stuff
mkdir .deps
g++-4.2 -DHAVE_CONFIG_H -I. -I..  -I.. -I../Assembly -I../Common -I../DataLayer -I../Graph  -Wall -Wextra -fopenmp -g -O2 -MT DAssembler-DAssembler.o -MD -MP -MF .deps/DAssembler-DAssembler.Tpo -c -o DAssembler-DAssembler.o DAssembler.cpp
mv -f .deps/DAssembler-DAssembler.Tpo .deps/DAssembler-DAssembler.Po
g++-4.2 -Wall -Wextra -fopenmp -g -O2   -o DAssembler DAssembler-DAssembler.o ../DataLayer/libdatalayer.a ../Common/libcommon.a ../Graph/libgraph.a -ldl -lm

# Remove intermediate files
rm -rf .deps
rm -f DAssembler-DAssembler.o
