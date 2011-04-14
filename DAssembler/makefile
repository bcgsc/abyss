CXX = g++44
CC = g++44
CXXFLAGS = -I. -I.. -I../DataLayer -I../Common -g -O2 -Wall -Wextra 

all:		DAssembler

RotatedRead.o:	RotatedRead.cpp RotatedRead.h

DAssembler.o:	DAssembler.cpp RotatedRead.h Rotation.h

DAssembler:	DAssembler.o RotatedRead.o ../DataLayer/libdatalayer.a ../Common/libcommon.a

test:
	./DAssembler test_set.fa
	
clean:
	rm -f *.o DAssembler

