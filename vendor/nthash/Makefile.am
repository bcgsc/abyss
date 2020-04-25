bin_PROGRAMS = nttest

nttest_CPPFLAGS = -I$(top_srcdir)/lib

nttest_SOURCES = \
	lib/BloomFilter.hpp \
	lib/city.cc \
	lib/city.h \
	lib/murmur.hpp \
	lib/seedgen.cpp \
	lib/seqgen.hpp \
	lib/xxhash.c \
	lib/xxhash.h \
	nthash.hpp \
	ntHashIterator.hpp \
	nttest.cpp

dist_doc_DATA = \
	ChangeLog \
	CITATION.bib \
	LICENSE \
	README.md

EXTRA_DIST=autogen.sh

noinst_HEADERS = stHashIterator.hpp

SUBDIRS = \
	unittest \
	vendor
