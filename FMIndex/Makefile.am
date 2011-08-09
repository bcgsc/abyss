noinst_PROGRAMS = fmconstruct fmsearch
lib_LIBRARIES = libfmindex.a

fmconstruct_CPPFLAGS = -I$(top_srcdir)/Common
fmconstruct_LDADD = libfmindex.a
fmconstruct_SOURCES = construct.cc

fmsearch_CPPFLAGS = -I$(top_srcdir)/Common
fmsearch_LDADD = libfmindex.a
fmsearch_SOURCES = search.cc

libfmindex_a_CPPFLAGS = -I$(top_srcdir)/Common

libfmindex_a_SOURCES = \
	BitArrays.h \
	bit_array.cc bit_array.h \
	FMIndex.h \
	sais.h
