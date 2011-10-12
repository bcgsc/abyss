noinst_LIBRARIES = libfmindex.a
noinst_PROGRAMS = abyss-dawg

libfmindex_a_CPPFLAGS = -I$(top_srcdir)/Common

libfmindex_a_SOURCES = \
	BitArrays.h \
	bit_array.cc bit_array.h \
	DAWG.h \
	FMIndex.h \
	sais.hxx

abyss_dawg_SOURCES = abyss-dawg.cc
abyss_dawg_LDADD = libfmindex.a
abyss_dawg_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common
