noinst_LIBRARIES = libfmindex.a
noinst_PROGRAMS = abyss-count abyss-dawg

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

abyss_count_SOURCES = count.cc
abyss_count_LDADD = libfmindex.a
abyss_count_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common
