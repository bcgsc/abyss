bin_PROGRAMS = abyss-fac abyss-tofastq
noinst_LIBRARIES = libdatalayer.a

abyss_fac_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common

abyss_fac_LDADD = libdatalayer.a \
	$(top_builddir)/Common/libcommon.a

abyss_fac_SOURCES = fac.cc

abyss_tofastq_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common

abyss_tofastq_LDADD = libdatalayer.a \
	$(top_builddir)/Common/libcommon.a

abyss_tofastq_SOURCES = abyss-tofastq.cc

libdatalayer_a_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common

libdatalayer_a_SOURCES = \
	FastaIndex.h \
	FastaInterleave.h \
	FastaReader.cpp FastaReader.h \
	FastaWriter.cpp FastaWriter.h \
	FastaConcat.h \
	Options.h
