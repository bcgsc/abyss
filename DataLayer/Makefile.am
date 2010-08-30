bin_PROGRAMS = abyss-tofastq
noinst_LIBRARIES = libdatalayer.a

abyss_tofastq_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common

abyss_tofastq_LDADD = libdatalayer.a \
	$(top_builddir)/Common/libcommon.a

abyss_tofastq_SOURCES = abyss-tofastq.cc

libdatalayer_a_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common

libdatalayer_a_SOURCES = \
	FastaReader.cpp FastaReader.h \
	FastaWriter.cpp FastaWriter.h \
	Options.h
