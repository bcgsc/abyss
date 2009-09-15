bin_PROGRAMS = kmerprint

kmerprint_CPPFLAGS = \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

kmerprint_LDADD = \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a

kmerprint_SOURCES = kmerprint.cc
