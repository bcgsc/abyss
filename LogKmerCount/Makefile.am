bin_PROGRAMS = logcounter

logcounter_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

logcounter_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

logcounter_LDADD = \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Align/libalign.a \
	$(top_builddir)/Common/libcommon.a

logcounter_SOURCES = logcounter.cc \
	BloomFilterBase.h \
	CountingBloomFilter.h \
	plc.h
