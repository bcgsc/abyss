bin_PROGRAMS = abyss-gapfiller2

abyss_gapfiller2_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

abyss_gapfiller2_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

abyss_gapfiller2_LDADD = \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Align/libalign.a \
	$(top_builddir)/Common/libcommon.a

abyss_gapfiller2_SOURCES = gapfiller2.cc \
	$(top_srcdir)/connectpairs/BloomFilter.h \
	$(top_srcdir)/connectpairs/CountingBloomFilter.h \
	$(top_srcdir)/connectpairs/DBGBloom.h \
	$(top_srcdir)/connectpairs/DBGBloomAlgorithms.h \
	$(top_srcdir)/connectpairs/connectpairs.h
