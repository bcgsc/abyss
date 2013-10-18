bin_PROGRAMS = abyss-connectpairs

abyss_connectpairs_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

abyss_connectpairs_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

abyss_connectpairs_LDADD = \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a

abyss_connectpairs_SOURCES = connectpairs.cc \
	BloomFilter.h \
	CountingBloomFilter.h \
	DBGBloom.h \
	DBGBloomAlgorithms.h
