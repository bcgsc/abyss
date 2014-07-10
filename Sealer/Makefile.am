bin_PROGRAMS = abyss-sealer

abyss_sealer_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

abyss_sealer_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

abyss_sealer_LDADD = \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Align/libalign.a \
	$(top_builddir)/Common/libcommon.a

abyss_sealer_SOURCES = sealer.cc \
	$(top_srcdir)/Bloom/BloomFilter.h \
	$(top_srcdir)/connectpairs/CountingBloomFilter.h \
	$(top_srcdir)/Konnector/DBGBloom.h \
	$(top_srcdir)/Konnector/DBGBloomAlgorithms.h \
	$(top_srcdir)/Konnector/konnector.h
