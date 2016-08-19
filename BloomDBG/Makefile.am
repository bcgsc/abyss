bin_PROGRAMS = abyss-bloom-dbg

abyss_bloom_dbg_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

abyss_bloom_dbg_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

abyss_bloom_dbg_LDADD = \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a

abyss_bloom_dbg_SOURCES = bloom-dbg.cc \
	bloom-dbg.h \
	MaskedKmer.h \
	SpacedSeed.h \
	HashAgnosticCascadingBloom.h \
	LightweightKmer.h \
	RollingBloomDBG.h \
	RollingHash.h \
	RollingHashIterator.h \
	$(top_srcdir)/lib/bloomfilter/BloomFilter.hpp \
	$(top_srcdir)/lib/rolling-hash/rolling.h
