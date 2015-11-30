bin_PROGRAMS = abyss-bloom-dbg

abyss_bloom_dbg_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

abyss_bloom_dbg_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

abyss_bloom_dbg_LDADD = \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a

abyss_bloom_dbg_SOURCES = bloom-dbg.cc \
	HashAgnosticCascadingBloom.h \
	RollingBloomDBG.h \
	RollingHash.h \
	RollingHashIterator.h
