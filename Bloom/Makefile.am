bin_PROGRAMS = abyss-bloom

abyss_bloom_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer \
	-I$(top_srcdir)/vendor

abyss_bloom_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

abyss_bloom_LDADD = \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Align/libalign.a \
	$(top_builddir)/Common/libcommon.a

abyss_bloom_SOURCES = bloom.cc \
	Bloom.h \
	BloomFilter.h \
	BloomFilterWindow.h \
	ConcurrentBloomFilter.h \
	CascadingBloomFilter.h \
	CascadingBloomFilterWindow.h \
	RollingBloomDBGVisitor.h \
	HashAgnosticCascadingBloom.h

CXXFLAGS += -Wno-deprecated-declarations
