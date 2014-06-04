noinst_LIBRARIES = libbloom.a
bin_PROGRAMS = abyss-bloom

libbloom_a_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

libbloom_a_SOURCES = \
	Bloom.h \
	BloomFilter.h \
	CascadingBloomFilter.h \
	CascadingBloomFilterWindow.h

libbloom_a_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

abyss_bloom_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

abyss_bloom_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

abyss_bloom_LDADD = \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Align/libalign.a \
	$(top_builddir)/Common/libcommon.a \
	$(top_builddir)/Bloom/libbloom.a

abyss_bloom_SOURCES = bloom.cc \
	Bloom.h \
	BloomFilter.h \
	CascadingBloomFilter.h \
	CascadingBloomFilterWindow.h