bin_PROGRAMS = abyss-sealer

dist_man1_MANS = abyss-sealer.1

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

# Convert the README.md to a man page using Pandoc
abyss-sealer.1: README.md
	-pandoc -s -o $@ $<
