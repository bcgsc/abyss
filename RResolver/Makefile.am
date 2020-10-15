bin_PROGRAMS = abyss-rresolver-short

abyss_rresolver_short_CXXFLAGS = \
	$(AM_CXXFLAGS) $(OPENMP_CXXFLAGS) \
	-I$(top_srcdir) -I$(top_srcdir)/Common -I$(top_srcdir)/DataLayer

abyss_rresolver_short_SOURCES = \
	RResolverShort.cpp \
	Contigs.cpp Contigs.h

LDFLAGS += -pthread

abyss_rresolver_short_LDADD = \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a \
	libralgorithmsshort.a

noinst_LIBRARIES = libralgorithmsshort.a

libralgorithmsshort_a_CXXFLAGS = \
  $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS) \
	-I$(top_srcdir) -I$(top_srcdir)/Common -I$(top_srcdir)/DataLayer

libralgorithmsshort_a_SOURCES = \
	RAlgorithmsShort.cpp RAlgorithmsShort.h \
	BloomFilters.cpp BloomFilters.h \
	Contigs.cpp Contigs.h \
	SequenceTree.cpp SequenceTree.h \
	RUtils.cpp RUtils.h

libralgorithmsshort_a_LIBADD = \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a

EXTRA_DIST = btllib
