bin_PROGRAMS = abyss-rresolver-short
AM_CXXFLAGS += $(PTHREAD_CFLAGS) -Wno-deprecated-declarations

abyss_rresolver_short_CXXFLAGS = \
	$(AM_CXXFLAGS) $(OPENMP_CXXFLAGS) \
	-I$(top_srcdir) -I$(top_srcdir)/Common -I$(top_srcdir)/DataLayer -I$(top_srcdir)/RResolver/btllib/include

abyss_rresolver_short_SOURCES = \
	RResolverShort.cpp \
	Contigs.cpp Contigs.h

LDFLAGS += -pthread

abyss_rresolver_short_LDADD = \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a \
	libralgorithmsshort.a \
	-lbtllib

noinst_LIBRARIES = libralgorithmsshort.a

libralgorithmsshort_a_CXXFLAGS = \
  $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS) \
	-I$(top_srcdir) -I$(top_srcdir)/Common -I$(top_srcdir)/DataLayer -I$(top_srcdir)/RResolver/btllib/include

libralgorithmsshort_a_SOURCES = \
	RAlgorithmsShort.cpp RAlgorithmsShort.h \
	BloomFilters.cpp BloomFilters.h \
	Contigs.cpp Contigs.h \
	SequenceTree.cpp SequenceTree.h \
	RUtils.cpp RUtils.h

libralgorithmsshort_a_LIBADD = \
	$(wildcard $(top_builddir)/DataLayer/*.o) \
	$(wildcard $(top_builddir)/Common/*.o)

