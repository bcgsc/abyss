bin_PROGRAMS = abyss-rresolver-short

abyss_rresolver_short_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

abyss_rresolver_short_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

abyss_rresolver_short_LDADD = \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a

abyss_rresolver_short_SOURCES = \
	RAlgorithmsShort.cpp \
	BloomFilters.cpp \
	Contigs.cpp \
	SequenceTree.cpp \
	RUtils.cpp \
	RResolverShort.cpp
