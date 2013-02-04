bin_PROGRAMS = abyss-gapfill

abyss_gapfill_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

abyss_gapfill_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

abyss_gapfill_LDADD = \
	$(top_builddir)/Graph/libgraph.a \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a

abyss_gapfill_SOURCES = gapfill.cpp
