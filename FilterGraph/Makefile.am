bin_PROGRAMS = abyss-filtergraph

abyss_filtergraph_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

abyss_filtergraph_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

abyss_filtergraph_LDADD = \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a

abyss_filtergraph_SOURCES = FilterGraph.cc
