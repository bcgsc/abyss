bin_PROGRAMS = abyss-filtergraph

abyss_filtergraph_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Align \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer \
	-I$(top_srcdir)/Graph

abyss_filtergraph_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

abyss_filtergraph_LDADD = \
	$(top_builddir)/Graph/libgraph.a \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Align/libalign.a \
	$(top_builddir)/Common/libcommon.a

abyss_filtergraph_SOURCES = FilterGraph.cpp
