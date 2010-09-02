noinst_LIBRARIES = libgraph.a
noinst_PROGRAMS = abyss-adjtodot

libgraph_a_CPPFLAGS = -I$(top_srcdir)/Common

libgraph_a_SOURCES = \
	AdjIO.h \
	ContigGraph.h \
	ContigGraphAlgorithms.h \
	DirectedGraph.cpp DirectedGraph.h \
	DotIO.h \
	Graph.h \
	MapGraph.h \
	SAMIO.h

abyss_adjtodot_CPPFLAGS = -I$(top_srcdir)/Common
abyss_adjtodot_LDFLAGS = -L.
abyss_adjtodot_LDADD = libgraph.a \
	$(top_builddir)/Common/libcommon.a
abyss_adjtodot_SOURCES = adjtodot.cc
