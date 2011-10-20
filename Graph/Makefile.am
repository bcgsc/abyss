noinst_LIBRARIES = libgraph.a
bin_PROGRAMS = abyss-adjtodot

libgraph_a_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/Common

libgraph_a_SOURCES = \
	AdjIO.h \
	AsqgIO.h \
	ConstrainedSearch.cpp ConstrainedSearch.h \
	ContigGraph.h \
	ContigGraphAlgorithms.h \
	DepthFirstSearch.h \
	DirectedGraph.cpp DirectedGraph.h \
	DotIO.h \
	FastaIO.h \
	Graph.h \
	GraphAlgorithms.h \
	GraphIO.h \
	GraphUtil.h \
	MapGraph.h \
	Options.h \
	PopBubbles.h \
	Properties.h \
	SAMIO.h

abyss_adjtodot_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/Common
abyss_adjtodot_LDFLAGS = -L.
abyss_adjtodot_LDADD = libgraph.a \
	$(top_builddir)/Common/libcommon.a
abyss_adjtodot_SOURCES = adjtodot.cc
