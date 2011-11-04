noinst_LIBRARIES = libgraph.a
bin_PROGRAMS = abyss-todot

libgraph_a_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/Common

libgraph_a_SOURCES = \
	AdjIO.h \
	AsqgIO.h \
	ConstrainedSearch.cpp ConstrainedSearch.h \
	ContigGraph.h \
	ContigGraphAlgorithms.h \
	DepthFirstSearch.h \
	DirectedGraph.h \
	DotIO.h \
	FastaIO.h \
	GraphAlgorithms.h \
	GraphIO.h \
	GraphUtil.h \
	Options.h \
	PopBubbles.h \
	Properties.h \
	SAMIO.h

abyss_todot_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/Common
abyss_todot_LDFLAGS = -L.
abyss_todot_LDADD = libgraph.a \
	$(top_builddir)/Common/libcommon.a
abyss_todot_SOURCES = todot.cc
