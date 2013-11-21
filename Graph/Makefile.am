bin_PROGRAMS = abyss-gc abyss-todot

abyss_todot_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/Common
abyss_todot_LDFLAGS = -L.
abyss_todot_LDADD = \
	$(top_builddir)/Common/libcommon.a
abyss_todot_SOURCES = todot.cc

EXTRA_DIST = \
	AdjIO.h \
	AsqgIO.h \
	Assemble.h \
	ConstrainedSearch.h \
	ContigGraph.h \
	ContigGraphAlgorithms.h \
	DepthFirstSearch.h \
	DirectedGraph.h \
	DistIO.h \
	DotIO.h \
	FastaIO.h \
	GraphAlgorithms.h \
	GraphIO.h \
	GraphUtil.h \
	Options.h \
	PopBubbles.h \
	Properties.h \
	SAMIO.h

abyss_gc_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/Common
abyss_gc_LDADD = $(top_builddir)/Common/libcommon.a
abyss_gc_SOURCES = gc.cc
