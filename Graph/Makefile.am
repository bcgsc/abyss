noinst_LIBRARIES = libgraph.a
noinst_PROGRAMS = adjtodot

libgraph_a_CPPFLAGS = -I$(top_srcdir)/Common

libgraph_a_SOURCES = \
	ContigGraph.h \
	ContigGraphAlgorithms.cpp ContigGraphAlgorithms.h \
	DirectedGraph.cpp DirectedGraph.h

adjtodot_CPPFLAGS = -I$(top_srcdir)/Common
adjtodot_LDFLAGS = -L.
adjtodot_LDADD = libgraph.a \
	$(top_builddir)/Common/libcommon.a
adjtodot_SOURCES = adjtodot.cc
