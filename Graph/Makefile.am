noinst_LIBRARIES = libgraph.a
noinst_PROGRAMS = adjtodot

libgraph_a_CPPFLAGS = \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/ExpPaired

libgraph_a_SOURCES = \
	ContigGraph.cpp ContigGraph.h \
	DirectedGraph.cpp DirectedGraph.h DirectedGraphImpl.h

adjtodot_CPPFLAGS = -I$(top_srcdir)/Common
adjtodot_LDFLAGS = -L.
adjtodot_LDADD = libgraph.a \
	$(top_builddir)/Common/libcommon.a
adjtodot_SOURCES = adjtodot.cc
