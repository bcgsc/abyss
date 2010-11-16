bin_PROGRAMS = SimpleGraph

SimpleGraph_CPPFLAGS = -pthread -I$(top_srcdir) \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/Graph

SimpleGraph_CXXFLAGS = -Wno-strict-aliasing

SimpleGraph_LDADD = \
	$(top_builddir)/Graph/libgraph.a \
	$(top_builddir)/Common/libcommon.a

SimpleGraph_LDFLAGS = -pthread

SimpleGraph_SOURCES = SimpleGraph.cpp
