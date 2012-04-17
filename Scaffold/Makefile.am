bin_PROGRAMS = abyss-scaffold abyss-junction
noinst_PROGRAMS = abyss-drawgraph

abyss_scaffold_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common

abyss_scaffold_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

abyss_scaffold_LDADD = \
	$(top_builddir)/Common/libcommon.a

abyss_scaffold_SOURCES = scaffold.cc

abyss_junction_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common

abyss_junction_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

abyss_junction_LDADD = \
	$(top_builddir)/Common/libcommon.a

abyss_junction_SOURCES = junction.cc

abyss_drawgraph_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/Common
abyss_drawgraph_LDADD = $(top_builddir)/Common/libcommon.a
abyss_drawgraph_SOURCES = drawgraph.cc
