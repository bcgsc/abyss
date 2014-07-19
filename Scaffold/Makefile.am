bin_PROGRAMS = abyss-scaffold abyss-junction abyss-longseqdist
noinst_PROGRAMS = abyss-drawgraph

abyss_scaffold_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common

abyss_scaffold_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

if HAVE_SQLITE
abyss_scaffold_LDADD = \
	$(top_builddir)/DataBase/libdb.a \
	$(SQLITE_LIBS) \
	$(top_builddir)/Common/libcommon.a
else
abyss_scaffold_LDADD = \
	$(top_builddir)/Common/libcommon.a
endif

abyss_scaffold_SOURCES = scaffold.cc

abyss_longseqdist_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common

abyss_longseqdist_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

abyss_longseqdist_LDADD = \
	$(top_builddir)/Common/libcommon.a

abyss_longseqdist_SOURCES = longseqdist.cpp

abyss_junction_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common

abyss_junction_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

abyss_junction_LDADD = \
	$(top_builddir)/Common/libcommon.a

abyss_junction_SOURCES = junction.cc

abyss_drawgraph_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/Common
abyss_drawgraph_LDADD = $(top_builddir)/Common/libcommon.a
abyss_drawgraph_SOURCES = drawgraph.cc
