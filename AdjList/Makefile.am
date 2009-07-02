bin_PROGRAMS = AdjList

AdjList_CPPFLAGS = \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer \
	-I$(top_srcdir)/ExpPaired \
	-I$(top_srcdir)/Graph

AdjList_LDADD = \
	$(top_builddir)/Common/libcommon.a \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/ExpPaired/libexppaired.a \
	$(top_builddir)/Graph/libgraph.a

AdjList_SOURCES = \
	AdjList.cpp
