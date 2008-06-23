bin_PROGRAMS = AdjList

AdjList_CPPFLAGS = -I$(top_srcdir)/Common -I$(top_srcdir)/DataLayer \
	-I$(top_srcdir)/Parallel

AdjList_LDFLAGS = -Wl,--start-group

AdjList_LDADD = \
	$(top_builddir)/Common/libcommon.a \
	$(top_builddir)/DataLayer/libdatalayer.a

AdjList_SOURCES = \
	AdjList.cpp
