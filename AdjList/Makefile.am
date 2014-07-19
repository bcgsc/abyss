bin_PROGRAMS = AdjList

AdjList_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

if HAVE_SQLITE
AdjList_LDADD = \
	$(top_builddir)/DataBase/libdb.a \
	$(SQLITE_LIBS) \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a
else
AdjList_LDADD = \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a
endif

AdjList_SOURCES = \
	AdjList.cpp
