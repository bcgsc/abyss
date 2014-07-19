bin_PROGRAMS = PathOverlap

PathOverlap_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

if HAVE_SQLITE
PathOverlap_LDADD = \
	$(top_builddir)/DataBase/libdb.a \
	$(SQLITE_LIBS) \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a
else
PathOverlap_LDADD = \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a
endif

PathOverlap_SOURCES = \
	PathOverlap.cpp
