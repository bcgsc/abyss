bin_PROGRAMS = ABYSS

ABYSS_CPPFLAGS = -I$(top_srcdir)

if HAVE_SQLITE
ABYSS_LDADD = \
	$(top_builddir)/DataBase/libdb.a \
	$(SQLITE_LIBS) \
	$(top_builddir)/Assembly/libassembly.a \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a
else
ABYSS_LDADD = \
	$(top_builddir)/Assembly/libassembly.a \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a
endif

ABYSS_SOURCES = Abyss.cpp
