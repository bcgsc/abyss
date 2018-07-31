bin_PROGRAMS = ABYSS

ABYSS_CPPFLAGS = -I$(top_srcdir)

ABYSS_LDADD = \
	$(top_builddir)/DataBase/libdb.a \
	$(SQLITE_LIBS) \
	$(top_builddir)/Assembly/libassembly.a \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a

ABYSS_SOURCES = abyss.cc
