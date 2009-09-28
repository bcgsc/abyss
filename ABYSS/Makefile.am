bin_PROGRAMS = ABYSS

ABYSS_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Assembly \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

ABYSS_LDADD = \
	$(top_builddir)/Assembly/libassembly.a \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a

ABYSS_SOURCES = Abyss.cpp
