bin_PROGRAMS = abyss-dbgfm

abyss_dbgfm_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

abyss_dbgfm_LDADD = -ldbgfm \
	$(top_builddir)/FMIndex/libfmindex.a \
	$(top_builddir)/Assembly/libassembly.a \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a

abyss_dbgfm_SOURCES = dbgfm.cc dbgfm.h DBGAlgorithms.h
