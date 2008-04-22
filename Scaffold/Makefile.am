bin_PROGRAMS = Scaffold

Scaffold_CPPFLAGS = \
	-I$(top_srcdir)/Align \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer \
	-I$(top_srcdir)/Parallel

Scaffold_LDADD = \
	$(top_builddir)/Align/libalign.a \
	$(top_builddir)/Common/libcommon.a \
	$(top_builddir)/DataLayer/libdatalayer.a

Scaffold_SOURCES = \
	Scaffold.cpp Scaffold.h \
	Stats.cpp Stats.h
