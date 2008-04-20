bin_PROGRAMS = Scaffold

Scaffold_CPPFLAGS = \
	-I$(top_srcdir)/Align \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer \
	-I$(top_srcdir)/Parallel

Scaffold_LDFLAGS = \
	-L$(top_builddir)/Align \
	-L$(top_builddir)/Common \
	-L$(top_builddir)/DataLayer

Scaffold_LDADD = -lalign -lcommon -ldatalayer

Scaffold_SOURCES = \
	Scaffold.cpp Scaffold.h \
	Stats.cpp Stats.h

