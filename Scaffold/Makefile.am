bin_PROGRAMS = Scaffold

Scaffold_CPPFLAGS = -I$(top_srcdir)/Common -I$(top_srcdir)/DataLayer \
	-I$(top_srcdir)/Parallel

Scaffold_LDFLAGS = -L$(top_builddir)/Common -L$(top_builddir)/DataLayer

Scaffold_LDADD = -lcommon -ldatalayer

Scaffold_SOURCES = Scaffold.cpp Scaffold.h
