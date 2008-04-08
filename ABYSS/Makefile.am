bin_PROGRAMS = ABYSS

ABYSS_CPPFLAGS = -I$(top_srcdir)/Common -I$(top_srcdir)/DataLayer \
	-I$(top_srcdir)/Parallel

ABYSS_LDFLAGS = -L$(top_builddir)/Common -L$(top_builddir)/DataLayer

ABYSS_LDADD = -lcommon -ldatalayer

ABYSS_SOURCES = \
	Abyss.cpp Abyss.h \
	Path.cpp Path.h \
	PathDriver.cpp PathDriver.h
