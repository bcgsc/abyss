bin_PROGRAMS = abyss

abyss_CPPFLAGS = -I$(top_srcdir)/Common -I$(top_srcdir)/DataLayer \
	-I$(top_srcdir)/Parallel

abyss_LDFLAGS = -L$(top_builddir)/Common -L$(top_builddir)/DataLayer

abyss_LDADD = -lcommon -ldatalayer

abyss_SOURCES = \
	Abyss.cpp Abyss.h \
	Path.cpp Path.h \
	PathDriver.cpp PathDriver.h
