bin_PROGRAMS = ABYSS

ABYSS_CPPFLAGS = -I$(top_srcdir)/Common -I$(top_srcdir)/DataLayer \
	-I$(top_srcdir)/Parallel

ABYSS_LDFLAGS = -Wl,--start-group

ABYSS_LDADD = \
	$(top_builddir)/Common/libcommon.a \
	$(top_builddir)/DataLayer/libdatalayer.a

ABYSS_SOURCES = \
	Abyss.cpp \
	Path.cpp Path.h \
	PathDriver.cpp PathDriver.h
