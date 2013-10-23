bin_PROGRAMS = abyss-gc abyss-todot

abyss_todot_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/Common
abyss_todot_LDFLAGS = -L.
abyss_todot_LDADD = \
	$(top_builddir)/Common/libcommon.a
abyss_todot_SOURCES = todot.cc

abyss_gc_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/Common
abyss_gc_LDADD = $(top_builddir)/Common/libcommon.a
abyss_gc_SOURCES = gc.cc
