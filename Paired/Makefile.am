bin_PROGRAMS = PairedAbyss

PairedAbyss_CPPFLAGS = \
	-I$(top_srcdir)/Align \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer \
	-I$(top_srcdir)/Parallel \
	-I$(top_srcdir)/ExpPaired \
	-I$(top_srcdir)/Graph

PairedAbyss_LDFLAGS = -Wl,--start-group

PairedAbyss_LDADD = \
	$(top_builddir)/Align/libalign.a \
	$(top_builddir)/Common/libcommon.a \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/ExpPaired/libexppaired.a \
	$(top_builddir)/Graph/libgraph.a

PairedAbyss_SOURCES = \
	PairedAbyss.cpp PairedAbyss.h

