bin_PROGRAMS = PathOverlap

PathOverlap_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

PathOverlap_LDADD = \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Graph/libgraph.a \
	$(top_builddir)/Common/libcommon.a

PathOverlap_SOURCES = \
	PathOverlap.cpp
