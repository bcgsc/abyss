bin_PROGRAMS = PopBubbles

PopBubbles_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Align \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer \
	-I$(top_srcdir)/Graph

PopBubbles_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

PopBubbles_LDADD = \
	$(top_builddir)/Graph/libgraph.a \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Align/libalign.a \
	$(top_builddir)/Common/libcommon.a

PopBubbles_SOURCES = PopBubbles.cpp
