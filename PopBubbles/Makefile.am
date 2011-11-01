bin_PROGRAMS = PopBubbles

PopBubbles_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Align \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

PopBubbles_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

PopBubbles_LDADD = \
	$(top_builddir)/Graph/libgraph.a \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Align/libalign.a \
	$(top_builddir)/Common/libcommon.a \
	$(top_builddir)/dialign/libdialign.a

PopBubbles_SOURCES = PopBubbles.cpp
