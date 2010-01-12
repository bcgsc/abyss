bin_PROGRAMS = PopBubbles

PopBubbles_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

PopBubbles_LDADD = \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a

PopBubbles_SOURCES = PopBubbles.cpp
