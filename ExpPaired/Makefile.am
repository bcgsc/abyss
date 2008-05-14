bin_PROGRAMS = ExpPaired

ExpPaired_CPPFLAGS = -I$(top_srcdir)/Common -I$(top_srcdir)/DataLayer \
	-I$(top_srcdir)/Parallel

ExpPaired_LDADD = \
	$(top_builddir)/Common/libcommon.a \
	$(top_builddir)/DataLayer/libdatalayer.a

ExpPaired_SOURCES = \
	ExpPaired.cpp ExpPaired.h \
	NodeScores.cpp NodeScores.h \
	ParentTree.cpp ParentTree.h

