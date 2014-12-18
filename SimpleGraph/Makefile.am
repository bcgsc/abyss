bin_PROGRAMS = SimpleGraph

SimpleGraph_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common

SimpleGraph_CXXFLAGS = $(AM_CXXFLAGS) -Wno-strict-aliasing

SimpleGraph_LDADD = \
	$(top_builddir)/DataBase/libdb.a \
	$(SQLITE_LIBS) \
	$(top_builddir)/Common/libcommon.a \
	-lpthread

SimpleGraph_SOURCES = SimpleGraph.cpp
