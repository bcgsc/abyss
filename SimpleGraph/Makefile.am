bin_PROGRAMS = SimpleGraph

SimpleGraph_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common

SimpleGraph_CXXFLAGS = $(AM_CXXFLAGS) -Wno-strict-aliasing

if HAVE_SQLITE
SimpleGraph_LDADD = \
	$(top_builddir)/DataBase/libdb.a \
	$(SQLITE_LIBS) \
	$(top_builddir)/Common/libcommon.a \
	-lpthread
else
SimpleGraph_LDADD = \
	$(top_builddir)/Common/libcommon.a \
	-lpthread
endif

SimpleGraph_SOURCES = SimpleGraph.cpp
