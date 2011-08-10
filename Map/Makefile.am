bin_PROGRAMS = abyss-index abyss-map

abyss_index_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Align \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer \
	-I$(top_srcdir)/FMIndex

abyss_index_LDADD = \
	$(top_builddir)/FMIndex/libfmindex.a \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a

abyss_index_SOURCES = index.cc

abyss_map_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Align \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer \
	-I$(top_srcdir)/FMIndex

abyss_map_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

abyss_map_LDADD = \
	$(top_builddir)/FMIndex/libfmindex.a \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a

abyss_map_SOURCES = map.cc
