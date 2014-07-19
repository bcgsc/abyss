bin_PROGRAMS = abyss-index abyss-map abyss-map-ssq abyss-overlap

abyss_index_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer \
	-I$(top_srcdir)/FMIndex

abyss_index_LDADD = \
	$(top_builddir)/FMIndex/libfmindex.a \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a

abyss_index_SOURCES = index.cc

abyss_map_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer \
	-I$(top_srcdir)/FMIndex

abyss_map_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

if HAVE_SQLITE
abyss_map_LDADD = \
	$(top_builddir)/FMIndex/libfmindex.a \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a \
	$(top_builddir)/DataBase/libdb.a \
	$(SQLITE_LIBS)
else
abyss_map_LDADD = \
	$(top_builddir)/FMIndex/libfmindex.a \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a
endif

abyss_map_SOURCES = map.cc

abyss_map_ssq_CPPFLAGS = $(abyss_map_CPPFLAGS) \
	-D SAM_SEQ_QUAL=1

abyss_map_ssq_CXXFLAGS = $(abyss_map_CXXFLAGS)

abyss_map_ssq_LDADD = $(abyss_map_LDADD)

abyss_map_ssq_SOURCES = $(abyss_map_SOURCES)

abyss_overlap_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer \
	-I$(top_srcdir)/FMIndex

abyss_overlap_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

abyss_overlap_LDADD = \
	$(top_builddir)/FMIndex/libfmindex.a \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a

abyss_overlap_SOURCES = overlap.cc
