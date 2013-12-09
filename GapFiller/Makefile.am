bin_PROGRAMS = abyss-gapfill

abyss_gapfill_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer \
	-I$(top_srcdir)/Align

abyss_gapfill_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

abyss_gapfill_LDADD = \
	$(top_builddir)/Align/libalign.a \
    $(top_builddir)/dialign/libdialign.a \
    $(top_builddir)/DataLayer/libdatalayer.a \
    $(top_builddir)/Common/libcommon.a

abyss_gapfill_SOURCES = gapfill.cpp gapfill.h
