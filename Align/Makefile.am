noinst_LIBRARIES = libalign.a

libalign_a_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/Common

libalign_a_SOURCES = \
	alignGlobal.cc alignGlobal.h \
	dialign.cpp dialign.h dna_diag_prob.cc \
	smith_waterman.cpp smith_waterman.h Options.h

bin_PROGRAMS = abyss-align abyss-mergepairs

abyss_align_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

abyss_align_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

abyss_align_LDADD = $(builddir)/libalign.a \
	$(top_builddir)/dialign/libdialign.a \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a

abyss_align_SOURCES = align.cc

abyss_mergepairs_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

abyss_mergepairs_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

abyss_mergepairs_LDADD = $(builddir)/libalign.a \
	$(top_builddir)/dialign/libdialign.a \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a

abyss_mergepairs_SOURCES = mergepairs.cc

