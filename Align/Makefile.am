noinst_LIBRARIES = libalign.a

libalign_a_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/Common

libalign_a_SOURCES = \
	Aligner.cpp Aligner.h Options.h SAM.h \
	dialign.cpp dialign.h dna_diag_prob.cc \
	needleman_wunsch.cpp needleman_wunsch.h \
	smith_waterman.cpp smith_waterman.h

bin_PROGRAMS = abyss-align

abyss_align_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

abyss_align_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

abyss_align_LDADD = $(builddir)/libalign.a \
	$(top_builddir)/dialign/libdialign.a \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a

abyss_align_SOURCES = align.cc
