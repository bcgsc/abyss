bin_PROGRAMS = MergePaths MergeContigs PathConsensus

MergePaths_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

MergePaths_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

MergePaths_LDADD = \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a \
	$(top_builddir)/Graph/libgraph.a

MergePaths_SOURCES = \
	MergePaths.cpp

MergeContigs_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Align \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

MergeContigs_LDADD = \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Graph/libgraph.a \
	$(top_builddir)/Align/libalign.a \
	$(top_builddir)/Common/libcommon.a

MergeContigs_SOURCES = MergeContigs.cpp

PathConsensus_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Align \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer \
	-I$(top_srcdir)/SimpleGraph

PathConsensus_LDADD = \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Graph/libgraph.a \
	$(top_builddir)/Align/libalign.a \
	$(top_builddir)/Common/libcommon.a \
	$(top_builddir)/dialign/libdialign.a

PathConsensus_SOURCES = \
	PathConsensus.cpp
