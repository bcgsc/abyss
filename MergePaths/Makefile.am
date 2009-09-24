bin_PROGRAMS = MergePaths MergeContigs

MergePaths_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer \
	-I$(top_srcdir)/ExpPaired \
	-I$(top_srcdir)/Graph

MergePaths_LDADD = \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a \
	$(top_builddir)/ExpPaired/libexppaired.a \
	$(top_builddir)/Graph/libgraph.a

MergePaths_SOURCES = \
	MergePaths.cpp

MergeContigs_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

MergeContigs_LDADD = \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a

MergeContigs_SOURCES = MergeContigs.cpp
