bin_PROGRAMS = MergePaths MergeContigs PathConsensus

MergePaths_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

MergePaths_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

if HAVE_SQLITE
MergePaths_LDADD = \
	$(top_builddir)/DataBase/libdb.a \
	$(SQLITE_LIBS) \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a
else
MergePaths_LDADD = \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a
endif

MergePaths_SOURCES = \
	MergePaths.cpp

MergeContigs_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Align \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

if HAVE_SQLITE
MergeContigs_LDADD = \
	$(top_builddir)/DataBase/libdb.a \
	$(SQLITE_LIBS) \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Align/libalign.a \
	$(top_builddir)/Common/libcommon.a
else
MergeContigs_LDADD = \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Align/libalign.a \
	$(top_builddir)/Common/libcommon.a
endif

MergeContigs_SOURCES = MergeContigs.cpp

PathConsensus_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Align \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer \
	-I$(top_srcdir)/SimpleGraph

if HAVE_SQLITE
PathConsensus_LDADD = \
	$(top_builddir)/DataBase/libdb.a \
	$(SQLITE_LIBS) \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Align/libalign.a \
	$(top_builddir)/Common/libcommon.a \
	$(top_builddir)/dialign/libdialign.a
else
PathConsensus_LDADD = \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Align/libalign.a \
	$(top_builddir)/Common/libcommon.a \
	$(top_builddir)/dialign/libdialign.a
endif

PathConsensus_SOURCES = \
	PathConsensus.cpp
