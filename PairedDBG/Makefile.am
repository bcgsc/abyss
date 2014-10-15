bin_PROGRAMS = abyss-paired-dbg
noinst_LIBRARIES = libpaireddbg.a

# abyss-paired-dbg

abyss_paired_dbg_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

if HAVE_SQLITE
libdb = $(top_builddir)/DataBase/libdb.a $(SQLITE_LIBS)
endif

abyss_paired_dbg_LDADD = libpaireddbg.a \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a \
	$(libdb)
abyss_paired_dbg_SOURCES = abyss-paired-dbg.cc

# libpaireddbg.a

libpaireddbg_a_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

libpaireddbg_a_SOURCES = \
	AssemblyAlgorithms.cc AssemblyAlgorithms.h \
	BranchGroup.cc BranchGroup.h \
	BranchRecord.cc BranchRecord.h \
	ISequenceCollection.h \
	KmerPair.cc KmerPair.h \
	KmerPairData.h \
	Options.cc Options.h \
	SequenceCollection.cc SequenceCollection.h

if HAVE_SQLITE
libpaireddbg_a_LIBADD = $(top_builddir)/Common/libcommon.a
endif
