if PAIRED_DBG
bin_PROGRAMS = abyss-paired-dbg
noinst_LIBRARIES = libpaireddbg.a
endif

# abyss-paired-dbg

abyss_paired_dbg_CPPFLAGS = -DPAIRED_DBG -I$(top_srcdir)

libdb = $(top_builddir)/DataBase/libdb.a $(SQLITE_LIBS)

abyss_paired_dbg_LDADD = \
	$(top_builddir)/Assembly/libassembly.a \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a \
	$(libdb)

abyss_paired_dbg_SOURCES = \
	abyss-paired-dbg.cc \
	BranchRecord.h \
	KmerPair.cc KmerPair.h \
	Dinuc.h \
	PairedDBGAlgorithms.h \
	SequenceCollection.h

libpaireddbg_a_SOURCES = KmerPair.cc KmerPair.h

CXXFLAGS += -Wno-deprecated-declarations
