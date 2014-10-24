if PAIRED_DBG
bin_PROGRAMS = abyss-paired-dbg
endif

# abyss-paired-dbg

abyss_paired_dbg_CPPFLAGS = -DPAIRED_DBG -I$(top_srcdir)

if HAVE_SQLITE
libdb = $(top_builddir)/DataBase/libdb.a $(SQLITE_LIBS)
endif

abyss_paired_dbg_LDADD = \
	$(top_builddir)/Assembly/libassembly.a \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a \
	$(libdb)

abyss_paired_dbg_SOURCES = \
	abyss-paired-dbg.cc \
	BranchRecord.h \
	KmerPair.cc KmerPair.h \
	SequenceCollection.h
