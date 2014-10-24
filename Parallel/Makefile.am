bin_PROGRAMS = ABYSS-P
if PAIRED_DBG
bin_PROGRAMS += abyss-paired-dbg-mpi
endif

if HAVE_SQLITE
libdb = $(top_builddir)/DataBase/libdb.a $(SQLITE_LIBS)
endif

ABYSS_P_CPPFLAGS = -I$(top_srcdir)

ABYSS_P_LDADD = \
	$(top_builddir)/Assembly/libassembly.a \
	$(top_builddir)/Common/libcommon.a \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(libdb) \
	$(MPI_LIBS)

ABYSS_P_SOURCES = \
	parallelAbyss.cpp \
	CommLayer.cpp CommLayer.h \
	NetworkSequenceCollection.cpp NetworkSequenceCollection.h \
	MessageBuffer.cpp MessageBuffer.h \
	Messages.cpp Messages.h

abyss_paired_dbg_mpi_CPPFLAGS = $(ABYSS_P_CPPFLAGS) -DPAIRED_DBG

abyss_paired_dbg_mpi_LDADD = $(ABYSS_P_LDADD)

abyss_paired_dbg_mpi_SOURCES = $(ABYSS_P_SOURCES)
