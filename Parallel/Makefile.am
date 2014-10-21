bin_PROGRAMS = ABYSS-P

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
