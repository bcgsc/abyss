bin_PROGRAMS = ABYSS-P

ABYSS_P_CPPFLAGS = -I$(top_srcdir)

if HAVE_SQLITE
ABYSS_P_LDADD = \
	$(top_builddir)/DataBase/libdb.a \
	$(SQLITE_LIBS) \
	$(top_builddir)/Assembly/libassembly.a \
	$(top_builddir)/Common/libcommon.a \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(MPI_LIBS)
else
ABYSS_P_LDADD = \
	$(top_builddir)/Assembly/libassembly.a \
	$(top_builddir)/Common/libcommon.a \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(MPI_LIBS)
endif

ABYSS_P_SOURCES = \
	parallelAbyss.cpp \
	CommLayer.cpp CommLayer.h \
	NetworkSequenceCollection.cpp NetworkSequenceCollection.h \
	MessageBuffer.cpp MessageBuffer.h \
	Messages.cpp Messages.h
