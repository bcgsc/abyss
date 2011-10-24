bin_PROGRAMS = ABYSS-P

ABYSS_P_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Assembly \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

ABYSS_P_LDADD = \
	$(top_builddir)/Assembly/libassembly.a \
	$(top_builddir)/Common/libcommon.a \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(MPI_LIBS)

ABYSS_P_SOURCES = \
	parallelAbyss.cpp \
	CommLayer.cpp CommLayer.h \
	NetworkSequenceCollection.cpp NetworkSequenceCollection.h \
	MessageBuffer.cpp MessageBuffer.h \
	Messages.cpp Messages.h 
