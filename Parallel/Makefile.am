bin_PROGRAMS = ABYSS-P

ABYSS_P_CPPFLAGS = -I$(top_srcdir)/Common -I$(top_srcdir)/DataLayer \
	-I$(top_srcdir)/Parallel

ABYSS_P_LDADD = \
	$(top_builddir)/Common/libcommon.a \
	$(top_builddir)/DataLayer/libdatalayer.a

ABYSS_P_SOURCES = \
	parallelAbyss.cpp parallelAbyss.h \
	CommLayer.cpp CommLayer.h \
	DistributedPhaseSpace.cpp DistributedPhaseSpace.h \
	NetworkDefs.h \
	NetworkSequenceCollection.cpp NetworkSequenceCollection.h \
	MessageBuffer.cpp MessageBuffer.h \
	Messages.cpp Messages.h 
