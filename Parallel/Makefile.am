bin_PROGRAMS = ABYSS-P

ABYSS_P_CPPFLAGS = -I$(top_srcdir)/Common -I$(top_srcdir)/DataLayer \
	-I$(top_srcdir)/Parallel

ABYSS_P_LDFLAGS = -L$(top_builddir)/Common -L$(top_builddir)/DataLayer

ABYSS_P_LDADD = -lcommon -ldatalayer \
	-llammpi++ -lmpi -llam \
	-ldl -lpthread

ABYSS_P_SOURCES = \
	parallelAbyss.cpp parallelAbyss.h \
	CommLayer.cpp CommLayer.h \
	DistributedPhaseSpace.cpp DistributedPhaseSpace.h \
	NetworkDefs.h \
	NetworkSequenceCollection.cpp NetworkSequenceCollection.h
