bin_PROGRAMS = konnector

konnector_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

konnector_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

konnector_LDADD = \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Align/libalign.a \
	$(top_builddir)/Common/libcommon.a

konnector_SOURCES = konnector.cc \
	DBGBloom.h \
	DBGBloomAlgorithms.h \
	konnector.h

CXXFLAGS += -Wno-deprecated-declarations
