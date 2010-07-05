bin_PROGRAMS = DistanceEst

DistanceEst_CPPFLAGS = \
	-I$(top_srcdir)/Align \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/Graph

DistanceEst_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

DistanceEst_LDADD = \
	$(top_builddir)/Align/libalign.a \
	$(top_builddir)/Common/libcommon.a \
	$(top_builddir)/Graph/libgraph.a

DistanceEst_SOURCES = DistanceEst.cpp MLE.cpp MLE.h
