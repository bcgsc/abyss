bin_PROGRAMS = DistanceEst

DistanceEst_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common

DistanceEst_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

DistanceEst_LDADD = \
	$(top_builddir)/Common/libcommon.a \
	$(top_builddir)/Graph/libgraph.a

DistanceEst_SOURCES = DistanceEst.cpp MLE.cpp MLE.h
