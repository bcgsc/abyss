bin_PROGRAMS = DistanceEst DistanceEst-ssq

DistanceEst_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common

DistanceEst_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

DistanceEst_LDADD = \
	$(top_builddir)/DataBase/libdb.a \
	$(SQLITE_LIBS) \
	$(top_builddir)/Common/libcommon.a

DistanceEst_SOURCES = DistanceEst.cpp MLE.cpp MLE.h

DistanceEst_ssq_CPPFLAGS = $(DistanceEst_CPPFLAGS) \
	-D SAM_SEQ_QUAL=1

DistanceEst_ssq_CXXFLAGS = $(DistanceEst_CXXFLAGS)

DistanceEst_ssq_LDADD = $(DistanceEst_LDADD)

DistanceEst_ssq_SOURCES = $(DistanceEst_SOURCES)
