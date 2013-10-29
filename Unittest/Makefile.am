UNIT_TESTS = common_stringutil
check_PROGRAMS = common_stringutil
common_stringutil_SOURCES = $(top_srcdir)/Common/test/StringUtilTest.cpp
common_stringutil_CPPFLAGS = -I$(top_srcdir)
common_stringutil_LDADD = $(GTEST_LIBS)

UNIT_TESTS += common_histogram
check_PROGRAMS += common_histogram
common_histogram_SOURCES = $(top_srcdir)/Common/test/HistogramTest.cpp
common_histogram_CPPFLAGS = -I$(top_srcdir)
common_histogram_LDADD = $(GTEST_LIBS)

UNIT_TESTS += common_bitutil
check_PROGRAMS += common_bitutil
common_bitutil_SOURCES = $(top_srcdir)/Common/test/BitUtilTest.cpp
common_bitutil_CPPFLAGS = -I$(top_srcdir)
common_bitutil_LDADD = $(GTEST_LIBS)

UNIT_TESTS += common_sequence
check_PROGRAMS += common_sequence
common_sequence_SOURCES = $(top_srcdir)/Common/test/Sequence.cc
common_sequence_CPPFLAGS = -I$(top_srcdir)
common_sequence_LDADD = $(top_builddir)/Common/libcommon.a $(GTEST_LIBS)

UNIT_TESTS += BloomFilter
check_PROGRAMS += BloomFilter
BloomFilter_SOURCES = $(top_srcdir)/connectpairs/test/BloomFilter.cc
BloomFilter_CPPFLAGS = -I$(top_srcdir)
BloomFilter_LDADD = $(top_builddir)/Common/libcommon.a $(GTEST_LIBS)

UNIT_TESTS += connectpairs_DBGBloomAlgorithms
check_PROGRAMS += connectpairs_DBGBloomAlgorithms
connectpairs_DBGBloomAlgorithms_SOURCES = $(top_srcdir)/connectpairs/test/DBGBloomAlgorithmsTest.cpp
connectpairs_DBGBloomAlgorithms_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/Common
connectpairs_DBGBloomAlgorithms_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)
connectpairs_DBGBloomAlgorithms_LDADD = $(top_builddir)/Common/libcommon.a $(GTEST_LIBS)

UNIT_TESTS += graph_ConstrainedBFSVisitor
check_PROGRAMS += graph_ConstrainedBFSVisitor
graph_ConstrainedBFSVisitor_SOURCES = $(top_srcdir)/Graph/test/ConstrainedBFSVisitorTest.cpp
graph_ConstrainedBFSVisitor_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/Common
graph_ConstrainedBFSVisitor_LDADD = $(top_builddir)/Common/libcommon.a $(GTEST_LIBS)

UNIT_TESTS += graph_BidirectionalBFS
check_PROGRAMS += graph_BidirectionalBFS
graph_BidirectionalBFS_SOURCES = $(top_srcdir)/Graph/test/BidirectionalBFSTest.cpp
graph_BidirectionalBFS_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/Common
graph_BidirectionalBFS_LDADD = $(top_builddir)/Common/libcommon.a $(GTEST_LIBS)

UNIT_TESTS += graph_AllPathsSearch
check_PROGRAMS += graph_AllPathsSearch
graph_AllPathsSearch_SOURCES = $(top_srcdir)/Graph/test/AllPathsSearchTest.cpp
graph_AllPathsSearch_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/Common
graph_AllPathsSearch_LDADD = $(top_builddir)/Common/libcommon.a $(GTEST_LIBS)

UNIT_TESTS += graph_HashGraph
check_PROGRAMS += graph_HashGraph
graph_HashGraph_SOURCES = $(top_srcdir)/Graph/test/HashGraphTest.cpp
graph_HashGraph_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/Common
graph_HashGraph_LDADD = $(top_builddir)/Common/libcommon.a $(GTEST_LIBS)

UNIT_TESTS += graph_ConstrainedBidiBFSVisitor
check_PROGRAMS += graph_ConstrainedBidiBFSVisitor
graph_ConstrainedBidiBFSVisitor_SOURCES = \
	$(top_srcdir)/Graph/test/ConstrainedBidiBFSVisitorTest.cpp
graph_ConstrainedBidiBFSVisitor_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/Common
graph_ConstrainedBidiBFSVisitor_LDADD = $(top_builddir)/Common/libcommon.a $(GTEST_LIBS)

TESTS = $(UNIT_TESTS)
