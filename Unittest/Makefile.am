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

UNIT_TESTS += DBGBloomAlgorithms
check_PROGRAMS += DBGBloomAlgorithms
DBGBloomAlgorithms_SOURCES = $(top_srcdir)/connectpairs/test/DBGBloomAlgorithmsTest.cpp
DBGBloomAlgorithms_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/Common
DBGBloomAlgorithms_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)
DBGBloomAlgorithms_LDADD = $(top_builddir)/Common/libcommon.a $(GTEST_LIBS)

UNIT_TESTS += ConstrainedBFSVisitor
check_PROGRAMS += ConstrainedBFSVisitor
ConstrainedBFSVisitor_SOURCES = $(top_srcdir)/Graph/test/ConstrainedBFSVisitorTest.cpp
ConstrainedBFSVisitor_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/Common
ConstrainedBFSVisitor_LDADD = $(top_builddir)/Common/libcommon.a $(GTEST_LIBS)

UNIT_TESTS += BidirectionalBFS
check_PROGRAMS += BidirectionalBFS
BidirectionalBFS_SOURCES = $(top_srcdir)/Graph/test/BidirectionalBFSTest.cpp
BidirectionalBFS_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/Common
BidirectionalBFS_LDADD = $(top_builddir)/Common/libcommon.a $(GTEST_LIBS)

TESTS = $(UNIT_TESTS)
