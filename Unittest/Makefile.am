GTEST_LIBS = $(top_builddir)/lib/gtest-1.7.0/libgtest_main.a
# -Wno-error is used here because there is no portable way
# to suppress warning: "argument unused during compilation: '-pthread'"
# for clang on OSX.
# See: http://stackoverflow.com/questions/17841140/os-x-clang-pthread
GTEST_CXXFLAGS = $(AM_CXXFLAGS) $(PTHREAD_CFLAGS) -Wno-error
GTEST_INCLUDES = -I$(top_srcdir) -I$(top_srcdir)/lib/gtest-1.7.0/include
GTEST_LDFLAGS = $(PTHREAD_LIBS)

check_PROGRAMS = common_stringutil
common_stringutil_SOURCES = Common/StringUtilTest.cpp
common_stringutil_CPPFLAGS = $(GTEST_INCLUDES)
common_stringutil_LDADD = $(GTEST_LIBS)
common_stringutil_CXXFLAGS = $(GTEST_CXXFLAGS)
common_stringutil_LDFLAGS = $(GTEST_LDFLAGS)

check_PROGRAMS += common_histogram
common_histogram_SOURCES = Common/HistogramTest.cpp
common_histogram_CPPFLAGS = $(GTEST_INCLUDES)
common_histogram_LDADD = $(GTEST_LIBS)
common_histogram_CXXFLAGS = $(GTEST_CXXFLAGS)
common_histogram_LDFLAGS = $(GTEST_LDFLAGS)

check_PROGRAMS += common_bitutil
common_bitutil_SOURCES = Common/BitUtilTest.cpp
common_bitutil_CPPFLAGS = $(GTEST_INCLUDES)
common_bitutil_LDADD = $(GTEST_LIBS)
common_bitutil_CXXFLAGS = $(GTEST_CXXFLAGS)
common_bitutil_LDFLAGS = $(GTEST_LDFLAGS)

check_PROGRAMS += common_kmer
common_kmer_SOURCES = Common/KmerTest.cpp
common_kmer_CPPFLAGS = $(GTEST_INCLUDES)
common_kmer_LDADD = $(top_builddir)/Common/libcommon.a $(GTEST_LIBS)
common_kmer_CXXFLAGS = $(GTEST_CXXFLAGS)
common_kmer_LDFLAGS = $(GTEST_LDFLAGS)

check_PROGRAMS += common_sequence
common_sequence_SOURCES = Common/Sequence.cc
common_sequence_CPPFLAGS = $(GTEST_INCLUDES)
common_sequence_LDADD = $(top_builddir)/Common/libcommon.a $(GTEST_LIBS)
common_sequence_CXXFLAGS = $(GTEST_CXXFLAGS)
common_sequence_LDFLAGS = $(GTEST_LDFLAGS)

check_PROGRAMS += common_KmerIterator
common_KmerIterator_SOURCES = Common/KmerIteratorTest.cpp
common_KmerIterator_CPPFLAGS = $(GTEST_INCLUDES)
common_KmerIterator_LDADD = $(top_builddir)/Common/libcommon.a $(GTEST_LIBS)
common_KmerIterator_CXXFLAGS = $(GTEST_CXXFLAGS)
common_KmerIterator_LDFLAGS = $(GTEST_LDFLAGS)

check_PROGRAMS += common_sam
common_sam_SOURCES = Common/SAM.cc
common_sam_CPPFLAGS = $(GTEST_INCLUDES)
common_sam_LDADD = $(top_builddir)/Common/libcommon.a $(GTEST_LIBS)
common_sam_CXXFLAGS = $(GTEST_CXXFLAGS)
common_sam_LDFLAGS = $(GTEST_LDFLAGS)

check_PROGRAMS += BloomFilter
BloomFilter_SOURCES = Konnector/BloomFilter.cc
BloomFilter_CPPFLAGS = $(GTEST_INCLUDES) -I$(top_srcdir)/Common
BloomFilter_LDADD = $(top_builddir)/Common/libcommon.a $(GTEST_LIBS)
BloomFilter_CXXFLAGS = $(GTEST_CXXFLAGS) $(OPENMP_CXXFLAGS)
BloomFilter_LDFLAGS = $(GTEST_LDFLAGS)

check_PROGRAMS += Konnector_DBGBloom
Konnector_DBGBloom_SOURCES = Konnector/DBGBloomTest.cpp
Konnector_DBGBloom_CPPFLAGS = $(GTEST_INCLUDES) -I$(top_srcdir)/Common
Konnector_DBGBloom_LDADD = $(top_builddir)/Common/libcommon.a $(GTEST_LIBS)
Konnector_DBGBloom_CXXFLAGS = $(GTEST_CXXFLAGS) $(OPENMP_CXXFLAGS)
Konnector_DBGBloom_LDFLAGS = $(GTEST_LDFLAGS)

check_PROGRAMS += Konnector_DBGBloomAlgorithms
Konnector_DBGBloomAlgorithms_SOURCES = Konnector/DBGBloomAlgorithmsTest.cpp
Konnector_DBGBloomAlgorithms_CPPFLAGS = $(GTEST_INCLUDES) -I$(top_srcdir)/Common
Konnector_DBGBloomAlgorithms_CXXFLAGS = $(GTEST_CXXFLAGS) $(OPENMP_CXXFLAGS)
Konnector_DBGBloomAlgorithms_LDADD = \
	$(top_builddir)/Common/libcommon.a \
	$(GTEST_LIBS)
Konnector_DBGBloomAlgorithms_LDFLAGS = $(GTEST_LDFLAGS)

check_PROGRAMS += graph_ConstrainedBFSVisitor
graph_ConstrainedBFSVisitor_SOURCES = Graph/ConstrainedBFSVisitorTest.cpp
graph_ConstrainedBFSVisitor_CPPFLAGS = $(GTEST_INCLUDES) -I$(top_srcdir)/Common
graph_ConstrainedBFSVisitor_LDADD = $(top_builddir)/Common/libcommon.a $(GTEST_LIBS)
graph_ConstrainedBFSVisitor_CXXFLAGS = $(GTEST_CXXFLAGS)
graph_ConstrainedBFSVisitor_LDFLAGS = $(GTEST_LDFLAGS)

check_PROGRAMS += graph_BidirectionalBFS
graph_BidirectionalBFS_SOURCES = Graph/BidirectionalBFSTest.cpp
graph_BidirectionalBFS_CPPFLAGS = $(GTEST_INCLUDES) -I$(top_srcdir)/Common
graph_BidirectionalBFS_LDADD = $(top_builddir)/Common/libcommon.a $(GTEST_LIBS)
graph_BidirectionalBFS_CXXFLAGS = $(GTEST_CXXFLAGS)
graph_BidirectionalBFS_LDFLAGS = $(GTEST_LDFLAGS)

check_PROGRAMS += graph_AllPathsSearch
graph_AllPathsSearch_SOURCES = Graph/AllPathsSearchTest.cpp
graph_AllPathsSearch_CPPFLAGS = $(GTEST_INCLUDES) -I$(top_srcdir)/Common
graph_AllPathsSearch_LDADD = $(top_builddir)/Common/libcommon.a $(GTEST_LIBS)
graph_AllPathsSearch_CXXFLAGS = $(GTEST_CXXFLAGS)
graph_AllPathsSearch_LDFLAGS = $(GTEST_LDFLAGS)

check_PROGRAMS += graph_HashGraph
graph_HashGraph_SOURCES = Graph/HashGraphTest.cpp
graph_HashGraph_CPPFLAGS = $(GTEST_INCLUDES) -I$(top_srcdir)/Common
graph_HashGraph_LDADD = $(top_builddir)/Common/libcommon.a $(GTEST_LIBS)
graph_HashGraph_CXXFLAGS = $(GTEST_CXXFLAGS)
graph_HashGraph_LDFLAGS = $(GTEST_LDFLAGS)

check_PROGRAMS += graph_ConstrainedBidiBFSVisitor
graph_ConstrainedBidiBFSVisitor_SOURCES = \
	Graph/ConstrainedBidiBFSVisitorTest.cpp
graph_ConstrainedBidiBFSVisitor_CPPFLAGS = $(GTEST_INCLUDES) -I$(top_srcdir)/Common
graph_ConstrainedBidiBFSVisitor_LDADD = $(top_builddir)/Common/libcommon.a $(GTEST_LIBS)
graph_ConstrainedBidiBFSVisitor_CXXFLAGS = $(GTEST_CXXFLAGS)
graph_ConstrainedBidiBFSVisitor_LDFLAGS = $(GTEST_LDFLAGS)

check_PROGRAMS += graph_ExtendPath
graph_ExtendPath_SOURCES = Graph/ExtendPathTest.cpp
graph_ExtendPath_CPPFLAGS = $(GTEST_INCLUDES) -I$(top_srcdir)/Common
graph_ExtendPath_LDADD = $(top_builddir)/Common/libcommon.a $(GTEST_LIBS)
graph_ExtendPath_CXXFLAGS = $(GTEST_CXXFLAGS)
graph_ExtendPath_LDFLAGS = $(GTEST_LDFLAGS)

check_PROGRAMS += Konnector_konnector
Konnector_konnector_SOURCES = \
	Konnector/konnectorTest.cpp
Konnector_konnector_CPPFLAGS = $(GTEST_INCLUDES) -I$(top_srcdir)/Common
Konnector_konnector_LDADD = \
	$(top_builddir)/Align/libalign.a \
	$(top_builddir)/Common/libcommon.a $(GTEST_LIBS)
Konnector_konnector_CXXFLAGS = $(GTEST_CXXFLAGS) $(OPENMP_CXXFLAGS)
Konnector_konnector_LDFLAGS = $(GTEST_LDFLAGS)

check_PROGRAMS += DBG_LoadAlgorithm
DBG_LoadAlgorithm_SOURCES = \
	DBG/LoadAlgorithmTest.cpp
DBG_LoadAlgorithm_CPPFLAGS = \
	$(GTEST_INCLUDES) \
	-I$(top_srcdir)/DataLayer \
	-I$(top_srcdir)/Common
DBG_LoadAlgorithm_LDADD = \
	$(top_builddir)/Assembly/libassembly.a \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a $(GTEST_LIBS)
DBG_LoadAlgorithm_CXXFLAGS = $(GTEST_CXXFLAGS) $(OPENMP_CXXFLAGS)
DBG_LoadAlgorithm_LDFLAGS = $(GTEST_LDFLAGS)

if PAIRED_DBG

check_PROGRAMS += PairedDBG_LoadAlgorithm
PairedDBG_LoadAlgorithm_SOURCES = \
	PairedDBG/LoadAlgorithmTest.cpp
PairedDBG_LoadAlgorithm_CPPFLAGS = \
	$(GTEST_INCLUDES) \
	-I$(top_srcdir)/DataLayer \
	-I$(top_srcdir)/Common
PairedDBG_LoadAlgorithm_LDADD = \
	$(top_builddir)/PairedDBG/libpaireddbg.a \
	$(top_builddir)/Assembly/libassembly.a \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a $(GTEST_LIBS)
PairedDBG_LoadAlgorithm_CXXFLAGS = $(GTEST_CXXFLAGS) $(OPENMP_CXXFLAGS)
PairedDBG_LoadAlgorithm_LDFLAGS = $(GTEST_LDFLAGS)

check_PROGRAMS += PairedDBG_KmerPair
PairedDBG_KmerPair_SOURCES = \
	PairedDBG/KmerPairTest.cc
PairedDBG_KmerPair_CPPFLAGS = \
	$(GTEST_INCLUDES) \
	-I$(top_srcdir)/DataLayer \
	-I$(top_srcdir)/Common
PairedDBG_KmerPair_LDADD = \
	$(top_builddir)/PairedDBG/libpaireddbg.a \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a $(GTEST_LIBS)
PairedDBG_KmerPair_CXXFLAGS = $(GTEST_CXXFLAGS) $(OPENMP_CXXFLAGS)
PairedDBG_KmerPair_LDFLAGS = $(GTEST_LDFLAGS)

check_PROGRAMS += PairedDBG_Dinuc
PairedDBG_Dinuc_SOURCES = \
	PairedDBG/DinucTest.cc
PairedDBG_Dinuc_CPPFLAGS = \
	$(GTEST_INCLUDES) \
	-I$(top_srcdir)/DataLayer \
	-I$(top_srcdir)/Common
PairedDBG_Dinuc_LDADD = \
	$(top_builddir)/PairedDBG/libpaireddbg.a \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a $(GTEST_LIBS)
PairedDBG_Dinuc_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)
PairedDBG_Dinuc_LDFLAGS = $(GTEST_LDFLAGS)

check_PROGRAMS += PairedDBG_BranchRecord
PairedDBG_BranchRecord_SOURCES = \
	PairedDBG/BranchRecordTest.cpp
PairedDBG_BranchRecord_CPPFLAGS = \
	$(GTEST_INCLUDES) \
	-I$(top_srcdir)/DataLayer \
	-I$(top_srcdir)/Common
PairedDBG_BranchRecord_LDADD = \
	$(top_builddir)/PairedDBG/libpaireddbg.a \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a $(GTEST_LIBS)
PairedDBG_BranchRecord_CXXFLAGS = $(GTEST_CXXFLAGS) $(OPENMP_CXXFLAGS)
PairedDBG_BranchRecord_LDFLAGS = $(GTEST_LDFLAGS)

endif # PAIRED_DBG

TESTS = $(check_PROGRAMS)
