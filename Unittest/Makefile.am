# -Wno-error is used here because there is no portable way
# to suppress warning: "argument unused during compilation: '-pthread'"
# for clang on OSX.
# See: http://stackoverflow.com/questions/17841140/os-x-clang-pthread
AM_CXXFLAGS += $(PTHREAD_CFLAGS) -Wno-error
AM_LDFLAGS = $(PTHREAD_LIBS)
AM_CPPFLAGS = \
	-I$(top_srcdir) \
	-I$(top_srcdir)/lib/gtest-1.7.0/include
LDADD = $(top_builddir)/lib/gtest-1.7.0/libgtest_main.a

check_PROGRAMS = common_stringutil
common_stringutil_SOURCES = Common/StringUtilTest.cpp

check_PROGRAMS += common_histogram
common_histogram_SOURCES = Common/HistogramTest.cpp

check_PROGRAMS += common_bitutil
common_bitutil_SOURCES = Common/BitUtilTest.cpp

check_PROGRAMS += common_kmer
common_kmer_SOURCES = Common/KmerTest.cpp
common_kmer_LDADD = $(top_builddir)/Common/libcommon.a $(LDADD)

check_PROGRAMS += common_sequence
common_sequence_SOURCES = Common/Sequence.cc
common_sequence_LDADD = $(top_builddir)/Common/libcommon.a $(LDADD)

check_PROGRAMS += common_KmerIterator
common_KmerIterator_SOURCES = Common/KmerIteratorTest.cpp
common_KmerIterator_LDADD = $(top_builddir)/Common/libcommon.a $(LDADD)

check_PROGRAMS += common_sam
common_sam_SOURCES = Common/SAM.cc
common_sam_LDADD = $(top_builddir)/Common/libcommon.a $(LDADD)

check_PROGRAMS += BloomFilter
BloomFilter_SOURCES = Konnector/BloomFilter.cc
BloomFilter_CPPFLAGS = $(AM_CPPFLAGS) -I$(top_srcdir)/Common
BloomFilter_LDADD = $(top_builddir)/Common/libcommon.a $(LDADD)
BloomFilter_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

check_PROGRAMS += Konnector_DBGBloom
Konnector_DBGBloom_SOURCES = Konnector/DBGBloomTest.cpp
Konnector_DBGBloom_CPPFLAGS = $(AM_CPPFLAGS) -I$(top_srcdir)/Common
Konnector_DBGBloom_LDADD = $(top_builddir)/Common/libcommon.a $(LDADD)
Konnector_DBGBloom_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

check_PROGRAMS += Konnector_DBGBloomAlgorithms
Konnector_DBGBloomAlgorithms_SOURCES = Konnector/DBGBloomAlgorithmsTest.cpp
Konnector_DBGBloomAlgorithms_CPPFLAGS = $(AM_CPPFLAGS) -I$(top_srcdir)/Common
Konnector_DBGBloomAlgorithms_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)
Konnector_DBGBloomAlgorithms_LDADD = $(top_builddir)/Common/libcommon.a $(LDADD)

check_PROGRAMS += graph_ConstrainedBFSVisitor
graph_ConstrainedBFSVisitor_SOURCES = Graph/ConstrainedBFSVisitorTest.cpp
graph_ConstrainedBFSVisitor_CPPFLAGS = $(AM_CPPFLAGS) -I$(top_srcdir)/Common
graph_ConstrainedBFSVisitor_LDADD = $(top_builddir)/Common/libcommon.a $(LDADD)

check_PROGRAMS += graph_BidirectionalBFS
graph_BidirectionalBFS_SOURCES = Graph/BidirectionalBFSTest.cpp
graph_BidirectionalBFS_CPPFLAGS = $(AM_CPPFLAGS) -I$(top_srcdir)/Common
graph_BidirectionalBFS_LDADD = $(top_builddir)/Common/libcommon.a $(LDADD)

check_PROGRAMS += graph_AllPathsSearch
graph_AllPathsSearch_SOURCES = Graph/AllPathsSearchTest.cpp
graph_AllPathsSearch_CPPFLAGS = $(AM_CPPFLAGS) -I$(top_srcdir)/Common
graph_AllPathsSearch_LDADD = $(top_builddir)/Common/libcommon.a $(LDADD)

check_PROGRAMS += graph_HashGraph
graph_HashGraph_SOURCES = Graph/HashGraphTest.cpp
graph_HashGraph_CPPFLAGS = $(AM_CPPFLAGS) -I$(top_srcdir)/Common
graph_HashGraph_LDADD = $(top_builddir)/Common/libcommon.a $(LDADD)

check_PROGRAMS += graph_ConstrainedBidiBFSVisitor
graph_ConstrainedBidiBFSVisitor_SOURCES = \
	Graph/ConstrainedBidiBFSVisitorTest.cpp
graph_ConstrainedBidiBFSVisitor_CPPFLAGS = $(AM_CPPFLAGS) -I$(top_srcdir)/Common
graph_ConstrainedBidiBFSVisitor_LDADD = $(top_builddir)/Common/libcommon.a $(LDADD)

check_PROGRAMS += graph_ExtendPath
graph_ExtendPath_SOURCES = Graph/ExtendPathTest.cpp
graph_ExtendPath_CPPFLAGS = $(AM_CPPFLAGS) -I$(top_srcdir)/Common
graph_ExtendPath_LDADD = $(top_builddir)/Common/libcommon.a $(LDADD)

check_PROGRAMS += Konnector_konnector
Konnector_konnector_SOURCES = \
	Konnector/konnectorTest.cpp
Konnector_konnector_CPPFLAGS = $(AM_CPPFLAGS) -I$(top_srcdir)/Common
Konnector_konnector_LDADD = \
	$(top_builddir)/Align/libalign.a \
	$(top_builddir)/Common/libcommon.a \
	$(LDADD)
Konnector_konnector_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

check_PROGRAMS += DBG_LoadAlgorithm
DBG_LoadAlgorithm_SOURCES = \
	DBG/LoadAlgorithmTest.cpp
DBG_LoadAlgorithm_CPPFLAGS = \
	$(AM_CPPFLAGS) \
	-I$(top_srcdir)/DataLayer \
	-I$(top_srcdir)/Common
DBG_LoadAlgorithm_LDADD = \
	$(top_builddir)/Assembly/libassembly.a \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a \
	$(LDADD)
DBG_LoadAlgorithm_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

if PAIRED_DBG

check_PROGRAMS += PairedDBG_LoadAlgorithm
PairedDBG_LoadAlgorithm_SOURCES = \
	PairedDBG/LoadAlgorithmTest.cpp
PairedDBG_LoadAlgorithm_CPPFLAGS = \
	$(AM_CPPFLAGS) \
	-I$(top_srcdir)/DataLayer \
	-I$(top_srcdir)/Common
PairedDBG_LoadAlgorithm_LDADD = \
	$(top_builddir)/PairedDBG/libpaireddbg.a \
	$(top_builddir)/Assembly/libassembly.a \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a \
	$(LDADD)
PairedDBG_LoadAlgorithm_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

check_PROGRAMS += PairedDBG_KmerPair
PairedDBG_KmerPair_SOURCES = \
	PairedDBG/KmerPairTest.cc
PairedDBG_KmerPair_CPPFLAGS = \
	$(AM_CPPFLAGS) \
	-I$(top_srcdir)/DataLayer \
	-I$(top_srcdir)/Common
PairedDBG_KmerPair_LDADD = \
	$(top_builddir)/PairedDBG/libpaireddbg.a \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a \
	$(LDADD)
PairedDBG_KmerPair_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

check_PROGRAMS += PairedDBG_Dinuc
PairedDBG_Dinuc_SOURCES = \
	PairedDBG/DinucTest.cc
PairedDBG_Dinuc_CPPFLAGS = \
	$(AM_CPPFLAGS) \
	-I$(top_srcdir)/DataLayer \
	-I$(top_srcdir)/Common
PairedDBG_Dinuc_LDADD = \
	$(top_builddir)/PairedDBG/libpaireddbg.a \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a \
	$(LDADD)
PairedDBG_Dinuc_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

check_PROGRAMS += PairedDBG_BranchRecord
PairedDBG_BranchRecord_SOURCES = \
	PairedDBG/BranchRecordTest.cpp
PairedDBG_BranchRecord_CPPFLAGS = \
	$(AM_CPPFLAGS) \
	-I$(top_srcdir)/DataLayer \
	-I$(top_srcdir)/Common
PairedDBG_BranchRecord_LDADD = \
	$(top_builddir)/PairedDBG/libpaireddbg.a \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a \
	$(LDADD)
PairedDBG_BranchRecord_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

endif # PAIRED_DBG

check_PROGRAMS += BloomDBG_RollingHash
BloomDBG_RollingHash_SOURCES = BloomDBG/RollingHashTest.cpp

check_PROGRAMS += BloomDBG_RollingHashIterator
BloomDBG_RollingHashIterator_SOURCES = BloomDBG/RollingHashIteratorTest.cpp

check_PROGRAMS += BloomDBG_HashAgnosticCascadingBloom
BloomDBG_HashAgnosticCascadingBloom_SOURCES = \
	BloomDBG/HashAgnosticCascadingBloomTest.cpp
BloomDBG_HashAgnosticCascadingBloom_CXXFLAGS = $(AM_CXXFLAGS) \
	$(OPENMP_CXXFLAGS)

check_PROGRAMS += BloomDBG_RollingBloomDBG
BloomDBG_RollingBloomDBG_SOURCES = BloomDBG/RollingBloomDBGTest.cpp
BloomDBG_RollingBloomDBG_CXXFLAGS = $(AM_CXXFLAGS) \
	$(OPENMP_CXXFLAGS)
BloomDBG_RollingBloomDBG_LDADD = \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a \
	$(LDADD)

TESTS = $(check_PROGRAMS)
