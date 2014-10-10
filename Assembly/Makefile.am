noinst_LIBRARIES = libassembly.a

libassembly_a_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

libassembly_a_SOURCES = \
	AssemblyAlgorithms.cc AssemblyAlgorithms.h \
	BranchGroup.cc BranchGroup.h \
	BranchRecord.cc BranchRecord.h \
	DotWriter.cc DotWriter.h \
	ISequenceCollection.h \
	KmerData.h \
	Options.cc Options.h \
	SequenceCollection.cc SequenceCollection.h \
	AdjacencyAlgorithm.h \
	AssembleAlgorithm.h \
	BubbleAlgorithm.h \
	CoverageAlgorithm.h \
	ErodeAlgorithm.h \
	LoadAlgorithm.h \
	SplitAlgorithm.h \
	TrimAlgorithm.h

if HAVE_SQLITE
libassembly_a_LIBADD = $(top_builddir)/Common/libcommon.a
endif
