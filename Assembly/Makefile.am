noinst_LIBRARIES = libassembly.a

libassembly_a_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

libassembly_a_SOURCES = \
	AssemblyAlgorithms.cc AssemblyAlgorithms.h \
	BranchGroup.h \
	BranchRecord.h \
	BranchRecordBase.h \
	DBG.h \
	DotWriter.h \
	Options.cc Options.h \
	SequenceCollection.h \
	VertexData.h \
	AdjacencyAlgorithm.h \
	AssembleAlgorithm.h \
	BubbleAlgorithm.h \
	CoverageAlgorithm.h \
	ErodeAlgorithm.h \
	LoadAlgorithm.h \
	SeqExt.h \
	SplitAlgorithm.h \
	TrimAlgorithm.h

if HAVE_SQLITE
libassembly_a_LIBADD = $(top_builddir)/Common/libcommon.a
endif
