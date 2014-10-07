noinst_LIBRARIES = libpaireddbg.a

libpaireddbg_a_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

libpaireddbg_a_SOURCES = \
	AssemblyAlgorithms.cc AssemblyAlgorithms.h \
	BranchGroup.cc BranchGroup.h \
	BranchRecord.cc BranchRecord.h \
	DotWriter.cc DotWriter.h \
	ISequenceCollection.h \
	KmerData.h \
	Options.cc Options.h \
	SequenceCollection.cc SequenceCollection.h

if HAVE_SQLITE
libpaireddbg_a_LIBADD = $(top_builddir)/Common/libcommon.a
endif
