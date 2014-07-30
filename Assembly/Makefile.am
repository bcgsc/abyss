noinst_LIBRARIES = libassembly.a

libassembly_a_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

libassembly_a_SOURCES = \
	AssemblyAlgorithms.cpp AssemblyAlgorithms.h \
	BranchGroup.cpp BranchGroup.h \
	BranchRecord.cpp BranchRecord.h \
	DotWriter.cpp DotWriter.h \
	ISequenceCollection.h \
	KmerData.h \
	Options.cpp Options.h \
	SequenceCollection.cpp SequenceCollection.h

if HAVE_SQLITE
libassembly_a_LIBADD = $(top_builddir)/Common/libcommon.a
endif
