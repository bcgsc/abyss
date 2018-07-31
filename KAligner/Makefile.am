bin_PROGRAMS = KAligner

KAligner_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

KAligner_LDADD = \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a \
	-lpthread

KAligner_SOURCES = KAligner.cpp Aligner.cpp Aligner.h Options.h \
	Pipe.h PipeMux.h Semaphore.h
