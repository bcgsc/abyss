bin_PROGRAMS = KAligner

KAligner_CPPFLAGS = -pthread -I$(top_srcdir) \
	-I$(top_srcdir)/Align \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

KAligner_LDADD = \
	$(top_builddir)/Align/libalign.a \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a

KAligner_LDFLAGS = -pthread

KAligner_SOURCES = KAligner.cpp Barrier.h Pipe.h
