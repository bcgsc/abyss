bin_PROGRAMS = KAligner

KAligner_CPPFLAGS = -pthread \
	-I$(top_srcdir)/Align \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer \
	-I$(top_srcdir)/ExpPaired

KAligner_LDADD = \
	$(top_builddir)/Align/libalign.a \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a \
	$(top_builddir)/ExpPaired/libexppaired.a

KAligner_LDFLAGS = -pthread

KAligner_SOURCES = \
	KAligner.cpp
