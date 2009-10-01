bin_PROGRAMS = PathOverlap

PathOverlap_CPPFLAGS = \
	-I$(top_srcdir)/Align \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer \
	-I$(top_srcdir)/ExpPaired

PathOverlap_LDADD = \
	$(top_builddir)/Align/libalign.a \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a \
	$(top_builddir)/ExpPaired/libexppaired.a

PathOverlap_SOURCES = \
	PathOverlap.cpp
