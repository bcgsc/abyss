bin_PROGRAMS = ParseAligns

ParseAligns_CPPFLAGS = \
	-I$(top_srcdir)/Align \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/ExpPaired

ParseAligns_LDADD = \
	$(top_builddir)/Align/libalign.a \
	$(top_builddir)/Common/libcommon.a \
	$(top_builddir)/ExpPaired/libexppaired.a

ParseAligns_SOURCES = \
	ParseAligns.cpp
