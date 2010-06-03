bin_PROGRAMS = ParseAligns

ParseAligns_CPPFLAGS = \
	-I$(top_srcdir)/Align \
	-I$(top_srcdir)/Common

ParseAligns_LDADD = \
	$(top_builddir)/Align/libalign.a \
	$(top_builddir)/Common/libcommon.a

ParseAligns_SOURCES = \
	ParseAligns.cpp
