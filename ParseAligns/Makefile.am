bin_PROGRAMS = abyss-fixmate ParseAligns

abyss_fixmate_CPPFLAGS= -I$(top_srcdir) \
	-I$(top_srcdir)/Align \
	-I$(top_srcdir)/Common

abyss_fixmate_LDADD= \
	$(top_builddir)/Align/libalign.a \
	$(top_builddir)/Common/libcommon.a

abyss_fixmate_SOURCES=abyss-fixmate.cc

ParseAligns_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Align \
	-I$(top_srcdir)/Common

ParseAligns_LDADD = \
	$(top_builddir)/Align/libalign.a \
	$(top_builddir)/Common/libcommon.a

ParseAligns_SOURCES = \
	ParseAligns.cpp
