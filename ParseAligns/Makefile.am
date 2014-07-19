bin_PROGRAMS = abyss-fixmate abyss-fixmate-ssq ParseAligns

abyss_fixmate_CPPFLAGS= -I$(top_srcdir) \
	-I$(top_srcdir)/Common

if HAVE_SQLITE
abyss_fixmate_LDADD= \
	$(top_builddir)/DataBase/libdb.a \
	$(SQLITE_LIBS) \
	$(top_builddir)/Common/libcommon.a
else
abyss_fixmate_LDADD= \
	$(top_builddir)/Common/libcommon.a
endif

abyss_fixmate_SOURCES=abyss-fixmate.cc

abyss_fixmate_ssq_CPPFLAGS = $(abyss_fixmate_CPPFLAGS) \
	-D SAM_SEQ_QUAL=1

abyss_fixmate_ssq_LDADD = $(abyss_fixmate_LDADD)

abyss_fixmate_ssq_SOURCES = $(abyss_fixmate_SOURCES)

ParseAligns_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common

ParseAligns_LDADD = \
	$(top_builddir)/Common/libcommon.a

ParseAligns_SOURCES = \
	ParseAligns.cpp
