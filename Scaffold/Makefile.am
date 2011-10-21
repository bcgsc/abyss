bin_PROGRAMS = abyss-scaffold

abyss_scaffold_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/Graph

abyss_scaffold_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

abyss_scaffold_LDADD = \
	$(top_builddir)/Common/libcommon.a

abyss_scaffold_SOURCES = scaffold.cc
