bin_PROGRAMS = Consensus

Consensus_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Align \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

Consensus_LDADD = \
	$(top_builddir)/Align/libalign.a \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a

Consensus_SOURCES = \
	Consensus.cpp
