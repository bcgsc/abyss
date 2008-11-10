noinst_PROGRAMS = fasta2psq

fasta2psq_CPPFLAGS = -I$(top_srcdir)/Common -I$(top_srcdir)/DataLayer

fasta2psq_LDADD = \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a

fasta2psq_SOURCES = \
	fasta2psq.cpp fasta2psq.h \
	fasta2psqOptions.cpp fasta2psqOptions.h
