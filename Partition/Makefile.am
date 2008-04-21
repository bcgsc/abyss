bin_PROGRAMS = Partition

Partition_CPPFLAGS = \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

Partition_LDFLAGS = \
	-L$(top_builddir)/Common \
	-L$(top_builddir)/DataLayer

Partition_LDADD = -lcommon -ldatalayer

Partition_SOURCES = \
	Partition.cpp Partition.h
