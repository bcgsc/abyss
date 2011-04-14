bin_PROGRAMS = DAssembler

DAssembler_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

DAssembler_LDADD = \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a

DAssembler_SOURCES = DAssembler.cpp \
	RotatedRead.cpp RotatedRead.h \
	Rotation.h
