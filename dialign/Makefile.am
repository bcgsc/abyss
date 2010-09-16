noinst_LIBRARIES = libdialign.a
noinst_PROGRAMS = dialign

libdialign_a_SOURCES = \
	alig.c \
	assemble.c \
	diag.c \
	io.c io.h \
	orf.c orf.h \
	parameters.c parameters.h \
	prob.c \
	struct.h \
	translate.c translate.h

dialign_SOURCES = museq.c

dialign_LDADD = libdialign.a
