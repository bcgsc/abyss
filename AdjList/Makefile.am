bin_PROGRAMS = AdjList

AdjList_CPPFLAGS = -I$(top_srcdir) \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

if PAIRED_DBG
AdjList_CPPFLAGS += -DPAIRED_DBG
endif

AdjList_LDADD = \
	$(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a

if PAIRED_DBG
AdjList_LDADD += \
	$(top_builddir)/PairedDBG/libpaireddbg.a
endif

if HAVE_SQLITE
AdjList_LDADD += \
	$(top_builddir)/DataBase/libdb.a \
	$(SQLITE_LIBS)
endif

AdjList_SOURCES = \
	AdjList.cpp
