noinst_LIBRARIES = libdb.a
libdb_a_SOURCES = DB.cc DB.h Options.h
libdb_a_CPPFLAGS = -I$(top_srcdir)
libdb_a_LIBADD = $(top_builddir)/Common/libcommon.a

if HAVE_SQLITE3
bin_PROGRAMS = abyss-db-csv
endif

abyss_db_csv_SOURCES = DB.cc DB.h db-csv.cc
abyss_db_csv_CPPFLAGS = -I$(top_srcdir)
abyss_db_csv_LDADD = -lsqlite3

dist_bin_SCRIPTS = abyss-db-txt
