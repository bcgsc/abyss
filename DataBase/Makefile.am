noinst_LIBRARIES = libdb.a
libdb_a_SOURCES = DB.cc DB.h StatQueries.h Options.h
libdb_a_CPPFLAGS = -I$(top_srcdir)
libdb_a_LIBADD = $(top_builddir)/Common/libcommon.a

bin_PROGRAMS = DB_config abyss-db-csv
DB_config_SOURCES = DB_config.cc
abyss_db_csv_SOURCES = DB.cc DB.h db-csv.cc
abyss_db_csv_CPPFLAGS = -I$(top_srcdir)
