UNIT_TESTS = common_stringutil
check_PROGRAMS = common_stringutil
common_stringutil_SOURCES = $(top_srcdir)/Common/test/StringUtilTest.cpp
common_stringutil_CPPFLAGS = -I$(top_srcdir)
common_stringutil_LDADD = -lgtest_main -lpthread

UNIT_TESTS += common_histogram
check_PROGRAMS += common_histogram
common_histogram_SOURCES = $(top_srcdir)/Common/test/HistogramTest.cpp
common_histogram_CPPFLAGS = -I$(top_srcdir)
common_histogram_LDADD = -lgtest_main -lpthread

TESTS = $(UNIT_TESTS)
