TESTS=
UNIT_TESTS=
check_PROGRAMS=

UNIT_TESTS += common_stringutil
check_PROGRAMS += common_stringutil
common_stringutil_SOURCES = $(top_srcdir)/Common/test/StringUtilTest.cpp
common_stringutil_CPPFLAGS = -I$(top_srcdir)
common_stringutil_LDADD = $(LIBRARY_PATH)/gtest_main.a -lpthread

UNIT_TESTS += common_histogram
check_PROGRAMS += common_histogram
common_histogram_SOURCES = $(top_srcdir)/Common/test/HistogramTest.cpp
common_histogram_CPPFLAGS = -I$(top_srcdir)
common_histogram_LDADD = $(LIBRARY_PATH)/gtest_main.a -lpthread

TESTS += $(UNIT_TESTS)
