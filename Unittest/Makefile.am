TESTS=
UNIT_TESTS=
check_PROGRAMS=

UNIT_TESTS += common_stringutil
check_PROGRAMS += common_stringutil
common_stringutil_SOURCES = $(top_srcdir)/Common/test/StringUtilTest.cpp
common_stringutil_CPPFLAGS = $(CPPFLAGS) 
common_stringutil_CXXFLAGS = $(CXXFLAGS)
common_stringutil_LDADD = $(LIBRARY_PATH)/gtest_main.a -lpthread

TESTS += $(UNIT_TESTS)

