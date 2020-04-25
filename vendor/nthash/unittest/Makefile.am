# program that will only be compiled with 'make check'
check_PROGRAMS = UNITTESTS

UNITTESTS_SOURCES = UnitTests.cpp
UNITTESTS_CPPFLAGS = -I$(top_srcdir)
UNITTESTS_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS) -std=c++11

# programs that will be run by 'make check' (after compiling)
TESTS = $(check_PROGRAMS)
