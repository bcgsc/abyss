# programs that will only be compiled with 'make check'
check_PROGRAMS = \
    BloomFilterTests \
    CountingBloomFilterTests

BloomFilterTests_SOURCES = BloomFilterTests.cpp
BloomFilterTests_CPPFLAGS = -I$(top_srcdir)
BloomFilterTests_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS) -std=c++11

CountingBloomFilterTests_SOURCES = CountingBloomFilterTests.cpp
CountingBloomFilterTests_CPPFLAGS = -I$(top_srcdir)
CountingBloomFilterTests_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS) -std=c++11

# programs that will be run by 'make check' (after compiling)
TESTS = $(check_PROGRAMS)
