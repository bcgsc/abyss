# programs that will only be compiled with 'make check'
check_PROGRAMS = BloomFilterTests \
	ParallelFilter

BloomFilterTests_SOURCES = BloomFilterTests.cpp
BloomFilterTests_CPPFLAGS = -I$(top_srcdir)
BloomFilterTests_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS) -std=c++11

ParallelFilter_SOURCES = ParallelFilter.cpp
ParallelFilter_CPPFLAGS = -I$(top_srcdir)
ParallelFilter_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS) -std=c++11

# programs that will be run by 'make check' (after compiling)
TESTS = $(check_PROGRAMS)
