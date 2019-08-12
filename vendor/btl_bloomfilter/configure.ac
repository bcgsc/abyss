#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(BLOOMFILTER, 1.2.1, cjustin@bcgsc.ca)
AC_CANONICAL_SYSTEM

AM_INIT_AUTOMAKE()

AC_CONFIG_HEADER([config.h])

AC_PROG_RANLIB


# Checks for programs.
AC_PROG_CXX

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Find the absolute paths.
my_abs_srcdir=$(cd $srcdir; pwd)

# Checks for header files.
AC_LANG([C++])	#sets the auto conf to look for c++ headers

# Check for OpenMP.
AC_OPENMP
if test -z $OPENMP_CXXFLAGS; then
	OPENMP_CXXFLAGS=-Wno-unknown-pragmas
fi

# Set compiler flags.
AC_SUBST(AM_CXXFLAGS, '-Wall -Wextra -Werror')

# Checks for typedefs, structures, and compiler characteristics.

AC_CONFIG_FILES([
		Makefile
		Tests/Unit/Makefile
		Tests/AdHoc/Makefile
])
AC_OUTPUT
