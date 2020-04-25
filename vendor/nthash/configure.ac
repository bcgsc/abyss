AC_PREREQ(2.62)
AC_INIT(NTHASH, 2.1.0, hmohamadi@bcgsc.ca, ntHash,
		https://github.com/bcgsc/ntHash)
# add 'foreign' to allow README.md instead of README
AM_INIT_AUTOMAKE([foreign subdir-objects])
AC_CONFIG_SRCDIR([nttest.cpp])
AC_CONFIG_HEADER([config.h])

if test -z $CXXFLAGS; then
    CXXFLAGS='-O3'
fi

if test -z $CCFLAGS; then
    CCFLAGS='-O3'
fi

# Checks for programs.
AC_PROG_AWK
AC_PROG_CC
AC_PROG_CPP
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_RANLIB

# Checks for header files.
AC_CHECK_HEADERS([dlfcn.h fcntl.h float.h limits.h \
	stddef.h stdint.h stdlib.h sys/param.h zlib.h])
AC_HEADER_STDBOOL
AC_HEADER_STDC

# Checks for library functions.
AC_CHECK_FUNCS([dup2 gethostname getopt_long getpagesize \
				memset strdup strerror strtoul])
AC_FUNC_FORK
AC_FUNC_MALLOC
AC_FUNC_MEMCMP
AC_FUNC_REALLOC
AC_FUNC_SETVBUF_REVERSED
AC_FUNC_VPRINTF

# Checks for typedefs, structures, and compiler characteristics.
AC_C_BIGENDIAN
AC_C_CONST
AC_C_INLINE
AC_CHECK_TYPES([ptrdiff_t])
AC_TYPE_MODE_T
AC_TYPE_PID_T
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T
AC_TYPE_INT64_T
AC_TYPE_UINT8_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T

# Check for the math library.
AC_CHECK_LIB([m], [sqrt])
AC_CHECK_FUNCS([pow sqrt])
AC_CHECK_FUNC(ceilf, [], AC_DEFINE(ceilf, [ceil],
			  [Define if the system does not provide ceilf]))

# Check for the dynamic linking library.
AC_CHECK_LIB([dl], [dlsym])

# zlib (optional)
AC_CHECK_LIB([z], [gzwrite])
AM_CONDITIONAL([HAVE_LIBZ],
	[test $ac_cv_header_zlib_h = yes -a $ac_cv_lib_z_gzwrite = yes])

# Check for OpenMP.
AC_LANG_PUSH([C++])
AC_OPENMP
AC_LANG_POP([C++])

# Find the absolute path to the source.
my_abs_srcdir=$(cd $srcdir; pwd)

# Set compiler flags
AC_SUBST(CPPFLAGS,
		 "-I$my_abs_srcdir $OPENMP_CXXFLAGS $CPPFLAGS")
AC_SUBST(LDFLAGS, "$OPENMP_CXXFLAGS $LDFLAGS")

AC_SUBST(AM_CXXFLAGS, '-Wall -Wextra -Werror')


AC_CONFIG_FILES([
	Makefile
	unittest/Makefile
	vendor/Makefile
])
AC_OUTPUT
