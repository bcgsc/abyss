AC_PREREQ(2.59)
AC_INIT(ABySS, 1.3.5, abyss-users@bcgsc.ca, abyss,
		http://www.bcgsc.ca/platform/bioinfo/software/abyss)
AM_INIT_AUTOMAKE(foreign)
AC_CONFIG_SRCDIR([ABYSS/Abyss.cpp])
AC_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_AWK
AC_PROG_CC
AC_PROG_CPP
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_RANLIB
AC_CHECK_TOOL(GHC, ghc)
AM_CONDITIONAL([HAVE_GHC], [test "$GHC"])

# Checks for header files.
AC_CHECK_HEADERS([dlfcn.h fcntl.h float.h limits.h \
	stddef.h stdint.h stdlib.h sys/param.h])
AC_HEADER_STDBOOL
AC_HEADER_STDC

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

# Checks for library functions.
AC_CHECK_FUNCS([dup2 gethostname getopt_long getpagesize \
				memset strdup strerror strtoul])
AC_FUNC_FORK
AC_FUNC_MALLOC
AC_FUNC_MEMCMP
AC_FUNC_REALLOC
AC_FUNC_SETVBUF_REVERSED
AC_FUNC_VPRINTF

# Checks for library constants.
AC_CHECK_DECL(HOST_NAME_MAX, [],
	AC_DEFINE(HOST_NAME_MAX, [_POSIX_HOST_NAME_MAX],
			  [Define if the system does not provide HOST_NAME_MAX]),
	[#include <limits.h>])

# Options to configure.
# Boost
AC_ARG_WITH(boost, AS_HELP_STRING([--with-boost=PATH],
			[specify directory for the boost header files]))
if test "$with_boost" -a -d "$with_boost"; then
	boost_cppflags="-I$with_boost"
fi
# MPI
AC_ARG_WITH(mpi, AS_HELP_STRING([--with-mpi=PATH],
	[specify prefix directory for the installed MPI parallel
	computing library]))
if test "$with_mpi" -a -d "$with_mpi"; then
	mpi_cppflags="-I$with_mpi/include"
	mpi_ldflags="-L$with_mpi/lib"
fi

AC_ARG_ENABLE(mpich, AS_HELP_STRING([--enable-mpich],
	[use MPICH (default is to use Open MPI)]))
AC_ARG_ENABLE(lammpi, AS_HELP_STRING([--enable-lammpi],
	[use LAM/MPI (default is to use Open MPI)]))

AC_ARG_ENABLE(fm, AS_HELP_STRING([--enable-fm],
	[specify the width of the FM-index in bits (default is 64-bit)]),
	[], [enable_fm=64])
AC_DEFINE_UNQUOTED(FMBITS, $enable_fm,
				   [Width of bits of the FM-index in bits])

AC_ARG_ENABLE(maxk, AS_HELP_STRING([--enable-maxk=N],
	[set the maximum k-mer length (default is 64)]),
	[], [enable_maxk=64])
AC_DEFINE_UNQUOTED(MAX_KMER, [$enable_maxk], [maximum k-mer length])

AC_ARG_ENABLE(samseqqual,
	AS_HELP_STRING([--enable-samseqqual],
		[enable SAM sequence and quality fields]),
	AC_DEFINE_UNQUOTED(SAM_SEQ_QUAL, 1,
		[Define to use SAM sequence and quality fields]))

# Set compiler flags.
AC_SUBST(CPPFLAGS,
		 "-I$srcdir $boost_cppflags $mpi_cppflags $CPPFLAGS")
AC_SUBST(LDFLAGS, "$mpi_ldflags $LDFLAGS")

# Check for the MPI parallel computing library.
libs="$LIBS"
AC_DEFINE(MPICH_SKIP_MPICXX, 1,
		  [Define to disable MPICH C++ bindings])
AC_DEFINE(OMPI_SKIP_MPICXX, 1,
		  [Define to disable OpenMPI C++ bindings])
AC_CHECK_HEADERS([mpi.h])
if test "$enable_mpich"; then
	AC_CHECK_LIB([pthread], [pthread_create])
	AC_CHECK_LIB([mpl], [MPL_env2int])
	AC_CHECK_LIB([mpich], [MPI_Init])
	ac_cv_lib_mpi_MPI_Init=$ac_cv_lib_mpich_MPI_Init
elif test "$enable_lammpi"; then
	AC_CHECK_LIB([pthread], [pthread_create])
	AC_CHECK_LIB([dl], [dlopen])
	AC_CHECK_LIB([lam], [lam_mutex_lock])
	AC_CHECK_LIB([mpi], [MPI_Init])
	AC_LANG_PUSH([C++])
	AC_CHECK_LIB([lammpi++], [main]) 
	AC_LANG_POP([C++])
else
	AC_CHECK_LIB([mpi], [MPI_Init])
fi
AM_CONDITIONAL([HAVE_LIBMPI],
	[test $ac_cv_header_mpi_h = yes -a $ac_cv_lib_mpi_MPI_Init = yes])
AC_SUBST(MPI_LIBS, "$LIBS")

# Check for the math library.
LIBS="$libs"
AC_CHECK_LIB([m], [sqrt])
AC_CHECK_FUNCS([pow sqrt])
AC_CHECK_FUNC(ceilf, [], AC_DEFINE(ceilf, [ceil],
			  [Define if the system does not provide ceilf]))

# Check for the dynamic linking library.
AC_CHECK_LIB([dl], [dlsym])

# Check for popcnt instruction.
AC_COMPILE_IFELSE(
	[AC_LANG_PROGRAM([[#include <stdint.h>],
					  [uint64_t x = 0;]],
					 [[__asm__("popcnt %1,%0" : "=r" (x) : "r" (x));]])],
	[AC_DEFINE([HAVE_POPCNT], 1, [Define to 1 if you have popcnt.])],
	[AC_DEFINE([HAVE_POPCNT], 0, [Define to 0 if you do not have popcnt.])])

# Check for the hash table implementation.
AC_LANG([C++])
AC_CHECK_HEADERS([ \
	boost/property_map/property_map.hpp \
	google/sparse_hash_map \
	unordered_map tr1/unordered_map \
])

# Check for Boost.
if test $ac_cv_header_boost_property_map_property_map_hpp != yes; then
	AC_MSG_ERROR([ABySS requires the Boost C++ libraries, which may
	be downloaded from here: http://www.boost.org/users/download/
	It is not necessary to compile Boost before installing it. The
	following commands will download and install Boost for ABySS:
	wget http://downloads.sourceforge.net/project/boost/boost/1.50.0/boost_1_50_0.tar.bz2
	tar jxf boost_1_50_0.tar.bz2
	ln -s boost_1_50_0/boost boost
	])
fi

# Check for OpenMP.
AC_OPENMP
if test -z $OPENMP_CXXFLAGS; then
	OPENMP_CXXFLAGS=-Wno-unknown-pragmas
fi

# Set compiler flags.
AC_SUBST(AM_CXXFLAGS, '-Wall -Wextra -Werror')

AC_CONFIG_FILES([
	Makefile
	ABYSS/Makefile
	Align/Makefile
	Assembly/Makefile
	Common/Makefile
	DataLayer/Makefile
	FMIndex/Makefile
	Graph/Makefile
	Parallel/Makefile
	bin/Makefile
	doc/Makefile
	dialign/Makefile
	kmerprint/Makefile
	AdjList/Makefile
	DAssembler/Makefile
	DistanceEst/Makefile
	Layout/Makefile
	Map/Makefile
	Misc/Makefile
	Overlap/Makefile
	PopBubbles/Makefile
	Scaffold/Makefile
	SimpleGraph/Makefile
	MergePaths/Makefile
	KAligner/Makefile
	ParseAligns/Makefile
	PathOverlap/Makefile
	Consensus/Makefile
	FilterGraph/Makefile
])
AC_OUTPUT

if test $ac_cv_header_google_sparse_hash_map != yes; then
	AC_MSG_WARN([ABySS should be compiled with Google sparsehash to
	reduce memory usage. It may be downloaded here:
	http://code.google.com/p/google-sparsehash])
fi
