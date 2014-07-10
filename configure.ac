AC_PREREQ(2.62)
AC_INIT(ABySS, 1.5.2, abyss-users@bcgsc.ca, abyss,
		http://www.bcgsc.ca/platform/bioinfo/software/abyss)
AM_INIT_AUTOMAKE(1.9.6 foreign subdir-objects)
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
AC_CHECK_TOOL(GHC_PKG, ghc-pkg)
AM_CONDITIONAL([HAVE_GHC_MMAP], ["$GHC_PKG" list |grep -q mmap])

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

# Check for std::hash and std::tr1::hash.
AC_LANG_PUSH([C++])
AC_CHECK_TYPE([std::hash<int>],
	AC_DEFINE(HAVE_STD_HASH, [1],
		[Define if the system provides std::hash]),
	[], [
#ifdef __FUJITSU
#include <string>
#else
#include <functional>
#endif])
AC_CHECK_TYPE([std::tr1::hash<int>],
	AC_DEFINE(HAVE_STD_TR1_HASH, [1],
		[Define if the system provides std::tr1::hash]),
	[], [#include <tr1/functional>])
AC_LANG_POP([C++])

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
	boost_cppflags="-isystem$with_boost -isystem$with_boost/include"
fi

# MPI
AC_ARG_WITH(mpi, AS_HELP_STRING([--with-mpi=PATH],
	[specify prefix directory for the installed MPI parallel
	computing library]))
if test "$with_mpi" -a -d "$with_mpi"; then
	mpi_cppflags="-isystem$with_mpi/include"
	mpi_ldflags="-L$with_mpi/lib"
fi

AC_ARG_ENABLE(mpich, AS_HELP_STRING([--enable-mpich],
	[use MPICH (default is to use Open MPI)]))
AC_ARG_ENABLE(lammpi, AS_HELP_STRING([--enable-lammpi],
	[use LAM/MPI (default is to use Open MPI)]))

# GTest
AC_ARG_WITH(gtest, AS_HELP_STRING([--with-gtest=PATH],
	[specify extracted folder for gtest]))
gtest_src="/usr/src/gtest"
if test "$with_gtest" -a -d "$with_gtest"; then
	gtest_cppflags="-isystem$with_gtest/include"
	gtest_src="$with_gtest"
fi
AC_SUBST(GTEST, "$gtest_src")

# SQLite
AC_ARG_WITH(sqlite, AS_HELP_STRING([--with-sqlite=PATH],
	[specify prefix directory for the installed sqlite library]))
if test "$with_sqlite" -a -d "$with_sqlite"; then
	sqlite_cppflags="-I$with_sqlite/include"
	sqlite_ldflags="-L$with_sqlite/lib -lsqlite3"
fi

# SparseHash
AC_ARG_WITH(sparsehash, AS_HELP_STRING([--with-sparsehash=PATH],
	[specify prefix directory for the installed spasehash library]))
if test "$with_sparsehash" -a "$with_sparsehash" != "no" -a -d "$with_sparsehash" ; then
	sparsehash_cppflags="-isystem$with_sparsehash/include"
	sparsehash_ldflags="-L$with_sparsehash/lib"
fi

AC_ARG_ENABLE(fm, AS_HELP_STRING([--enable-fm],
	[specify the width of the FM-index in bits (default is 64-bit)]),
	[], [enable_fm=64])
AC_DEFINE_UNQUOTED(FMBITS, $enable_fm,
				   [Width of bits of the FM-index in bits])

AC_ARG_ENABLE(maxk, AS_HELP_STRING([--enable-maxk=N],
	[set the maximum k-mer length (default is 64)]),
	[], [enable_maxk=64])
AC_DEFINE_UNQUOTED(MAX_KMER, [$enable_maxk], [maximum k-mer length])

# Find the absolute path to the source.
my_abs_srcdir=$(cd $srcdir; pwd)

# Set compiler flags.
boost_ver=1.56.0
boost_ver_dir=boost_1_56_0
AC_SUBST(CPPFLAGS,
		 "-I$my_abs_srcdir $boost_cppflags $mpi_cppflags $gtest_cppflags $sqlite_cppflags $sparsehash_cppflags $CPPFLAGS -isystem$my_abs_srcdir/$boost_ver_dir")
AC_SUBST(LDFLAGS, "$mpi_ldflags $sqlite_ldflags $sparsehash_ldflags $LDFLAGS")

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
LIBS="$libs"

# Check for the math library.
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
	functional \
	tr1/functional \
	boost/functional/hash.hpp \
	boost/property_map/property_map.hpp \
	google/sparse_hash_map \
	unordered_map tr1/unordered_map \
	unordered_set tr1/unordered_set \
	boost/unordered_set.hpp \
])

# Check for Boost.
if test $ac_cv_header_boost_property_map_property_map_hpp != yes; then
	AC_MSG_ERROR([ABySS requires the Boost C++ libraries, which may
	be downloaded from here: http://www.boost.org/users/download/
	It is not necessary to compile Boost before installing it. The
	following commands will download and install Boost for ABySS:
	cd $my_abs_srcdir
	wget http://downloads.sourceforge.net/project/boost/boost/$boost_ver/$boost_ver_dir.tar.bz2
	tar jxf $boost_ver_dir.tar.bz2
	cd -])
fi

# Check for gtest
libs="$LIBS"
AC_CHECK_HEADERS([gtest/gtest.h])
AC_CHECK_LIB([pthread], [pthread_create])
AM_CONDITIONAL([HAVE_GTEST], [test $ac_cv_header_gtest_gtest_h = yes])
AC_SUBST(GTEST_LIBS, "$LIBS")
LIBS=$libs

# Check for SQLite
libs="$LIBS"
AC_CHECK_HEADERS([sqlite3.h])
AC_CHECK_LIB([sqlite3],[main])
AM_CONDITIONAL([HAVE_SQLITE], [test $ac_cv_header_sqlite3_h = yes -a $ac_cv_lib_sqlite3_main = yes])
if (test $ac_cv_header_sqlite3_h = yes -a $ac_cv_lib_sqlite3_main = yes); then
	AC_DEFINE(_SQL, 1, [Define to 1 if you have sqlite lib/header])
	AC_SUBST(LDFLAGS, "-lsqlite3 $LDFLAGS")
fi
AC_SUBST(SQLITE_LIBS, "$LIBS")
LIBS=$libs

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
	Konnector/Makefile
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
	GapFiller/Makefile
	Sealer/Makefile
	Unittest/Makefile
	LogKmerCount/Makefile
	Bloom/Makefile
	DataBase/Makefile
])

if test "$with_sparsehash" != "no" -a "$ac_cv_header_google_sparse_hash_map" != "yes"; then
	AC_MSG_ERROR([ABySS should be compiled with Google sparsehash to
	reduce memory usage. It may be downloaded here:
	http://code.google.com/p/google-sparsehash

	If you do not wish to use sparsehash, specify --without-sparsehash.])
fi

AC_OUTPUT
