AC_PREREQ(2.59)
AC_INIT(abyss, 1.0.4, abyss@bcgsc.ca)
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

# Checks for header files.
AC_CHECK_HEADERS([stdint.h stdlib.h sys/param.h])
AC_HEADER_STDC

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([gethostname memset pow setlocale sqrt \
	strdup strerror])
AC_FUNC_MEMCMP
AC_FUNC_SETVBUF_REVERSED
AC_FUNC_VPRINTF

# Checks for library constants.
AC_CHECK_DECL(HOST_NAME_MAX, [],
	AC_DEFINE(HOST_NAME_MAX, [_POSIX_HOST_NAME_MAX],
			  [Define if the system does not provide HOST_NAME_MAX]),
	[#include <limits.h>])

# Options to configure.
AC_ARG_WITH(mpi, AS_HELP_STRING([--with-mpi=PATH],
	[specify prefix directory for the installed MPI parallel
	computing library]))
if test "$with_mpi" -a -d "$with_mpi"; then
	mpi_cppflags="-I$with_mpi/include"
	mpi_ldflags="-L$with_mpi/lib"
fi

AC_ARG_ENABLE(lammpi, AS_HELP_STRING([--enable-lammpi],
	[use LAM/MPI (default is to use Open MPI)]))

AC_ARG_ENABLE(maxk, AS_HELP_STRING([--enable-maxk],
	[set the maximum k-mer length (default is 64)]),
	[max_kmer="$enableval"], [max_kmer=64])
AC_DEFINE_UNQUOTED(MAX_KMER, [$max_kmer], [maximum k-mer length]) 

# Set compiler flags.
AC_SUBST(CPPFLAGS, "$CPPFLAGS $mpi_cppflags")
AC_SUBST(LDFLAGS, "$LDFLAGS $mpi_ldflags")

# Check for the MPI parallel computing library.
if test "$enable_lammpi"; then
	AC_CHECK_LIB([pthread], [pthread_create])
	AC_CHECK_LIB([dl], [dlopen])
	AC_CHECK_LIB([lam], [lam_mutex_lock])
	AC_CHECK_LIB([mpi], [MPI_Init])
	AC_LANG_PUSH([C++])
	AC_CHECK_LIB([lammpi++], [main]) 
	AC_LANG_POP([C++])
else
	AC_CHECK_LIB([mpi], [MPI_Init])
	AC_DEFINE(OMPI_SKIP_MPICXX, 1,
			  [Define to disable OpenMPI C++ bindings])
fi
AM_CONDITIONAL([HAVE_LIBMPI], [test $ac_cv_lib_mpi_MPI_Init = yes])
AC_SUBST(MPI_LIBS, "$LIBS")
unset LIBS

# Check for the hash table implementation.
AC_LANG([C++])
AC_CHECK_HEADERS([ \
	google/sparse_hash_set \
	unordered_map tr1/unordered_map ext/hash_map \
	unordered_set tr1/unordered_set ext/hash_set \
])

# Set compiler flags.
AC_SUBST(AM_CXXFLAGS, '-Wall -Wextra -Werror')

AC_CONFIG_FILES([
	Makefile
	ABYSS/Makefile
	Align/Makefile
	Assembly/Makefile
	Common/Makefile
	DataLayer/Makefile
	ExpPaired/Makefile
	Graph/Makefile
	Parallel/Makefile
	bin/Makefile
	doc/Makefile
	fasta2psq/Makefile
	AdjList/Makefile
	DistanceEst/Makefile
	Overlap/Makefile
	SimpleGraph/Makefile
	MergePaths/Makefile
	KAligner/Makefile
	ParseAligns/Makefile
])
AC_OUTPUT
