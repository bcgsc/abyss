AC_PREREQ(2.59)
AC_INIT(abyss, 0.4, abyss@bcgsc.ca)
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
AC_CHECK_FUNCS([memset mkdir pow sqrt])
AC_FUNC_MEMCMP
AC_FUNC_VPRINTF

# Options to configure.
AC_ARG_WITH(mpi, AS_HELP_STRING([--with-mpi=PATH],
	[specify prefix directory for the installed MPI parallel
	computing library]))
if test "$with_mpi" -a -d "$with_mpi"; then
	mpi_cppflags="-I$with_mpi/include"
	mpi_ldflags="-L$with_mpi/lib"
fi

AC_ARG_ENABLE(openmpi, AS_HELP_STRING([--enable-openmpi],
	[use Open MPI (default is to use LAM/MPI)]))

AC_ARG_ENABLE(maxk, AS_HELP_STRING([--enable-maxk],
	[set the maximum k-mer length (default is 64)]),
	[max_kmer="$enableval"], [max_kmer=64])
AC_DEFINE_UNQUOTED(MAX_KMER, [$max_kmer], [maximum k-mer length]) 

# Set compiler flags.
AC_SUBST(CPPFLAGS, "$CPPFLAGS $mpi_cppflags")
AC_SUBST(LDFLAGS, "$LDFLAGS $mpi_ldflags")

# Check for the LAM/MPI parallel computing library.
AC_LANG([C++])
if test "$enable_openmpi"; then
	AC_CHECK_LIB([mpi], [MPI_Init])
	AC_CHECK_LIB([mpi_cxx], [main])
else
	AC_CHECK_LIB([pthread], [pthread_create])
	AC_CHECK_LIB([dl], [dlopen])
	AC_CHECK_LIB([lam], [lam_mutex_lock])
	AC_CHECK_LIB([mpi], [MPI_Init])
	AC_CHECK_LIB([lammpi++], [main]) 
fi
AM_CONDITIONAL([HAVE_LIBMPI], [test $ac_cv_lib_mpi_MPI_Init = yes])

# Check for the Google sparse hash set.
AC_CHECK_HEADERS([google/sparse_hash_set])

# Set compiler flags.
AC_SUBST(AM_CXXFLAGS, '-Wall -Wextra -Werror')

AC_CONFIG_FILES([
	Makefile
	ABYSS/Makefile
	Align/Makefile
	Common/Makefile
	DataLayer/Makefile
	ExpPaired/Makefile
	Graph/Makefile
	Parallel/Makefile
	bin/Makefile
	fasta2psq/Makefile
	AdjList/Makefile
	DistanceEst/Makefile
	SimpleGraph/Makefile
	MergePaths/Makefile
	KAligner/Makefile
	ParseAligns/Makefile
])
AC_OUTPUT
