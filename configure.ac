AC_PREREQ(2.59)
AC_INIT(abyss, 0.4, jsimpson@bcgsc.ca)
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
	[specify prefix directory for the installed LAM/MPI parallel
	computing library]))
if test "$with_mpi" -a -d "$with_mpi"; then
	mpi_cppflags="-I$with_mpi/include"
	mpi_ldflags="-L$with_mpi/lib"
fi

# Set compiler flags.
AC_SUBST(CPPFLAGS, "$CPPFLAGS $mpi_cppflags")
AC_SUBST(LDFLAGS, "$LDFLAGS $mpi_ldflags")

# Check for the LAM/MPI parallel computing library.
AC_LANG([C++])
AC_CHECK_LIB([pthread], [pthread_create])
AC_CHECK_LIB([dl], [dlopen])
AC_CHECK_LIB([lam], [lam_mutex_lock])
AC_CHECK_LIB([mpi], [MPI_Init])
AC_CHECK_LIB([lammpi++], [main]) 
AC_CHECK_LIB([lammpio], [MPI_File_open]) 
#AC_CHECK_LIB([mpiP], [main])
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
	Scaffold/Makefile
	bin/Makefile
	fasta2psq/Makefile
])
AC_OUTPUT
