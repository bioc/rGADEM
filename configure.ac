AC_PREREQ(2.59)
AC_INIT([rGADEM], [1.1.1], [http://www.rglab.org])

AC_CONFIG_SRCDIR([DESCRIPTION])

# determine RHOME
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi

## pick up compiler as will be used by R CMD INSTALL/SHLIB
CC=`"${R_HOME}/bin/R" CMD config CC`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
AC_PROG_CC
AC_PROG_CPP


dnl # Checks for libraries.
dnl # FIXME: Replace `main' with a function in `-lm':
dnl AC_CHECK_LIB([m], [pow])

# Checks for header files.
AC_SUBST(OPENMP_CFLAGS)
AC_CHECK_HEADERS(dispatch/dispatch.h,
  :,
  [
    AC_CHECK_HEADERS(omp.h,
      [
        LIBS+=-lgomp
        OPENMP_CFLAGS=-fopenmp
      ],
      [
        AC_MSG_NOTICE([no parallelism])
      ])
  ]
)

dnl AC_HEADER_STDC
dnl AC_CHECK_HEADERS([limits.h malloc.h stdlib.h string.h])


dnl # Checks for typedefs, structures, and compiler characteristics.
dnl AC_C_CONST
dnl AC_TYPE_SIZE_T

dnl AC_CHECK_SIZEOF(long)

dnl AM_CONDITIONAL(LARGERINT, [test "$ac_cv_sizeof_long" = "8"])

dnl # Checks for library functions.
dnl AC_CHECK_FUNCS([pow fabs qsort floor rint sqrt])
dnl AC_FUNC_MALLOC
dnl AC_FUNC_REALLOC

AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_FILES([src/Makevars])

AC_OUTPUT
