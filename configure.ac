AC_INIT(odpc, 2.0.2)

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi

RBIN="${R_HOME}/bin/R"
CC=`"${RBIN}" CMD config CC`
CPP=`"${RBIN}" CMD config CC -E`
CXX=`"${RBIN}" CMD config CXX`
CXXCPP=`"${RBIN}" CMD config CXX -E`
CXXFLAGS=`"${RBIN}" CMD config CXXFLAGS`
CPPFLAGS=`"${RBIN}" CMD config CPPFLAGS`

AC_PROG_CC
AC_PROG_CPP
AC_PROG_CXX

################################################################################
##
## Check for necessary support by C compiler
##
################################################################################
AC_LANG(C)
AC_CHECK_HEADERS(stdint.h inttypes.h limits.h float.h)

##
## check for -restrict CFLAGS support
##
oldCflags="$CFLAGS"
CFLAGS="$CFLAGS -restrict"
AC_MSG_CHECKING([whether C supports -restrict])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([])],
    [AC_MSG_RESULT([yes])]
    [RESTRICT_CFLAGS=-restrict],
    [AC_MSG_RESULT([no])]
)
CFLAGS="$oldCflags"
AC_SUBST([RESTRICT_CFLAGS])

##
## Check for OpenMP support
##
AC_OPENMP

# since some systems have broken OMP libraries
# we also check that the actual package will work
ac_pkg_openmp_c=no
if test -n "${OPENMP_CFLAGS}"; then
  AC_MSG_CHECKING([if OpenMP works in C])
  AC_LANG_CONFTEST(
  [AC_LANG_PROGRAM([[#include <omp.h>]], [[ return omp_get_num_threads (); ]])])
  PKG_CFLAGS="${OPENMP_CFLAGS}" PKG_LIBS="${OPENMP_CFLAGS}" "$RBIN" CMD SHLIB conftest.c 1>&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD && "$RBIN" --vanilla -q -e "dyn.load(paste('conftest',.Platform\$dynlib.ext,sep=''))" 1>&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD && ac_pkg_openmp_c=yes
  AC_MSG_RESULT([${ac_pkg_openmp_c}])
fi

# if ${ac_pkg_openmp_c} = "yes" then we have OMP, otherwise it will be "no"
if test "${ac_pkg_openmp_c}" = no; then
  OPENMP_CFLAGS=''
else
  AC_DEFINE([HAVE_OPENMP_C], [1], [Define to 1 if you have a C compiler with OpenMP support])
fi

AC_SUBST(OPENMP_CFLAGS)

################################################################################
##
## Check for necessary support by C++ compiler
##
################################################################################
AC_LANG(C++)
AC_CHECK_HEADERS(stdint.h inttypes.h climits cfloat)

# Check for "restrict" support
restrict_support="no"
AC_DEFINE([RESTRICT], [], [Set to the "restrict" keyword the C++ compiler understands.])

AC_MSG_CHECKING([whether C++ supports restrict keyword])
AC_COMPILE_IFELSE([
    AC_LANG_PROGRAM([[]],
                [[double *restrict ptr;]])
],
[
    AC_MSG_RESULT([yes])
    AC_DEFINE([RESTRICT],[restrict])
    restrict_support=yes
],
[AC_MSG_RESULT([no])]
)

if test $restrict_support != "yes"; then
    AC_MSG_CHECKING([whether C++ supports __restrict__ keyword])
    AC_COMPILE_IFELSE([
        AC_LANG_PROGRAM([[]],
                    [[double *__restrict__ ptr;]])
    ],
    [
        AC_MSG_RESULT([yes])
        AC_DEFINE([RESTRICT],[__restrict__])
        restrict_support="yes"
    ],
    [AC_MSG_RESULT([no])]
    )
fi

if test $restrict_support != "yes"; then
    AC_MSG_CHECKING([whether C++ supports __restrict keyword])
    AC_COMPILE_IFELSE([
        AC_LANG_PROGRAM([[]],
                    [[double *__restrict ptr;]])
    ],
    [
        AC_MSG_RESULT([yes])
        AC_DEFINE([RESTRICT],[__restrict])
        restrict_support=yes
    ],
    [AC_MSG_RESULT([no])]
    )
fi

##
## check for -restrict CXXFLAGS support
##
oldCXXflags="$CXXFLAGS"
CXXFLAGS="$CXXFLAGS -restrict"
AC_MSG_CHECKING([whether C++ supports -restrict])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([])],
    [AC_MSG_RESULT([yes])]
    [RESTRICT_CXXFLAGS=-restrict],
    [AC_MSG_RESULT([no])]
)
CXXFLAGS="$oldCXXflags"
AC_SUBST([RESTRICT_CXXFLAGS])


##
## Check for OpenMP support
##
AC_OPENMP

# since some systems have broken OMP libraries
# we also check that the actual package will work
ac_pkg_openmp_cxx=no
if test -n "${OPENMP_CXXFLAGS}"; then
  AC_MSG_CHECKING([if OpenMP works in C++])
  AC_LANG_CONFTEST(
  [AC_LANG_PROGRAM([[#include <omp.h>]], [[ return omp_get_num_threads (); ]])])
  PKG_CFLAGS="${OPENMP_CXXFLAGS}" PKG_LIBS="${OPENMP_CXXFLAGS}" "$RBIN" CMD SHLIB conftest.c 1>&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD && "$RBIN" --vanilla -q -e "dyn.load(paste('conftest',.Platform\$dynlib.ext,sep=''))" 1>&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD && ac_pkg_openmp_cxx=yes
  AC_MSG_RESULT([${ac_pkg_openmp_cxx}])
fi

# if ${ac_pkg_openmp_cxx} = "yes" then we have OMP, otherwise it will be "no"
if test "${ac_pkg_openmp_cxx}" = no; then
  OPENMP_CXXFLAGS=''
else
    AC_DEFINE([HAVE_OPENMP_C], [1], [Define to 1 if you have a C++ compiler with OpenMP support])
fi

AC_SUBST(OPENMP_CXXFLAGS)


AC_CONFIG_HEADERS([src/autoconfig.h])
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
