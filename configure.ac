AC_PREREQ([2.50])
AC_REVISION([$Id$])
AC_INIT([eGFRD],[0.9a],[gfrd@amolf.nl])
AC_CONFIG_SRCDIR([autogen.sh])
AC_DISABLE_STATIC
AC_LANG_CPLUSPLUS
AC_CANONICAL_TARGET([])
AM_INIT_AUTOMAKE([subdir-objects])

dnl checks for programs
AM_PROG_LIBTOOL
AC_PROG_PDFLATEX()
AC_PROG_CXX
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_RANLIB
AM_PATH_PYTHON(2.6)

dnl checks for libraries
AX_PATH_GSL([1.11],,AC_MSG_ERROR([could not find required version of GSL.]))
AC_CHECK_LIB(m,exp,,AC_MSG_ERROR([could not find libm.]))
AC_CHECK_LIB(python${PYTHON_VERSION},main,,AC_MSG_ERROR([could not find libpython.]))

DEBUG=
AC_ARG_ENABLE([debug],
  AC_HELP_STRING([--enable-debug],
                 [inhibit any optimization options]),
  [ DEBUG=1 ]
)

if test -n "$DEBUG"; then
  CXXFLAGS="-DDEBUG=1"
else
  CXXFLAGS="-g -O2"
fi
AC_SUBST(DEBUG)

CXXFLAGS="$CXXFLAGS -std=gnu++98 -Wall"

AX_BOOST_BASE([1.37],,AC_MSG_ERROR([could not find required version of BOOST.]))

BOOST_PYTHON_LIBNAME=boost_python

AC_ARG_WITH(boost-python-libname,
  AC_HELP_STRING([--with-boost-python-libname=LIBNAME],
                 [specify the library name of Boost.Python [[boost_python]]]
                ),
  [BOOST_PYTHON_LIBNAME=$withval]
)
AC_CHECK_LIB($BOOST_PYTHON_LIBNAME,main,,AC_MSG_ERROR([could not find libboost_python.]),-lpython${PYTHON_VERSION})
AC_SUBST(BOOST_PYTHON_LIBNAME)

BOOST_REGEX_LIBNAME=boost_regex

AC_ARG_WITH(boost-regex-libname,
  AC_HELP_STRING([--with-boost-date-time-libname=LIBNAME],
                 [specify the library name of Boost.Python [[boost_regex]]]
                ),
  [BOOST_REGEX_LIBNAME=$withval]
)
AC_CHECK_LIB($BOOST_REGEX_LIBNAME,main,,AC_MSG_ERROR([could not find libboost_regex.]))
AC_SUBST(BOOST_REGEX_LIBNAME)

AM_CHECK_PYMOD(numpy,,,[AC_MSG_ERROR([could not find Python numpy module.])])
AM_CHECK_PYMOD(scipy,,,[AC_MSG_ERROR([could not find Python scipy module.])])


dnl checks for header files
AC_CHECK_HEADERS([limits.h unistd.h])
dnl AC_CHECK_HEADER_STDBOOL

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_CHECK_FUNCS([memmove pow sqrt sincos isfinite])
AH_TEMPLATE(HAVE_SINCOS)
AH_TEMPLATE(HAVE_INLINE)
ECELL_CHECK_NUMPY
ECELL_CHECK_NUMPY_ARRAY_DESCR
ECELL_CHECK_LOGGING_MODULE
AM_CHECK_PYTHON_HEADERS()

ac_save_CPPFLAGS="$CPPFLAGS"
CPPFLAGS="$CPPFLAGS $PYTHON_INCLUDES"
AC_CHECK_TYPES([PyBaseExceptionObject],[],[],[
#include <Python.h>
#include <pyerrors.h>
])
CPPFLAGS="$ac_save_CPPFLAGS"
AC_CHECK_FUNCS([PyInt_FromSize_t],[],[],[
#include <Python.h>
])

AC_CHECK_HEADERS([unordered_map boost/unordered_map.hpp boost/functional/hash.hpp])
AC_CHECK_HEADER([tr1/unordered_map], [
  AC_LANG_SAVE
  AC_LANG_CPLUSPLUS
  AC_MSG_CHECKING([for a stupid bug in tr1/unordered_map])
  AC_TRY_COMPILE([
#include <tr1/unordered_map>
  ], [
std::tr1::unordered_map<int, int> a, b(a);
  ], [
    AC_MSG_RESULT([ok])
    AC_DEFINE([HAVE_TR1_UNORDERED_MAP], 1, [Define to 1 if you have usable <tr1/unordered_map> header file])
  ], [
    AC_MSG_RESULT([buggy])
  ])
  AC_LANG_RESTORE
], [
  AC_MSG_RESULT([buggy])
])

AC_CHECK_TYPES([std::hash<int>],[],[],[
#include <functional>
])
AC_CHECK_HEADERS([tr1/functional], [
  AC_CHECK_TYPES([std::tr1::hash<int>],[],[],[
#include <tr1/functional>
  ])
])

dnl checks for compiler characteristics
AC_C_INLINE
if test "$ac_cv_c_inline" != no ; then
   AC_DEFINE(HAVE_INLINE,1)
   AC_SUBST(HAVE_INLINE)
fi

dnl extract LTDL_SHLIB_EXT
rm -f conftest
./libtool --config > conftest
. ./conftest
rm -f conftest

dnl compatibility for libtool 1.5.6
LTDL_SHLIB_EXT=""
if test -n "$shrext_cmds"; then
    LTDL_SHLIB_EXT=$shrext_cmds
    AC_SUBST(LTDL_SHLIB_EXT)
dnl compatibility for libtool 1.5.0
elif test -n "$shrext"; then
    LTDL_SHLIB_EXT=$shrext
    AC_SUBST(LTDL_SHLIB_EXT)
fi

dnl checks for library functions
AC_CHECK_DECLS([INFINITY, HUGE_VAL],,,[
#include <math.h>
])

dnl no checks for system services
AM_CONFIG_HEADER([config.h])
AC_CONFIG_FILES([Makefile binding/Makefile test/Makefile doc/Makefile samples/benchmark/Makefile])
AC_OUTPUT([])
