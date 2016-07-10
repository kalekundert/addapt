#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([sgrna_monte_carlo], [0.0], [kale.kundert@ucsf.edu])
AC_PREFIX_DEFAULT([~/.local])
AC_CONFIG_AUX_DIR([build/aux])
AC_CONFIG_MACRO_DIR([build/m4])
AC_CONFIG_FILES([Makefile])
AM_INIT_AUTOMAKE([foreign subdir-objects -Wall -Werror])
AM_SILENT_RULES([yes])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_INSTALL
LT_INIT

# Checks for C++11.
AC_LANG([C++])
AX_CXX_COMPILE_STDCXX_11([noext], [mandatory])

# Checks for library functions.
AC_OPENMP
AC_SUBST(OPENMP_CXXFLAGS)

AC_CHECK_LIB(
    [RNA], [vrna_pf],
    [], [AC_MSG_ERROR([missing the ViennaRNA library])],
    [$OPENMP_CXXFLAGS])

AX_CXX_CHECK_LIB(
    [docopt], [docopt::value],
    [], [AC_MSG_ERROR([missing the docopt.cpp library])])

# Checks for header files.
AC_CHECK_HEADERS(
    [ViennaRNA/part_func.h],
    [], [AC_MSG_ERROR([missing the ViennaRNA headers])])

AC_CHECK_HEADERS(
    [docopt/docopt.h],
    [], [AC_MSG_ERROR([missing the docopt.cpp headers])])

AC_CHECK_HEADERS(
    [boost/format.hpp boost/algorithm/string.hpp],
    [], [AC_MSG_ERROR([missing the boost headers])])

AC_CHECK_HEADERS(
    [algorithm cmath iostream iterator list memory regex string utility vector],
    [], [AC_MSG_ERROR([missing the C++11 headers])])

# Make the build scripts.
AC_OUTPUT

echo \
"--------------------------------------------------

  $PACKAGE_NAME $PACKAGE_VERSION

  Prefix.........: $prefix
  C++ Compiler...: $CXX $CXXFLAGS $CPPFLAGS

--------------------------------------------------"
