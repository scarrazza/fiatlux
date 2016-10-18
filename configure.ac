#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT(fiatlux, 0.0.1, stefano.carrazza@cern.ch)
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([autoconfig])
#m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

LT_INIT()

AC_CONFIG_HEADERS([src/fiatlux/config.h])
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([bin/Makefile])
AC_CONFIG_FILES([bin/fiatlux-config])
AC_CONFIG_FILES([src/Makefile])
AC_CONFIG_FILES([examples/Makefile])
AM_INIT_AUTOMAKE([subdir-objects])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for dependencies
AC_SEARCH_LHAPDF

# Check for libraries.
AX_CXX_COMPILE_STDCXX_11([noext])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h sys/time.h quadmath.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([gettimeofday memset mkdir pow sqrt strstr])
AC_OUTPUT
