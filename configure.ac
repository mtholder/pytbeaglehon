#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_INIT(PYTBEAGLEHON, 0.0.1, mtholder@gmail.com)
AC_CONFIG_HEADERS([config.h])

AC_PREREQ(2.57)


################################################################################
#	According to (http://www.mail-archive.com/autoconf@gnu.org/msg14232.html)
#		this macro should be after AC_INIT but before AM_INIT_AUTOMAKE
################################################################################
AC_CONFIG_MACRO_DIR(config)

AM_INIT_AUTOMAKE([1.6 foreign dist-zip tar-ustar filename-length-max=299])

AC_CONFIG_SRCDIR([setup.py])

# Checks for programs.
AC_PROG_LN_S
AC_PROG_CPP
INP_CFLAGS="$CFLAGS"
AC_PROG_CC
CFLAGS="$INP_CFLAGS"

# required because we are building a library 
LT_INIT

AC_PROG_INSTALL
AC_PROG_LIBTOOL

NCL_INC_DIR="/usr/local/include"
NCL_LIB_DIR="/usr/local/lib"
NCL_BIN_DIR="/usr/local/bin"

AC_ARG_WITH(
	[ncl],
	AC_HELP_STRING(
		[--with-ncl=DIR],
		[Specify the root directory for the ncl version 2.1 or greater library (parent of the include/ncl and lib directories).  Omit this flag and configure will look for a system level NCL installation.  You might use the build_garli.sh script to automate the build process.]
		),
	[
	if ! test "$withval" = "yes" -o "$withval" = "no" ; then
		NCL_INC_DIR="$withval/include"
		NCL_LIB_DIR="$withval/lib"
		NCL_BIN_DIR="$withval/bin"
	fi
	])


BEAGLE_PREFIX="/usr/local"

AC_ARG_WITH(
	[beaglelib],
	AC_HELP_STRING(
		[--with-beaglelib=DIR],
		[Specify the prefix directory for the beaglelib (the parent of the include and lib directories).]
		),
	[
	if ! test "$withval" = "yes" -o "$withval" = "no" ; then
		BEAGLE_PREFIX="$withval"
	fi
	])

if ! test -d "$BEAGLE_PREFIX/include/libhmsbeagle-1" ; then
	AC_MSG_ERROR([beaglelib version 1.0 or higher is a prerequisite for building pytbeaglehon. Expecting a directory called $BEAGLE_PREFIX/include/libhmsbeagle-1, verify that the argument to --with-beaglelib is correct.])
else
	if ! test -f "$NCL_INC_DIR/ncl/nxsdefs.h" ; then
    	AC_MSG_ERROR([NCL version 2.1 or higher is a prerequisite for building pytbeaglehon. Expecting a file at $NCL_INC_DIR/ncl/nxsdefs.h, verify that the argument to --with-ncl is correct.])
	fi
fi

LIBS="$LIBS -lncl -lhmsbeagle"
LDFLAGS="$LDFLAGS -L$BEAGLE_PREFIX/lib"
CPPFLAGS="-I$BEAGLE_PREFIX/include/libhmsbeagle-1"
AC_SUBST(BEAGLE_PREFIX)
AC_SUBST(CFLAGS)
AC_SUBST(CPPFLAGS)

AC_HEADER_STDC
AC_CHECK_HEADERS([float.h stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE

# Checks for library functions.
AC_FUNC_MKTIME
AC_FUNC_STRTOD
AC_CHECK_FUNCS([memmove strchr strtol])

AC_PROG_MAKE_SET

AC_CONFIG_FILES([				\
	Makefile 					\
	pytbeaglehon/Makefile		\
	pytbeaglehon/ccore/Makefile	\
	c_example/Makefile	\
	])
AC_OUTPUT
