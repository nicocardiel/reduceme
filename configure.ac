# Process this file with autoconf to produce a configure script.

AC_INIT([reduceme],[5.0],[cardiel@ucm.es])
AM_INIT_AUTOMAKE([foreign])

AC_PROG_CC
AC_PROG_F77
AC_PROG_RANLIB

# Define directories to search for libraries
AC_MSG_CHECKING([uname -s for detecting host operating system])
OSNAME=`uname -s`
case "$OSNAME" in
  Darwin*)
    AC_MSG_RESULT(Darwin)
    LDFLAGS="$LDFLAGS -L/opt/local/lib"
  ;;
  Linux*)
    AC_MSG_RESULT(Linux)
    LDFLAGS="$LDFLAGS -L/usr/local/pgplot"
  ;;
esac
AC_MSG_RESULT(LDFLAGS set to... $LDFLAGS)

# Checks for libraries
#PKG_CHECK_MODULES([CFITSIO], [cfitsio])
AC_CHECK_LIB([cfitsio], main,,AC_MSG_ERROR(cfitsio library not found))
AC_LANG(Fortran 77)
AC_CHECK_LIB(pgplot, pgopen, [], [AC_MSG_ERROR([no pgplot found])])

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_LINKS([red/futils.inc:futils/futils.inc
		fspec/futils.inc:futils/futils.inc
		src/futils.inc:futils/futils.inc
		src/redlib.inc:red/redlib.inc
		fspec/redlib.inc:red/redlib.inc
		src/iofile.inc:red/iofile.inc])
AC_CONFIG_FILES([Makefile 
		 button/Makefile
		 fspec/Makefile
		 red/Makefile
		 red/redlib.inc
		 futils/Makefile
		 files/Makefile
		 src/Makefile
		 scripts/Makefile]) 

# Definitions
DEFAULT_NSMAX="2400"
NSMAX=${NSMAX="$DEFAULT_NSMAX"}
AC_ARG_VAR(NSMAX,[Maximum image dimension in the spatial direction])

DEFAULT_NCMAX="2400"
NCMAX=${NCMAX="$DEFAULT_NCMAX"}
AC_ARG_VAR(NCMAX,[Maximum image dimension in the spectral direction])

AC_OUTPUT

