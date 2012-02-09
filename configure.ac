# Process this file with autoconf to produce a configure script.

AC_INIT([reduceme],[4.2.03])
AM_INIT_AUTOMAKE([foreign])

AC_PROG_CC
AC_PROG_F77
AC_PROG_RANLIB

PKG_CHECK_MODULES([CFITSIO], [cfitsio])
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
		 futils/Makefile
		 src/Makefile
		 scripts/Makefile]) 
AC_OUTPUT

