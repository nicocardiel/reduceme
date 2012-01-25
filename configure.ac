# Process this file with autoconf to produce a configure script.

AC_INIT([reduceme],[4.2.03])
AM_INIT_AUTOMAKE([foreign])
LT_INIT([disable-shared])

AC_PROG_F77
AC_PROG_LIBTOOL

PKG_CHECK_MODULES([CFITSIO], [cfitsio])

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_FILES([Makefile 
		 button/Makefile
		 futils/Makefile])
AC_OUTPUT

