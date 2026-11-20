# Process this file with autoconf to produce a configure script.

AC_INIT([reduceme],[5.1],[cardiel@ucm.es])
AM_INIT_AUTOMAKE([foreign])

AC_PROG_CC
AC_PROG_FC
AC_PROG_RANLIB

# Define directories to search for libraries
AC_MSG_CHECKING([uname -s for detecting host operating system])
OSNAME=`uname -s`
case "$OSNAME" in
  Darwin*)
    AC_MSG_RESULT(Darwin)
    # MacPorts library path (only if present, to avoid linker warnings)
    if test -d /opt/local/lib; then
      LDFLAGS="$LDFLAGS -L/opt/local/lib"
    fi
  ;;
  Linux*)
    AC_MSG_RESULT(Linux)
    # PGPLOT library path (only if present)
    if test -d /usr/local/pgplot; then
      LDFLAGS="$LDFLAGS -L/usr/local/pgplot"
    fi
  ;;
esac
AC_MSG_RESULT(LDFLAGS set to... $LDFLAGS)

# Checks for libraries
#PKG_CHECK_MODULES([CFITSIO], [cfitsio])
AC_CHECK_LIB([cfitsio], main,,AC_MSG_ERROR(cfitsio library not found))

# Optional GNU readline (file name completion in READC_FILE).
# --with-readline=DIR looks for DIR/include and DIR/lib.
AC_ARG_WITH([readline],
  [AS_HELP_STRING([--with-readline@<:@=DIR@:>@],
    [use readline for file name completion @<:@default=check@:>@])],
  [], [with_readline=check])
have_readline=no
AS_IF([test "x$with_readline" != xno], [
  rl_save_CPPFLAGS=$CPPFLAGS
  rl_save_LDFLAGS=$LDFLAGS
  rl_save_LIBS=$LIBS
  case "$with_readline" in
    yes|check)
      # Homebrew installs GNU readline as keg-only (not linked into the
      # default search paths, where the libedit headers of the SDK are found);
      # MacPorts installs it under /opt/local
      case "$OSNAME" in
        Darwin*)
          for rl_dir in /opt/homebrew/opt/readline /usr/local/opt/readline \
                        /opt/local; do
            if test -f "$rl_dir/include/readline/readline.h"; then
              CPPFLAGS="$CPPFLAGS -I$rl_dir/include"
              LDFLAGS="$LDFLAGS -L$rl_dir/lib"
              break
            fi
          done
        ;;
      esac
    ;;
    *)
      CPPFLAGS="$CPPFLAGS -I$with_readline/include"
      LDFLAGS="$LDFLAGS -L$with_readline/lib"
    ;;
  esac
  AC_LANG_PUSH([C])
  AC_CHECK_HEADER([readline/readline.h], [rl_header=yes], [rl_header=no],
    [#include <stdio.h>])
  AS_IF([test "x$rl_header" = xyes], [
    # some readline builds are not linked against a termcap library
    for rl_extra in "" -ltinfo -lncurses -ltermcap; do
      AC_MSG_CHECKING([for readline with "$rl_extra"])
      LIBS="-lreadline $rl_extra $rl_save_LIBS"
      AC_LINK_IFELSE([AC_LANG_CALL([], [readline])],
        [have_readline=yes])
      AC_MSG_RESULT([$have_readline])
      test "x$have_readline" = xyes && break
    done
  ])
  AS_IF([test "x$have_readline" = xyes], [
    AC_DEFINE([HAVE_READLINE], [1],
      [Define to 1 if readline is available.])
    AC_CHECK_HEADERS([readline/history.h], [], [], [#include <stdio.h>])
  ], [
    CPPFLAGS=$rl_save_CPPFLAGS
    LDFLAGS=$rl_save_LDFLAGS
    LIBS=$rl_save_LIBS
    AS_IF([test "x$with_readline" != xcheck],
      [AC_MSG_ERROR([readline requested but not found])])
  ])
  AC_LANG_POP([C])
])
AC_MSG_NOTICE([file name completion with readline: $have_readline])

AC_LANG(Fortran)
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
		 futils/redlib_m.f90
		 futils/Makefile
		 files/Makefile
		 src/Makefile
		 scripts/Makefile
		 tests/Makefile]) 

# Definitions
DEFAULT_NSMAX="2400"
NSMAX=${NSMAX="$DEFAULT_NSMAX"}
AC_ARG_VAR(NSMAX,[Maximum image dimension in the spatial direction])

DEFAULT_NCMAX="2400"
NCMAX=${NCMAX="$DEFAULT_NCMAX"}
AC_ARG_VAR(NCMAX,[Maximum image dimension in the spectral direction])

AC_OUTPUT

