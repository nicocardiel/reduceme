#!/bin/sh
# Check that the Fortran sources conform to the Fortran 2018 standard, i.e.
# that no compiler extensions are used (called by "make check-std").
# Requires gfortran. Obsolescent (but standard) features are not reported.
#
# Usage: check_std.sh [FC [TOP_SRCDIR [TOP_BUILDDIR]]]
FC=${1:-gfortran}
SRC=${2:-.}
BUILD=${3:-.}
# -Werror=intrinsics-std: GNU intrinsic procedures (GETENV, SYSTEM, RAND...)
# would otherwise be silently treated as external procedures. (-w cannot be
# used to hide the obsolescence warnings: it would also hide these errors.)
FLAGS="-std=f2018 -Werror=intrinsics-std -fsyntax-only"
# module files written by the check go to a temporary directory; the modules
# of the package are read from the build tree (run after make)
MODDIR=`mktemp -d "${TMPDIR:-/tmp}/check_std.XXXXXX"` || exit 1
trap 'rm -rf "$MODDIR"' EXIT
nfiles=0
nbad=0
for d in futils red button fspec src; do
  for f in "$SRC"/$d/*.f90; do
    nfiles=`expr $nfiles + 1`
    if ! out=`$FC $FLAGS -J"$MODDIR" -I"$MODDIR" -I"$BUILD/futils" \
        -I"$BUILD/$d" -I"$SRC/$d" "$f" 2>&1`; then
      nbad=`expr $nbad + 1`
      # print only the diagnostic blocks that end with an error
      printf '%s\n' "$out" | awk '
        /^Error:/   { printf "%s%s\n", block, $0; block = ""; next }
        /^Warning:/ { block = ""; next }
        { block = block $0 "\n" }'
    fi
  done
done
echo "check-std: $nfiles files checked, $nbad with non-standard code"
test $nbad -eq 0
