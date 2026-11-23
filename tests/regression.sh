#!/bin/sh
# Run the REDUCEME regression tests (called by "make check"): first as they
# are, then with FITS files instead of REDUCEME files (--fits).
# top_builddir and top_srcdir are exported by AM_TESTS_ENVIRONMENT.
if ! command -v python3 >/dev/null 2>&1; then
  echo "python3 not found: regression tests skipped"
  exit 77
fi
status=0
python3 "${top_srcdir:-..}/tests/run_tests.py" \
  --builddir "${top_builddir:-..}" --srcdir "${top_srcdir:-..}" || status=1
python3 "${top_srcdir:-..}/tests/run_tests.py" --fits \
  --builddir "${top_builddir:-..}" --srcdir "${top_srcdir:-..}" || status=1
exit $status
