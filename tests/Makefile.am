# Regression tests: "make check" runs REDUCEME programs on synthetic data and
# compares the results with the reference results (see run_tests.py).
TESTS = regression.sh
AM_TESTS_ENVIRONMENT = top_builddir='$(abs_top_builddir)'; \
	top_srcdir='$(abs_top_srcdir)'; export top_builddir top_srcdir;
EXTRA_DIST = regression.sh run_tests.py mkdata.py check_std.sh README \
	MANUAL_CHECKS cases references
