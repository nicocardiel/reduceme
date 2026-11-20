SUBDIRS = files button futils fspec red src scripts tests
ACLOCAL_AMFLAGS = -I m4

# Check that the Fortran sources conform to the Fortran 2018 standard (no
# compiler extensions); requires gfortran. See tests/check_std.sh.
check-std:
	@$(SHELL) $(top_srcdir)/tests/check_std.sh '$(FC)' '$(top_srcdir)' \
	  '$(top_builddir)'

.PHONY: check-std
