noinst_LIBRARIES = libred.a
# the module file includes the files with the procedures
libred_a_SOURCES = red_m.f90 pgplot.inc redlib.inc

RED_ROUTINES = autoplot.f90 fftcorrzoom.f90 fftprep.f90 guessef.f90 infilex.f90 \
	my_pgend.f90 outfilex.f90 palette.f90 pgiden_red.f90 pidegter.f90 \
	showhlp.f90 welcome.f90 zscale.f90
EXTRA_DIST = $(RED_ROUTINES)
red_m.$(OBJEXT): $(RED_ROUTINES)

AM_FCFLAGS = -I. -I../futils -I../fspec
CLEANFILES = *.mod
