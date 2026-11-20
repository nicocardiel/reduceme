noinst_LIBRARIES = libbutton.a
# the module file includes the files with the procedures
libbutton_a_SOURCES = button_m.f90

BUTTON_ROUTINES = buttmis.f90 button.f90 buttqbr.f90 buttqcf.f90 buttqch.f90 buttqex.f90 \
	buttqit.f90 buttqpr.f90 buttqxb.f90 buttqyb.f90 buttqytext.f90 \
	buttsbr.f90 buttscf.f90 buttsch.f90 buttsex.f90 buttsit.f90 \
	buttspr.f90 buttsxb.f90 buttsyb.f90 buttsytext.f90 ifbutton.f90 \
	rpgband.f90 rpgbegin.f90 rpgbegok.f90 rpgenv.f90 rpgeras.f90 \
	rpgerasb.f90 rpgerasw.f90
EXTRA_DIST = $(BUTTON_ROUTINES)
button_m.$(OBJEXT): $(BUTTON_ROUTINES)

AM_FCFLAGS = -I. -I../futils
CLEANFILES = *.mod
