noinst_LIBRARIES = libfspec.a
# the module file includes the files with the procedures
libfspec_a_SOURCES = fspec_m.f90
nodist_libfspec_a_SOURCES = redlib.inc

FSPEC_ROUTINES = bicubspl.f90 bicubsplx.f90 binsearch.f90 broaden.f90 cauchyfit.f90 \
	cfft.f90 cfftd.f90 chequea_fileindex.f90 chrebin.f90 combpf.f90 \
	cubspl.f90 cubsplx.f90 downhill.f90 factorialpf.f90 fatmext.f90 \
	fft2power.f90 fftcorrel.f90 fftcosbell.f90 fftkfilter.f90 \
	fintgauss.f90 fmean0.f90 fmean1.f90 fmean2.f90 fmedian1.f90 \
	fpercent.f90 fpoly.f90 gauscfit.f90 gauscfit_movel.f90 gauss2afit.f90 \
	gauss2bfit.f90 gaussfit.f90 gaussfitamp.f90 integtab.f90 lagrange.f90 \
	linefit.f90 lininterp.f90 ludcmp.f90 lusolv.f90 mideind.f90 \
	ordena1f.f90 ordena1f1i.f90 ordena1f2i.f90 ordena1i.f90 ordena2f.f90 \
	ordena2i.f90 pacent.f90 polfit.f90 polfitsig.f90 pseudofit.f90 \
	ranred.f90 rebining.f90 rinterp.f90 rvrebin.f90 scompext.f90 \
	selbands.f90 selindex.f90 sellines.f90 sexag.f90 shindex.f90 \
	splfit.f90 subprece.f90 ulogreb.f90
EXTRA_DIST = $(FSPEC_ROUTINES)
fspec_m.$(OBJEXT): $(FSPEC_ROUTINES)

AM_FCFLAGS = -I. -I../futils
CLEANFILES = *.mod
