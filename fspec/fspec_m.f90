!------------------------------------------------------------------------------
! Version 20-November-2026                                    File: fspec_m.f90
!------------------------------------------------------------------------------
! Copyright N. Cardiel & J. Gorgas, Departamento de Astrofisica
! Universidad Complutense de Madrid, 28040-Madrid, Spain
! E-mail: ncl@astrax.fis.ucm.es or fjg@astrax.fis.ucm.es
!------------------------------------------------------------------------------
! This file is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the Free
! Software Foundation; either version 2 of the License, or (at your option) any
! later version. See the file gnu-public-license.txt for details.
!------------------------------------------------------------------------------
! Module FSPEC_M: spectral analysis routines.
! The procedures are kept in their own files, included below. There is
! no IMPLICIT NONE at module level because a few procedures (in
! linefit.f90 and polfit.f90) rely on implicit typing; the other
! procedures declare IMPLICIT NONE.
!------------------------------------------------------------------------------
        MODULE FSPEC_M
        USE FUTILS_M
        CONTAINS
        INCLUDE 'bicubspl.f90'
        INCLUDE 'bicubsplx.f90'
        INCLUDE 'binsearch.f90'
        INCLUDE 'broaden.f90'
        INCLUDE 'cauchyfit.f90'
        INCLUDE 'cfft.f90'
        INCLUDE 'cfftd.f90'
        INCLUDE 'chequea_fileindex.f90'
        INCLUDE 'chrebin.f90'
        INCLUDE 'combpf.f90'
        INCLUDE 'cubspl.f90'
        INCLUDE 'cubsplx.f90'
        INCLUDE 'downhill.f90'
        INCLUDE 'factorialpf.f90'
        INCLUDE 'fatmext.f90'
        INCLUDE 'fft2power.f90'
        INCLUDE 'fftcorrel.f90'
        INCLUDE 'fftcosbell.f90'
        INCLUDE 'fftkfilter.f90'
        INCLUDE 'fintgauss.f90'
        INCLUDE 'fmean0.f90'
        INCLUDE 'fmean1.f90'
        INCLUDE 'fmean2.f90'
        INCLUDE 'fmedian1.f90'
        INCLUDE 'fpercent.f90'
        INCLUDE 'fpoly.f90'
        INCLUDE 'gauscfit.f90'
        INCLUDE 'gauscfit_movel.f90'
        INCLUDE 'gauss2afit.f90'
        INCLUDE 'gauss2bfit.f90'
        INCLUDE 'gaussfit.f90'
        INCLUDE 'gaussfitamp.f90'
        INCLUDE 'integtab.f90'
        INCLUDE 'lagrange.f90'
        INCLUDE 'linefit.f90'
        INCLUDE 'lininterp.f90'
        INCLUDE 'ludcmp.f90'
        INCLUDE 'lusolv.f90'
        INCLUDE 'mideind.f90'
        INCLUDE 'ordena1f.f90'
        INCLUDE 'ordena1f1i.f90'
        INCLUDE 'ordena1f2i.f90'
        INCLUDE 'ordena1i.f90'
        INCLUDE 'ordena2f.f90'
        INCLUDE 'ordena2i.f90'
        INCLUDE 'pacent.f90'
        INCLUDE 'polfit.f90'
        INCLUDE 'polfitsig.f90'
        INCLUDE 'pseudofit.f90'
        INCLUDE 'ranred.f90'
        INCLUDE 'rebining.f90'
        INCLUDE 'rinterp.f90'
        INCLUDE 'rvrebin.f90'
        INCLUDE 'scompext.f90'
        INCLUDE 'selbands.f90'
        INCLUDE 'selindex.f90'
        INCLUDE 'sellines.f90'
        INCLUDE 'sexag.f90'
        INCLUDE 'shindex.f90'
        INCLUDE 'splfit.f90'
        INCLUDE 'subprece.f90'
        INCLUDE 'ulogreb.f90'
        END MODULE FSPEC_M
