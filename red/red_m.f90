!------------------------------------------------------------------------------
! Version 24-November-2026                                      File: red_m.f90
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
! Module RED_M: REDUCEME image input/output and common routines.
! The procedures are kept in their own files, included below; each one
! declares IMPLICIT NONE.
!
! FITS files are read and written through scratch files in REDUCEME format
! (INFILEX, OUTFILEX, CLOSEFILEX); the private variables below keep the
! state shared by those routines.
!------------------------------------------------------------------------------
        MODULE RED_M
        USE FUTILS_M
        USE FSPEC_M
!
! FITS output files pending: units opened by OUTFILEX as scratch files that
! CLOSEFILEX must convert into the FITS file with the given name
        INTEGER, PARAMETER, PRIVATE :: NMAX_FITSOUT=100
        INTEGER, PRIVATE :: NFITSOUT=0
        INTEGER, PRIVATE :: IUNIT_FITSOUT(NMAX_FITSOUT)
        CHARACTER*255, PRIVATE :: NAME_FITSOUT(NMAX_FITSOUT)
! header cards of the first FITS data file read (not an error frame), copied
! into the FITS output files; NCARDS_FITSHDR=-1 until such a file is read
        INTEGER, PRIVATE :: NCARDS_FITSHDR=-1
        CHARACTER*80, ALLOCATABLE, PRIVATE :: CARDS_FITSHDR(:)
!
        CONTAINS
        INCLUDE 'autoplot.f90'
        INCLUDE 'closefilex.f90'
        INCLUDE 'fftcorrzoom.f90'
        INCLUDE 'fftprep.f90'
        INCLUDE 'fits2scratch.f90'
        INCLUDE 'fitserr.f90'
        INCLUDE 'guessef.f90'
        INCLUDE 'infilex.f90'
        INCLUDE 'isfitsfile.f90'
        INCLUDE 'isfitsname.f90'
        INCLUDE 'lkeepcard.f90'
        INCLUDE 'my_pgend.f90'
        INCLUDE 'outfilex.f90'
        INCLUDE 'palette.f90'
        INCLUDE 'pgiden_red.f90'
        INCLUDE 'pidegter.f90'
        INCLUDE 'scratch2fits.f90'
        INCLUDE 'showhlp.f90'
        INCLUDE 'welcome.f90'
        INCLUDE 'zscale.f90'
        END MODULE RED_M
