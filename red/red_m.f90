!------------------------------------------------------------------------------
! Version 20-November-2026                                      File: red_m.f90
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
!------------------------------------------------------------------------------
        MODULE RED_M
        USE FUTILS_M
        USE FSPEC_M
        CONTAINS
        INCLUDE 'autoplot.f90'
        INCLUDE 'fftcorrzoom.f90'
        INCLUDE 'fftprep.f90'
        INCLUDE 'guessef.f90'
        INCLUDE 'infilex.f90'
        INCLUDE 'my_pgend.f90'
        INCLUDE 'outfilex.f90'
        INCLUDE 'palette.f90'
        INCLUDE 'pgiden_red.f90'
        INCLUDE 'pidegter.f90'
        INCLUDE 'showhlp.f90'
        INCLUDE 'welcome.f90'
        INCLUDE 'zscale.f90'
        END MODULE RED_M
