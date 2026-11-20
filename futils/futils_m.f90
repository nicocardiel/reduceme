!------------------------------------------------------------------------------
! Version 20-November-2026                                   File: futils_m.f90
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
! Module FUTILS_M: terminal input/output and character utilities.
! The procedures are kept in their own files, included below; each one
! declares IMPLICIT NONE.
!------------------------------------------------------------------------------
        MODULE FUTILS_M
        CONTAINS
        INCLUDE 'avoid_warkings.f90'
        INCLUDE 'chlower.f90'
        INCLUDE 'chupper.f90'
        INCLUDE 'endprogram.f90'
        INCLUDE 'findmm.f90'
        INCLUDE 'findmml.f90'
        INCLUDE 'findmmlmask.f90'
        INCLUDE 'lrunx.f90'
        INCLUDE 'read2i.f90'
        INCLUDE 'readc.f90'
        INCLUDE 'readc_file.f90'
        INCLUDE 'readf.f90'
        INCLUDE 'readi.f90'
        INCLUDE 'readilim.f90'
        INCLUDE 'retorna_ticks.f90'
        INCLUDE 'rmblank.f90'
        INCLUDE 'showperc.f90'
        INCLUDE 'systemfunction.f90'
        INCLUDE 'truebeg.f90'
        INCLUDE 'truelen.f90'
        END MODULE FUTILS_M
