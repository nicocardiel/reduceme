!------------------------------------------------------------------------------
! Version 23-November-2026                               File: index_data_m.f90
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
! Module INDEX_DATA_M: image and error frame measured by MIDEIND, shared
! with the programs index and snratio (they were in COMMON /BLKINDEX0/),
! and response curves FLUX(channel,curve) (COMMON /BLKINDEX1/). The
! programs allocate them with the actual image dimensions (FLUX with
! NCHAN channels and NFMAX curves). This module is kept apart from FSPEC_M
! so that the generic names A and ERR are only visible in the program
! units that USE it.
!------------------------------------------------------------------------------
        MODULE INDEX_DATA_M
        IMPLICIT NONE
        REAL, ALLOCATABLE :: A(:,:),ERR(:,:)
        REAL, ALLOCATABLE :: FLUX(:,:)
        END MODULE INDEX_DATA_M
