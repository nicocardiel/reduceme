!------------------------------------------------------------------------------
! Version 20-November-2026                             File: avoid_warnings.f90
!------------------------------------------------------------------------------
! Copyright N. Cardiel & J. Gorgas, Departamento de Astrofisica
! Universidad Complutense de Madrid, 28040-Madrid, Spain
! E-mail: ncl@astrax.fis.ucm.es or fjg@astrax.fis.ucm.es
!------------------------------------------------------------------------------
! This routine is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the Free
! Software Foundation; either version 2 of the License, or (at your option) any
! later version. See the file gnu-public-license.txt for details.
!------------------------------------------------------------------------------
!Comment
!
! SUBROUTINE AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
!
! Input: STWV,DISP,NSCAN,NCHAN
! Output: (nothing)
!
! Dummy function to avoid compilation warnings "unused variable".
!
!Comment
!------------------------------------------------------------------------------
        SUBROUTINE AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
        IMPLICIT NONE
        REAL STWV
        REAL DISP
        INTEGER NSCAN
        INTEGER NCHAN
!------------------------------------------------------------------------------
        STWV=STWV
        DISP=DISP
        NSCAN=NSCAN
        NCHAN=NCHAN
        END
