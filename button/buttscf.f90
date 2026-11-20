!------------------------------------------------------------------------------
! Version 20-November-2026                                    File: buttscf.f90
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
! SUBROUTINE BUTTSCF(FONT)
!
! Input : FONT
!
! Set the character font type in buttons.
!
! INTEGER FONT -> the current font number (in range 1-4)
!
!Comment
!------------------------------------------------------------------------------
!
        SUBROUTINE BUTTSCF(FONT)
        IMPLICIT NONE
        INTEGER FONT
!------------------------------------------------------------------------------
        IF((FONT.LT.1).OR.(FONT.GT.4))THEN
          WRITE(*,101)'ERROR: invalid font type in subroutine BUTTSCF.'
          WRITE(*,101)'=> previous font type unchanged'
          RETURN
        END IF
        PGSCF_BUTT=FONT
101     FORMAT(A)
        END
