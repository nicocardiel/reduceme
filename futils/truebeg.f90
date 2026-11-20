!------------------------------------------------------------------------------
! Version 20-November-2026                                    File: truebeg.f90
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
! INTEGER FUNCTION TRUEBEG(CADENA)
!
! Input: CADENA
! Output: TRUEBEG (function)
!
! Return the position of the first non-blank character in CADENA (ignoring
! also control characters with ASCII value < 32)
!
! CHARACTER*(*) CADENA -> input character string
!
!Comment
!------------------------------------------------------------------------------
        INTEGER FUNCTION TRUEBEG(CADENA)
        IMPLICIT NONE
        CHARACTER*(*) CADENA
!
        INTEGER I,L
!------------------------------------------------------------------------------
        L=LEN(CADENA)
!
        DO I=1,L
          IF(ICHAR(CADENA(I:I)).GT.32)THEN
            TRUEBEG=I
            RETURN
          END IF
        END DO
        TRUEBEG=0
        END
