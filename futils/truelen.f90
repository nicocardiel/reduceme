!------------------------------------------------------------------------------
! Version 25-November-1996                                      File: truelen.f
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
! INTEGER FUNCTION TRUELEN(CADENA)
!
! Input: CADENA
! Output: TRUELEN (function)
!
! Return the position of the last non-blank character in CADENA (ignoring also
! control characters with ASCII value < 32)
!
! CHARACTER*(*) CADENA -> input character string
!
!Comment
!------------------------------------------------------------------------------
        INTEGER FUNCTION TRUELEN(CADENA)
        IMPLICIT NONE
        CHARACTER*(*) CADENA
!
        INTEGER I,L
!------------------------------------------------------------------------------
        L=LEN(CADENA)
!
        DO I=L,1,-1
          IF(ICHAR(CADENA(I:I)).GT.32)THEN
            TRUELEN=I
            RETURN
          END IF
        END DO
        TRUELEN=0
        END
