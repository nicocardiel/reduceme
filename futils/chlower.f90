!------------------------------------------------------------------------------
! Version 20-November-2026                                    File: chlower.f90
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
! SUBROUTINE CHLOWER(CADENA)
!
! Input: CADENA
! Output: CADENA
!
! Upper case characters in CADENA are transformed to lower case
!
! CHARACTER*(*) CADENA -> character string to be transformed
!
!Comment
!------------------------------------------------------------------------------
        SUBROUTINE CHLOWER(CADENA)
        IMPLICIT NONE
        CHARACTER*(*) CADENA
!
        INTEGER I,N
!------------------------------------------------------------------------------
        DO I=1,LEN(CADENA)
          N=ICHAR(CADENA(I:I))
          IF((N.GE.65).AND.(N.LE.90)) CADENA(I:I)=CHAR(N+32)
        END DO
        END
