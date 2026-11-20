!------------------------------------------------------------------------------
! Version 20-November-2026                                    File: rmblank.f90
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
! SUBROUTINE RMBLANK(C1,C2,L)
!
! Input: C1
! Output: C2,L
!
! Remove blanks in character string C1, returning C2 with a true length L
!
! CHARACTER*(*) C1 -> input character string
! CHARACTER*(*) C2 -> output character string (C1 without blanks)
! INTEGER       L -> true len of C2
!
!Comment
!------------------------------------------------------------------------------
        SUBROUTINE RMBLANK(C1,C2,L)
        IMPLICIT NONE
        INTEGER L
        CHARACTER*(*) C1,C2
!
        INTEGER I,K
!------------------------------------------------------------------------------
        K=0
        DO I=1,LEN(C1)
          IF(C1(I:I).NE.CHAR(32))THEN
            K=K+1
            C2(K:K)=C1(I:I)
          END IF
        END DO
        L=K
        END
