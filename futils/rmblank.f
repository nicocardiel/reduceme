C------------------------------------------------------------------------------
C Version 25-November-1996                                      File: rmblank.f
C------------------------------------------------------------------------------
C Copyright N. Cardiel & J. Gorgas, Departamento de Astrofisica
C Universidad Complutense de Madrid, 28040-Madrid, Spain
C E-mail: ncl@astrax.fis.ucm.es or fjg@astrax.fis.ucm.es
C------------------------------------------------------------------------------
C This routine is free software; you can redistribute it and/or modify it
C under the terms of the GNU General Public License as published by the Free
C Software Foundation; either version 2 of the License, or (at your option) any
C later version. See the file gnu-public-license.txt for details.
C------------------------------------------------------------------------------
Comment
C
C SUBROUTINE RMBLANK(C1,C2,L)
C
C Input: C1
C Output: C2,L
C
C Remove blanks in character string C1, returning C2 with a true length L
C
C CHARACTER*(*) C1 -> input character string
C CHARACTER*(*) C2 -> output character string (C1 without blanks)
C INTEGER       L -> true len of C2
C
Comment
C------------------------------------------------------------------------------
        SUBROUTINE RMBLANK(C1,C2,L)
        IMPLICIT NONE
        INTEGER L
        CHARACTER*(*) C1,C2
C
        INTEGER I,K
C------------------------------------------------------------------------------
        K=0
        DO I=1,LEN(C1)
          IF(C1(I:I).NE.CHAR(32))THEN
            K=K+1
            C2(K:K)=C1(I:I)
          END IF
        END DO
        L=K
        END
