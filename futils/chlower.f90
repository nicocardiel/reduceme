C------------------------------------------------------------------------------
C Version 25-November-1996                                      File: chlower.f
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
C SUBROUTINE CHLOWER(CADENA)
C
C Input: CADENA
C Output: CADENA
C
C Upper case characters in CADENA are transformed to lower case
C
C CHARACTER*(*) CADENA -> character string to be transformed
C
Comment
C------------------------------------------------------------------------------
        SUBROUTINE CHLOWER(CADENA)
        IMPLICIT NONE
        CHARACTER*(*) CADENA
C
        INTEGER I,N
C------------------------------------------------------------------------------
        DO I=1,LEN(CADENA)
          N=ICHAR(CADENA(I:I))
          IF((N.GE.65).AND.(N.LE.90)) CADENA(I:I)=CHAR(N+32)
        END DO
        END
