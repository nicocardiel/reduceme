C------------------------------------------------------------------------------
C Version 8-October-1998                                         File: fmean0.f
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
C REAL FUNCTION FMEAN0(N,X,SIGMA)
C
C Input: N,X
C Output: FMEAN0 (function),SIGMA
C
C Calculate the mean value of X(N) and its r.m.s.
C
C INTEGER N -> no. of elements
C REAL    X(N) -> input matrix
C REAL SIGMA -> r.m.s. around the mean value
C
Comment
C------------------------------------------------------------------------------
        REAL FUNCTION FMEAN0(N,X,SIGMA)
        IMPLICIT NONE
        INTEGER N
        REAL X(N)
        REAL SIGMA
C
        INTEGER I
        DOUBLE PRECISION DSUM
C------------------------------------------------------------------------------
        IF(N.LE.0) STOP 'FATAL ERROR: in function FMEAN0: N.LE.0'
        DSUM=0.D0
        DO I=1,N
          DSUM=DSUM+DBLE(X(I))
        END DO
        DSUM=DSUM/DBLE(N)
        FMEAN0=REAL(DSUM)
C
        IF(N.EQ.1)THEN
          SIGMA=0.
        ELSE
          DSUM=0.D0
          DO I=1,N
            DSUM=DSUM+(DBLE(X(I))-DBLE(FMEAN0))*
     +       (DBLE(X(I))-DBLE(FMEAN0))
          END DO
          DSUM=DSUM/DBLE(N-1)
          SIGMA=SQRT(DSUM)
        END IF
C
        END
