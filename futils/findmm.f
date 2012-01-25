C------------------------------------------------------------------------------
C Version 25-November-1996                                       File: findmm.f
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
C SUBROUTINE FINDMM(N,X,XMIN,XMAX)
C
C Input: N,X
C Output: XMIN,XMAX
C
C Return the maximum and minimum value of matrix X of N elements
C
C INTEGER N -> no. of elements of matrix X
C REAL    X(N) -> data matrix
C REAL    XMIN -> minimum value of X()
C REAL    XMAX -> maximum value of X()
C
Comment
C------------------------------------------------------------------------------
        SUBROUTINE FINDMM(N,X,XMIN,XMAX)
        IMPLICIT NONE
        INTEGER N
        REAL X(N)
        REAL XMIN,XMAX
C
        INTEGER I
C------------------------------------------------------------------------------
        XMIN=X(1)
        XMAX=XMIN
        DO I=2,N
          IF(X(I).LT.XMIN) XMIN=X(I)
          IF(X(I).GT.XMAX) XMAX=X(I)
        END DO
        END
