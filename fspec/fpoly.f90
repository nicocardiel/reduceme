C------------------------------------------------------------------------------
C Version 26-November-2007                                        File: fpoly.f
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
C REAL FUNCTION FPOLY(NDEG,COEFF,X)
C
C Input: NDEG,COEFF,X
C Output: FPOLY (function)
C
C Evaluate the polynomial of degree NDEG and coefficients COEFF at X.
C
C INTEGER NDEG -> polynomial degree
C REAL    COEFF(NDEG+1) -> polynomial coefficients
C REAL    X -> abscissa at which the polynomial is going to be evaluated
C
Comment
C------------------------------------------------------------------------------
        REAL FUNCTION FPOLY(NDEG,COEFF,X)
        IMPLICIT NONE
        INTEGER NDEG
        REAL COEFF(NDEG+1)
        REAL X
C
        INTEGER K
        DOUBLE PRECISION DSUM
C------------------------------------------------------------------------------
        DSUM=DBLE(COEFF(NDEG+1))
        IF(NDEG.GT.0)THEN
          DO K=NDEG,1,-1
            DSUM=DSUM*DBLE(X)+DBLE(COEFF(K))
          END DO
        END IF
C
        FPOLY=REAL(DSUM)
        END
