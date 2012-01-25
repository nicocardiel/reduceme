C------------------------------------------------------------------------------
C Version 26-April-1999                                        File: fpercent.f
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
C REAL FUNCTION FPERCENT(N,X,PERCENTILE)
C
C Input: N,X,PERCENTILE
C Output: FPERCENT (function)
C
C Calculate a fixed percentile of X(N)
C
C INTEGER N -> no. of elements
C REAL    X(N) -> input matrix
C REAL    PERCENTILE -> percentile to be computed
C
Comment
C------------------------------------------------------------------------------
        REAL FUNCTION FPERCENT(N,X,PERCENTILE)
        IMPLICIT NONE
        INTEGER N
        REAL X(N),PERCENTILE
C
        INTEGER NMAX
        PARAMETER (NMAX=10000)
C
        INTEGER I
        INTEGER N1,N2
        INTEGER IFLAG
        REAL XSORTED(NMAX),XNUM(NMAX)
        REAL FRACTION
        REAL LININTERP
C------------------------------------------------------------------------------
        IF(N.GT.NMAX)THEN
          WRITE(*,101) 'FATAL ERROR in subroutine FPERCENT:'
          WRITE(*,101) 'N.GT.NMAX.'
          WRITE(*,100) 'Press <CR> to continue...'
          READ(*,*)
          FPERCENT=0.
          RETURN
        END IF
C
        IF(N.EQ.1)THEN
          FPERCENT=X(1)
          RETURN
        END IF
C
        DO I=1,N
          XNUM(I)=REAL(I)
          XSORTED(I)=X(I)
        END DO
C
        IF((PERCENTILE.LT.0.0).OR.(PERCENTILE.GT.100.0))THEN
          WRITE(*,101) 'FATAL ERROR in subroutine FPERCENT:'
          WRITE(*,100) 'PERCENTILE= '
          WRITE(*,*) PERCENTILE
          WRITE(*,101) 'PERCENTILE out of range.'
          WRITE(*,100) 'Press <CR> to continue...'
          READ(*,*)
          FPERCENT=0.
          RETURN
        END IF
C
        CALL ORDENA1F(N,XSORTED)
C
        FRACTION=1.+REAL(N-1)*PERCENTILE/100.
        FPERCENT=LININTERP(N,XNUM,XSORTED,FRACTION,IFLAG,N1,N2)
C
100     FORMAT(A,$)
101     FORMAT(A)
        END
