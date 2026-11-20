!------------------------------------------------------------------------------
! Version 20-November-2026                                   File: fpercent.f90
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
! REAL FUNCTION FPERCENT(N,X,PERCENTILE)
!
! Input: N,X,PERCENTILE
! Output: FPERCENT (function)
!
! Calculate a fixed percentile of X(N)
!
! INTEGER N -> no. of elements
! REAL    X(N) -> input matrix
! REAL    PERCENTILE -> percentile to be computed
!
!Comment
!------------------------------------------------------------------------------
        REAL FUNCTION FPERCENT(N,X,PERCENTILE)
        IMPLICIT NONE
        INTEGER N
        REAL X(N),PERCENTILE
!
        INTEGER NMAX
        PARAMETER (NMAX=10000)
!
        INTEGER I
        INTEGER N1,N2
        INTEGER IFLAG
        REAL XSORTED(NMAX),XNUM(NMAX)
        REAL FRACTION
!------------------------------------------------------------------------------
        IF(N.GT.NMAX)THEN
          WRITE(*,101) 'FATAL ERROR in subroutine FPERCENT:'
          WRITE(*,101) 'N.GT.NMAX.'
          WRITE(*,100,ADVANCE='NO') 'Press <CR> to continue...'
          READ(*,*)
          FPERCENT=0.
          RETURN
        END IF
!
        IF(N.EQ.1)THEN
          FPERCENT=X(1)
          RETURN
        END IF
!
        DO I=1,N
          XNUM(I)=REAL(I)
          XSORTED(I)=X(I)
        END DO
!
        IF((PERCENTILE.LT.0.0).OR.(PERCENTILE.GT.100.0))THEN
          WRITE(*,101) 'FATAL ERROR in subroutine FPERCENT:'
          WRITE(*,100,ADVANCE='NO') 'PERCENTILE= '
          WRITE(*,*) PERCENTILE
          WRITE(*,101) 'PERCENTILE out of range.'
          WRITE(*,100,ADVANCE='NO') 'Press <CR> to continue...'
          READ(*,*)
          FPERCENT=0.
          RETURN
        END IF
!
        CALL ORDENA1F(N,XSORTED)
!
        FRACTION=1.+REAL(N-1)*PERCENTILE/100.
        FPERCENT=LININTERP(N,XNUM,XSORTED,FRACTION,IFLAG,N1,N2)
!
100     FORMAT(A)
101     FORMAT(A)
        END
