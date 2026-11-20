!------------------------------------------------------------------------------
! Version 8-October-1998                                         File: fmean0.f
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
! REAL FUNCTION FMEAN0(N,X,SIGMA)
!
! Input: N,X
! Output: FMEAN0 (function),SIGMA
!
! Calculate the mean value of X(N) and its r.m.s.
!
! INTEGER N -> no. of elements
! REAL    X(N) -> input matrix
! REAL SIGMA -> r.m.s. around the mean value
!
!Comment
!------------------------------------------------------------------------------
        REAL FUNCTION FMEAN0(N,X,SIGMA)
        IMPLICIT NONE
        INTEGER N
        REAL X(N)
        REAL SIGMA
!
        INTEGER I
        DOUBLE PRECISION DSUM
!------------------------------------------------------------------------------
        IF(N.LE.0) STOP 'FATAL ERROR: in function FMEAN0: N.LE.0'
        DSUM=0.D0
        DO I=1,N
          DSUM=DSUM+DBLE(X(I))
        END DO
        DSUM=DSUM/DBLE(N)
        FMEAN0=REAL(DSUM)
!
        IF(N.EQ.1)THEN
          SIGMA=0.
        ELSE
          DSUM=0.D0
          DO I=1,N
            DSUM=DSUM+(DBLE(X(I))-DBLE(FMEAN0))*&
     &       (DBLE(X(I))-DBLE(FMEAN0))
          END DO
          DSUM=DSUM/DBLE(N-1)
          SIGMA=SQRT(DSUM)
        END IF
!
        END
