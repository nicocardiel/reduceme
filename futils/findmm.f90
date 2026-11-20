!------------------------------------------------------------------------------
! Version 20-November-2026                                     File: findmm.f90
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
! SUBROUTINE FINDMM(N,X,XMIN,XMAX)
!
! Input: N,X
! Output: XMIN,XMAX
!
! Return the maximum and minimum value of matrix X of N elements
!
! INTEGER N -> no. of elements of matrix X
! REAL    X(N) -> data matrix
! REAL    XMIN -> minimum value of X()
! REAL    XMAX -> maximum value of X()
!
!Comment
!------------------------------------------------------------------------------
        SUBROUTINE FINDMM(N,X,XMIN,XMAX)
        IMPLICIT NONE
        INTEGER N
        REAL X(N)
        REAL XMIN,XMAX
!
        INTEGER I
!------------------------------------------------------------------------------
        XMIN=X(1)
        XMAX=XMIN
        DO I=2,N
          IF(X(I).LT.XMIN) XMIN=X(I)
          IF(X(I).GT.XMAX) XMAX=X(I)
        END DO
        END
