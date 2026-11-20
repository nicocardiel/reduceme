!------------------------------------------------------------------------------
! Version 20-November-2026                                File: findmmlmask.f90
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
! SUBROUTINE FINDMMLMASK(N,N1,N2,X,LMASK,XMIN,XMAX)
!
! Input: N,N1,N2,X,LMASK
! Output: XMIN,XMAX
!
! Return the maximum and minimum value of matrix X of N elements (in the
! range from N1 to N2 exclusively, making use of the boolean mask LMASK)
!
! INTEGER N -> no. of elements of matrix X
! INTEGER N1 -> first element of X() to search minimum/maximum
! INTEGER N2 -> last element of X() to search minimum/maximum
! REAL    X(N) -> data matrix
! LOGICAL LMASK(N) -> boolean mask
! REAL    XMIN -> minimum value of X()
! REAL    XMAX -> maximum value of X()
!
!Comment
!------------------------------------------------------------------------------
        SUBROUTINE FINDMMLMASK(N,N1,N2,X,LMASK,XMIN,XMAX)
        IMPLICIT NONE
        INTEGER N,N1,N2
        REAL X(N)
        LOGICAL LMASK(N)
        REAL XMIN,XMAX
!
        INTEGER I
        LOGICAL LOOP,LANY
!------------------------------------------------------------------------------
        IF((N1.LT.1).OR.(N2.GT.N).OR.(N2.LT.N1))THEN
          WRITE(*,101)'ERROR: limits out of range in FINDMML'
          WRITE(*,101)'=> Returned values: XMIN = XMAX = 0'
          XMIN=0.
          XMAX=0.
          RETURN
        END IF
!
        I=N1-1
        LANY=.FALSE.
        LOOP=.TRUE.
        DO WHILE(LOOP)
          I=I+1
          IF(I.GT.N2)THEN
            LOOP=.FALSE.
          ELSE
            IF(LMASK(I))THEN
              LANY=.TRUE.
              LOOP=.FALSE.
            END IF
          END IF
        END DO
!
        IF(.NOT.LANY)THEN
          WRITE(*,101)'ERROR: no data available in FINDMMLMASK'
          WRITE(*,101)'=> Returned values: XMIN = XMAX = 0'
          XMIN=0.
          XMAX=0.
          RETURN
        END IF
!
        XMIN=X(I)
        XMAX=XMIN
        DO I=N1,N2
          IF(LMASK(I))THEN
            IF(X(I).LT.XMIN) XMIN=X(I)
            IF(X(I).GT.XMAX) XMAX=X(I)
          END IF
        END DO
101     FORMAT(A)
        END
