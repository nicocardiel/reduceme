!------------------------------------------------------------------------------
! Version 20-November-2026                                   File: lagrange.f90
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
! SUBROUTINE LAGRANGE(N,X,Y,X0,Y0)
!
! Input: N,X,Y,X0
! Output: Y0
!
! Calculate the Lagrangian polynomial and evaluate such polynomial at X=X0. 
! We do not assume uniform spacing between the x-values, nor do we need the 
! x-values arranged in a particular order. However, the x-values must all be 
! distinct. We follow the algorithm described by B.P. Demidovich and I.A. Maron
! in Calculo Numerico Fundamental, Paraninfo 1988, pag. 593.
!
! INTEGER N -> no. of elements
! REAL    X(N) -> input x-matrix
! REAL    Y(N) -> input y-matrix
! REAL    X0 -> x-value where the Lagrangian polynomial is evaluated
! REAL    Y0 -> polynomial value at X=X0
!
!Comment
!------------------------------------------------------------------------------
        SUBROUTINE LAGRANGE(N,X,Y,X0,Y0)
        IMPLICIT NONE
!
        INTEGER N
        REAL X(N),Y(N)
        REAL X0,Y0
!
        INTEGER NMAX
        PARAMETER(NMAX=100)                  !numero maximo de puntos permitido
!
        INTEGER I,J
        REAL XX(NMAX),XX0
        REAL XMIN,XMAX
        REAL CX1,CX2
        REAL D(NMAX),DIAGONAL,SUM
!------------------------------------------------------------------------------
        IF(N.LT.1)THEN
          WRITE(*,100,ADVANCE='NO')&
     &       'FATAL ERROR: in subroutine LAGRANGE.'
          WRITE(*,101)'Number of points too small.'
          STOP
        END IF
        IF(N.GT.NMAX)THEN
          WRITE(*,100,ADVANCE='NO')&
     &       'FATAL ERROR: in subroutine LAGRANGE.'
          WRITE(*,101)'Number of points too large.'
          STOP
        END IF
!------------------------------------------------------------------------------
! verificamos si el resultado es inmediato
        DO I=1,N
          IF(X(I).EQ.X0)THEN
            Y0=Y(I)
            RETURN
          END IF
        END DO
!------------------------------------------------------------------------------
! cambiamos el recorrido de la variable X al intervalo [-1,+1]
        XMIN=X(1)
        XMAX=XMIN
        DO I=2,N
          IF(X(I).LT.XMIN) XMIN=X(I)
          IF(X(I).GT.XMAX) XMAX=X(I)
        END DO
        IF(XMIN.EQ.XMAX)THEN
          CX1=1.
          CX2=0.
        ELSE
          CX1=2./(XMAX-XMIN)
          CX2=(XMAX+XMIN)/(XMAX-XMIN)
        END IF
        DO I=1,N
          XX(I)=X(I)*CX1-CX2
        END DO
        XX0=X0*CX1-CX2
!------------------------------------------------------------------------------
        DIAGONAL=1.
        DO I=1,N
          DIAGONAL=DIAGONAL*(XX0-XX(I))
        END DO
!------------------------------------------------------------------------------
        SUM=0.
        DO I=1,N
          D(I)=1.
          DO J=1,N
            IF(J.EQ.I)THEN
              D(I)=D(I)*(XX0-XX(I))
            ELSE
              D(I)=D(I)*(XX(I)-XX(J))
            END IF
          END DO
          SUM=SUM+Y(I)/D(I)
        END DO
!------------------------------------------------------------------------------
        Y0=DIAGONAL*SUM
!------------------------------------------------------------------------------
100     FORMAT(A)
101     FORMAT(A)
        END
