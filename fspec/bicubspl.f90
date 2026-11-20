!------------------------------------------------------------------------------
! Version 20-November-2026                                   File: bicubspl.f90
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
! SUBROUTINE BICUBSPL(X,Y,Z,NX,NY,NXDIM,NYDIM,AA,BB,CC)
!
! Input: X,Y,Z,NX,NY,NXDIM,NYDIM
! Output: AA,BB,CC
!
! This subroutine computes the spline coefficients required by the subroutine
! BICUBSPLX to compute a surface through bicubic spline interpolation.
!
! REAL X(NX) -> X-values to be fitted
! REAL Y(NY) -> Y-values to be fitted
! REAL Z(NX,NY) -> Z-values to be fitted
! INTEGER NX -> logical dimension of X and logical first dimension of Z
! INTEGER NY -> logical dimension of Y and logical second dimension of Z
! INTEGER NXDIM -> physical dimension of X and physical first dimension of Z
! INTEGER NYDIM -> pyshical dimension of Y and physical second dimension of Z
! REAL AA(NX,NY) -> spline coefficients of the one-dimensional cubic spline
!                   fits to the rows of Z
! REAL BB(NX,NY) -> spline coefficients of the one-dimensional cubic spline
!                   fits to the rows of Z
! REAL CC(NX,NY) -> spline coefficients of the one-dimensional cubic spline
!                   fits to the rows of Z
!
!Comment
!------------------------------------------------------------------------------
!?      SUBROUTINE BICUBSPL(X,Y,Z,NX,NY,NXDIM,NYDIM,AA,BB,CC)
        SUBROUTINE BICUBSPL(Y,Z,NX,NY,NXDIM,NYDIM,AA,BB,CC)
        IMPLICIT NONE
!
        INTEGER NX,NY,NXDIM,NYDIM
!?      REAL X(NX)
        REAL Y(NY)
        REAL Z(NXDIM,NYDIM)
        REAL AA(NXDIM,NYDIM),BB(NXDIM,NYDIM),CC(NXDIM,NYDIM)
! local parameters
        INTEGER NMAX
        PARAMETER (NMAX=100)
! local variables
        INTEGER I,J
        REAL YLOCAL(NMAX)
        REAL S(NMAX),A(NMAX),B(NMAX),C(NMAX)
!------------------------------------------------------------------------------
        IF(NY.GT.NMAX)THEN
          WRITE(*,100,ADVANCE='NO')&
     &       'FATAL ERROR in subroutine BICUBSPL: '
          WRITE(*,100,ADVANCE='NO')' no. of data points too large, NY: '
          WRITE(*,*)NY
          STOP
        END IF
!------------------------------------------------------------------------------
        DO I=1,NX
! tomamos una fila completa y la introducimos en YLOCAL
          DO J=1,NY
            YLOCAL(J)=Z(I,J)
          END DO
! ajustamos splines unidimensionales con IMODE=1 (splines naturales)
          CALL CUBSPL(Y,YLOCAL,NY,1,S,A,B,C)
          DO J=1,NY
            AA(I,J)=A(J)
            BB(I,J)=B(J)
            CC(I,J)=C(J)
          END DO
        END DO
!------------------------------------------------------------------------------
100     FORMAT(A)
        END
