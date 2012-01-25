C------------------------------------------------------------------------------
C Version 18-June-1998                                         File: bicubspl.f
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
C SUBROUTINE BICUBSPL(X,Y,Z,NX,NY,NXDIM,NYDIM,AA,BB,CC)
C
C Input: X,Y,Z,NX,NY,NXDIM,NYDIM
C Output: AA,BB,CC
C
C This subroutine computes the spline coefficients required by the subroutine
C BICUBSPLX to compute a surface through bicubic spline interpolation.
C
C REAL X(NX) -> X-values to be fitted
C REAL Y(NY) -> Y-values to be fitted
C REAL Z(NX,NY) -> Z-values to be fitted
C INTEGER NX -> logical dimension of X and logical first dimension of Z
C INTEGER NY -> logical dimension of Y and logical second dimension of Z
C INTEGER NXDIM -> physical dimension of X and physical first dimension of Z
C INTEGER NYDIM -> pyshical dimension of Y and physical second dimension of Z
C REAL AA(NX,NY) -> spline coefficients of the one-dimensional cubic spline
C                   fits to the rows of Z
C REAL BB(NX,NY) -> spline coefficients of the one-dimensional cubic spline
C                   fits to the rows of Z
C REAL CC(NX,NY) -> spline coefficients of the one-dimensional cubic spline
C                   fits to the rows of Z
C
Comment
C------------------------------------------------------------------------------
!?      SUBROUTINE BICUBSPL(X,Y,Z,NX,NY,NXDIM,NYDIM,AA,BB,CC)
        SUBROUTINE BICUBSPL(Y,Z,NX,NY,NXDIM,NYDIM,AA,BB,CC)
        IMPLICIT NONE
C
        INTEGER NX,NY,NXDIM,NYDIM
!?      REAL X(NX)
        REAL Y(NY)
        REAL Z(NXDIM,NYDIM)
        REAL AA(NXDIM,NYDIM),BB(NXDIM,NYDIM),CC(NXDIM,NYDIM)
C local parameters
        INTEGER NMAX
        PARAMETER (NMAX=100)
C local variables
        INTEGER I,J
        REAL YLOCAL(NMAX)
        REAL S(NMAX),A(NMAX),B(NMAX),C(NMAX)
C------------------------------------------------------------------------------
        IF(NY.GT.NMAX)THEN
          WRITE(*,100)'FATAL ERROR in subroutine BICUBSPL: '
          WRITE(*,100)' no. of data points too large, NY: '
          WRITE(*,*)NY
          STOP
        END IF
C------------------------------------------------------------------------------
        DO I=1,NX
C tomamos una fila completa y la introducimos en YLOCAL
          DO J=1,NY
            YLOCAL(J)=Z(I,J)
          END DO
C ajustamos splines unidimensionales con IMODE=1 (splines naturales)
          CALL CUBSPL(Y,YLOCAL,NY,1,S,A,B,C)
          DO J=1,NY
            AA(I,J)=A(J)
            BB(I,J)=B(J)
            CC(I,J)=C(J)
          END DO
        END DO
C------------------------------------------------------------------------------
100     FORMAT(A,$)
        END
