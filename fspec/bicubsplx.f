C------------------------------------------------------------------------------
C Version 18-June-1998                                        File: bicubsplx.f
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
C SUBROUTINE BICUBSPLX(X,Y,Z,NX,NY,NXDIM,NYDIM,AA,BB,CC,X0,Y0,Z0)
C
C Input: X,Y,Z,NX,NY,NXDIM,NYDIM,AA,BB,CC,X0,Y0
C Output: Z0
C
C This subroutine computes the bicubic spline Z0 at X0,Y0, using the spline
C coefficients computed with BICUBSPL.
C
C REAL X(NX) -> X-values fitted with BICUBSPL
C REAL Y(NY) -> Y-values fitted with BICUBSPL
C REAL Z(NX,NY) -> Z-values fitted with BICUBSPL
C INTEGER NX -> logical dimension of X and logical first dimension of Z
C INTEGER NY -> logical dimension of Y and logical second dimension of Z
C INTEGER NXDIM -> physical dimension of X and physical first dimension of Z
C INTEGER NYDIM -> pyshical dimension of Y and physical second dimension of Z
C REAL AA(NX,NY) -> coefficients of the one-dimensional cubic spline fits
C                  to the rows of Z computed with BICUBSPL
C REAL BB(NX,NY) -> coefficients of the one-dimensional cubic spline fits
C                  to the rows of Z computed with BICUBSPL
C REAL CC(NX,NY) -> coefficients of the one-dimensional cubic spline fits
C                  to the rows of Z computed with BICUBSPL
C REAL X0 -> X-value where the bicubic spline will be evaluated
C REAL Y0 -> Y-value where the bicubic spline will be evaluated
C REAL Z0 -> bicubic spline value at X0,Y0
C
Comment
C------------------------------------------------------------------------------
        SUBROUTINE BICUBSPLX(X,Y,Z,NX,NY,NXDIM,NYDIM,AA,BB,CC,X0,Y0,Z0)
        IMPLICIT NONE
C
        INTEGER NX,NY,NXDIM,NYDIM
        REAL X(NX),Y(NY)
        REAL Z(NXDIM,NYDIM)
        REAL AA(NXDIM,NYDIM),BB(NXDIM,NYDIM),CC(NXDIM,NYDIM)
        REAL X0,Y0,Z0
C local parameters
        INTEGER NMAX
        PARAMETER (NMAX=100)
C local variables
        INTEGER I,J,I0
        REAL A(NMAX),B(NMAX),C(NMAX),S(NMAX)
        REAL YLOCAL(NMAX),YYLOCAL(NMAX)
C------------------------------------------------------------------------------
        I0=1
        DO I=1,NX
          DO J=1,NY
            YLOCAL(J)=Z(I,J)
            A(J)=AA(I,J)
            B(J)=BB(I,J)
            C(J)=CC(I,J)
          END DO
          CALL CUBSPLX(Y,YLOCAL,A,B,C,NY,I0,Y0,YYLOCAL(I))
        END DO
        I0=1
        CALL CUBSPL(X,YYLOCAL,NX,1,S,A,B,C)                            !IMODE=1
        CALL CUBSPLX(X,YYLOCAL,A,B,C,NX,I0,X0,Z0)
C
        END
