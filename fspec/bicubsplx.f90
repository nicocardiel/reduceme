!------------------------------------------------------------------------------
! Version 18-June-1998                                        File: bicubsplx.f
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
! SUBROUTINE BICUBSPLX(X,Y,Z,NX,NY,NXDIM,NYDIM,AA,BB,CC,X0,Y0,Z0)
!
! Input: X,Y,Z,NX,NY,NXDIM,NYDIM,AA,BB,CC,X0,Y0
! Output: Z0
!
! This subroutine computes the bicubic spline Z0 at X0,Y0, using the spline
! coefficients computed with BICUBSPL.
!
! REAL X(NX) -> X-values fitted with BICUBSPL
! REAL Y(NY) -> Y-values fitted with BICUBSPL
! REAL Z(NX,NY) -> Z-values fitted with BICUBSPL
! INTEGER NX -> logical dimension of X and logical first dimension of Z
! INTEGER NY -> logical dimension of Y and logical second dimension of Z
! INTEGER NXDIM -> physical dimension of X and physical first dimension of Z
! INTEGER NYDIM -> pyshical dimension of Y and physical second dimension of Z
! REAL AA(NX,NY) -> coefficients of the one-dimensional cubic spline fits
!                  to the rows of Z computed with BICUBSPL
! REAL BB(NX,NY) -> coefficients of the one-dimensional cubic spline fits
!                  to the rows of Z computed with BICUBSPL
! REAL CC(NX,NY) -> coefficients of the one-dimensional cubic spline fits
!                  to the rows of Z computed with BICUBSPL
! REAL X0 -> X-value where the bicubic spline will be evaluated
! REAL Y0 -> Y-value where the bicubic spline will be evaluated
! REAL Z0 -> bicubic spline value at X0,Y0
!
!Comment
!------------------------------------------------------------------------------
        SUBROUTINE BICUBSPLX(X,Y,Z,NX,NY,NXDIM,NYDIM,AA,BB,CC,X0,Y0,Z0)
        IMPLICIT NONE
!
        INTEGER NX,NY,NXDIM,NYDIM
        REAL X(NX),Y(NY)
        REAL Z(NXDIM,NYDIM)
        REAL AA(NXDIM,NYDIM),BB(NXDIM,NYDIM),CC(NXDIM,NYDIM)
        REAL X0,Y0,Z0
! local parameters
        INTEGER NMAX
        PARAMETER (NMAX=100)
! local variables
        INTEGER I,J,I0
        REAL A(NMAX),B(NMAX),C(NMAX),S(NMAX)
        REAL YLOCAL(NMAX),YYLOCAL(NMAX)
!------------------------------------------------------------------------------
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
!
        END
