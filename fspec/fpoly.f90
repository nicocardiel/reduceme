!------------------------------------------------------------------------------
! Version 20-November-2026                                      File: fpoly.f90
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
! REAL FUNCTION FPOLY(NDEG,COEFF,X)
!
! Input: NDEG,COEFF,X
! Output: FPOLY (function)
!
! Evaluate the polynomial of degree NDEG and coefficients COEFF at X.
!
! INTEGER NDEG -> polynomial degree
! REAL    COEFF(NDEG+1) -> polynomial coefficients
! REAL    X -> abscissa at which the polynomial is going to be evaluated
!
!Comment
!------------------------------------------------------------------------------
        REAL FUNCTION FPOLY(NDEG,COEFF,X)
        IMPLICIT NONE
        INTEGER NDEG
        REAL COEFF(NDEG+1)
        REAL X
!
        INTEGER K
        DOUBLE PRECISION DSUM
!------------------------------------------------------------------------------
        DSUM=DBLE(COEFF(NDEG+1))
        IF(NDEG.GT.0)THEN
          DO K=NDEG,1,-1
            DSUM=DSUM*DBLE(X)+DBLE(COEFF(K))
          END DO
        END IF
!
        FPOLY=REAL(DSUM)
        END
