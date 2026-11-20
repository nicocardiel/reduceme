!------------------------------------------------------------------------------
! Version 20-November-2026                                  File: pseudofit.f90
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
! SUBROUTINE PSEUDOFIT(XF,YF,NF,NTERMS,YRMSTOL,WEIGHT,POWER,LUP,A)
!
! Input: XF,YF,NF,NTERMS,YRMSTOL,WEIGHT,POWER,LUP
! Output: A
!
! Calculate the polynomial fit to the upper/lower side of a set of data
! points.
!
! REAL XF(NF),YF(NF) -> data points to be fitted
! INTEGER NF -> number of data points
! INTEGER NTERMS -> number of coeffcients
! REAL YRMSTOL -> stopping criterion for DOWNHILL
! REAL WEIGHT -> weighting factor to enhance one side of the fit
! REAL POWER -> power to be used to compute distances
! LOGICAL LUP -> .TRUE.: fit upper side
!                .FALSE.: fit lower side
! REAL A(NTERMS) -> fitted coefficients
!
!Comment
!------------------------------------------------------------------------------
        SUBROUTINE PSEUDOFIT(XF,YF,NF,NTERMS,YRMSTOL,WEIGHT,POWER,LUP,A)
        USE REDLIB_M
        IMPLICIT NONE
        INTEGER NF
        REAL XF(NF),YF(NF)
        INTEGER NTERMS
        REAL YRMSTOL
        REAL WEIGHT
        REAL POWER
        LOGICAL LUP
        REAL A(20)
!
        INCLUDE 'redlib.inc'
!
        EXTERNAL YFUNK_PSEUDO
        REAL YFUNK_PSEUDO
!
        INTEGER NNF,NNTERMS
        INTEGER J,K
        INTEGER NEVAL
        REAL WWEIGHT
        REAL PPOWER
        !las siguientes variables se dimensionan a NMAXFFT porque, en
        !general, este parametro sera mayor que NCMAX
        REAL XXF(NMAXFFT),YYF(NMAXFFT)
        REAL SIGMAY(NMAXFFT)
        REAL X0(20),DX0(20),X(20),DX(20)
        REAL CHISQR
        LOGICAL LLUP
!
        COMMON/BLKFUNKPSEUDO0/NNF,NNTERMS
        COMMON/BLKFUNKPSEUDO1/XXF,YYF
        COMMON/BLKFUNKPSEUDO2/WWEIGHT,PPOWER
        COMMON/BLKFUNKPSEUDO3/LLUP
!------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
! protecciones
        IF(NF.GT.NMAXFFT)THEN
          WRITE(*,100,ADVANCE='NO') 'NF, NMAXFFT: '
          WRITE(*,*) NF,NMAXFFT
          STOP 'FATAL ERROR: NF.GT.NMAXFFT in PSEUDOFIT.'
        END IF
        IF(NTERMS.GT.20)THEN
          WRITE(*,100,ADVANCE='NO') 'NTERMS: '
          WRITE(*,*) NTERMS
          STOP 'FATAL ERROR: NTERMS.GT.20 in PSEUDOFIT.'
        END IF
! inicializacion (duplicamos los argumentos de entrada de la subrutina para 
! poder pasar la información mediante COMMON blocks a la funcion a minimizar)
        NNF=NF
        NNTERMS=NTERMS
        WWEIGHT=WEIGHT
        PPOWER=POWER
        LLUP=LUP
        DO J=1,NF
          XXF(J)=XF(J)
          YYF(J)=YF(J)
          SIGMAY(J)=0.0
        END DO
!------------------------------------------------------------------------------
! Primero hacemos un ajuste tradicional para obtener una primera estimacion
        CALL POLFIT(XF,YF,SIGMAY,NF,NTERMS,0,A,CHISQR)
!------------------------------------------------------------------------------
! Usamos DOWNHILL para calcular el ajuste final
        DO K=1,NTERMS
          X0(K)=A(K)
          IF(A(K).NE.0.0)THEN
            DX0(K)=0.01*A(K)
          ELSE
            DX0(K)=1.0
          END IF
        END DO
        CALL DOWNHILL(NTERMS,X0,DX0,YFUNK_PSEUDO,1.0,0.5,2.0,YRMSTOL,&
     &   X,DX,NEVAL)
        DO K=1,NTERMS
          A(K)=X(K)
        END DO
!------------------------------------------------------------------------------
100     FORMAT(A)
        END
!
!******************************************************************************
!
        REAL FUNCTION YFUNK_PSEUDO(X)
        USE REDLIB_M
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        REAL X(20)
!
        REAL FPOLY
!
        INTEGER J
        INTEGER NF,NTERMS
        INTEGER NDEG
        REAL XF(NMAXFFT),YF(NMAXFFT)
        REAL WEIGHT,POWER
        REAL W1,W2
        REAL YPOL
        DOUBLE PRECISION DSUM
        LOGICAL LUP
!
        COMMON/BLKFUNKPSEUDO0/NF,NTERMS
        COMMON/BLKFUNKPSEUDO1/XF,YF
        COMMON/BLKFUNKPSEUDO2/WEIGHT,POWER
        COMMON/BLKFUNKPSEUDO3/LUP
!------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
!------------------------------------------------------------------------------
        DSUM=0.D0
        NDEG=NTERMS-1
!
        IF(LUP)THEN
          W1=1.0
          W2=WEIGHT
        ELSE
          W1=WEIGHT
          W2=1.0
        END IF
!
        DO J=1,NF
          YPOL=FPOLY(NDEG,X,XF(J))
          IF(YPOL.GE.YF(J))THEN
            DSUM=DSUM+W1*((YPOL-YF(J))**POWER)
          ELSE
            DSUM=DSUM+W2*((YF(J)-YPOL)**POWER)
          END IF
        END DO
!
        YFUNK_PSEUDO=REAL(DSUM)
!------------------------------------------------------------------------------
        END
