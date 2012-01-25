C------------------------------------------------------------------------------
C Version 26-November-2007                                    File: pseudofit.f
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
C SUBROUTINE PSEUDOFIT(XF,YF,NF,NTERMS,YRMSTOL,WEIGHT,POWER,LUP,A)
C
C Input: XF,YF,NF,NTERMS,YRMSTOL,WEIGHT,POWER,LUP
C Output: A
C
C Calculate the polynomial fit to the upper/lower side of a set of data
C points.
C
C REAL XF(NF),YF(NF) -> data points to be fitted
C INTEGER NF -> number of data points
C INTEGER NTERMS -> number of coeffcients
C REAL YRMSTOL -> stopping criterion for DOWNHILL
C REAL WEIGHT -> weighting factor to enhance one side of the fit
C REAL POWER -> power to be used to compute distances
C LOGICAL LUP -> .TRUE.: fit upper side
C                .FALSE.: fit lower side
C REAL A(NTERMS) -> fitted coefficients
C
Comment
C------------------------------------------------------------------------------
        SUBROUTINE PSEUDOFIT(XF,YF,NF,NTERMS,YRMSTOL,WEIGHT,POWER,LUP,A)
        IMPLICIT NONE
        INTEGER NF
        REAL XF(NF),YF(NF)
        INTEGER NTERMS
        REAL YRMSTOL
        REAL WEIGHT
        REAL POWER
        LOGICAL LUP
        REAL A(20)
C
        INCLUDE 'redlib.inc'
C
        EXTERNAL YFUNK_PSEUDO
        REAL YFUNK_PSEUDO
C
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
C
        COMMON/BLKFUNKPSEUDO0/NNF,NNTERMS
        COMMON/BLKFUNKPSEUDO1/XXF,YYF
        COMMON/BLKFUNKPSEUDO2/WWEIGHT,PPOWER
        COMMON/BLKFUNKPSEUDO3/LLUP
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C protecciones
        IF(NF.GT.NMAXFFT)THEN
          WRITE(*,100) 'NF, NMAXFFT: '
          WRITE(*,*) NF,NMAXFFT
          STOP 'FATAL ERROR: NF.GT.NMAXFFT in PSEUDOFIT.'
        END IF
        IF(NTERMS.GT.20)THEN
          WRITE(*,100) 'NTERMS: '
          WRITE(*,*) NTERMS
          STOP 'FATAL ERROR: NTERMS.GT.20 in PSEUDOFIT.'
        END IF
C inicializacion (duplicamos los argumentos de entrada de la subrutina para 
C poder pasar la información mediante COMMON blocks a la funcion a minimizar)
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
C------------------------------------------------------------------------------
C Primero hacemos un ajuste tradicional para obtener una primera estimacion
        CALL POLFIT(XF,YF,SIGMAY,NF,NTERMS,0,A,CHISQR)
C------------------------------------------------------------------------------
C Usamos DOWNHILL para calcular el ajuste final
        DO K=1,NTERMS
          X0(K)=A(K)
          IF(A(K).NE.0.0)THEN
            DX0(K)=0.01*A(K)
          ELSE
            DX0(K)=1.0
          END IF
        END DO
        CALL DOWNHILL(NTERMS,X0,DX0,YFUNK_PSEUDO,1.0,0.5,2.0,YRMSTOL,
     +   X,DX,NEVAL)
        DO K=1,NTERMS
          A(K)=X(K)
        END DO
C------------------------------------------------------------------------------
100     FORMAT(A,$)
        END
C
C******************************************************************************
C
        REAL FUNCTION YFUNK_PSEUDO(X)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        REAL X(20)
C
        REAL FPOLY
C
        INTEGER J
        INTEGER NF,NTERMS
        INTEGER NDEG
        REAL XF(NMAXFFT),YF(NMAXFFT)
        REAL WEIGHT,POWER
        REAL W1,W2
        REAL YPOL
        DOUBLE PRECISION DSUM
        LOGICAL LUP
C
        COMMON/BLKFUNKPSEUDO0/NF,NTERMS
        COMMON/BLKFUNKPSEUDO1/XF,YF
        COMMON/BLKFUNKPSEUDO2/WEIGHT,POWER
        COMMON/BLKFUNKPSEUDO3/LUP
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C------------------------------------------------------------------------------
        DSUM=0.D0
        NDEG=NTERMS-1
C
        IF(LUP)THEN
          W1=1.0
          W2=WEIGHT
        ELSE
          W1=WEIGHT
          W2=1.0
        END IF
C
        DO J=1,NF
          YPOL=FPOLY(NDEG,X,XF(J))
          IF(YPOL.GE.YF(J))THEN
            DSUM=DSUM+W1*((YPOL-YF(J))**POWER)
          ELSE
            DSUM=DSUM+W2*((YF(J)-YPOL)**POWER)
          END IF
        END DO
C
        YFUNK_PSEUDO=REAL(DSUM)
C------------------------------------------------------------------------------
        END
