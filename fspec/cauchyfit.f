C------------------------------------------------------------------------------
C Version 13-October-2007                                     File: cauchyfit.f
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
C SUBROUTINE CAUCHYFIT(X0,SIGMA,AMP,EEX0,EESIGMA,EEAMP,YRMSTOL)
C
C Input: YRMSTOL
C Input (COMMON): NP,XF,YF
C Output: X0,SIGMA,AMP,EEX0,EESIGMA,EEAMP
C
C Fit numerically a Cauchy function (using DOWNHILL):
C Y=AMP/[SIGMA^2+(X-X0)^2]
C
C REAL X0 -> center of the fitted Cauchy function
C REAL SIGMA -> sigma value of the fitted Cauchy function
C REAL AMP -> maximum of the fitted Cauchy function
C REAL EEX0,EESIGMA,EEAMP -> errors in X0,SIGMA,AMP (rms from DOWNHILL)
C REAL YRMSTOL -> stopping criterion for DOWNHILL
C
Comment
C------------------------------------------------------------------------------
        SUBROUTINE CAUCHYFIT(X0,SIGMA,AMP,EEX0,EESIGMA,EEAMP,YRMSTOL)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INTEGER NP
        REAL XF(NCMAX),YF(NCMAX)
        COMMON/BLKFITG1/NP
        COMMON/BLKFITG2/XF,YF
        REAL X0,SIGMA,AMP
        REAL EEX0,EESIGMA,EEAMP
        REAL YRMSTOL
C
        INTEGER I
        REAL X0INI
        REAL XMIN,XMAX
        REAL A(3),CHISQR
        DOUBLE PRECISION MEAN,DISPER
C
        INTEGER NDIM,NEVAL
        REAL XX0(3),DXX0(3)
        REAL XX(3),DXX(3)
        EXTERNAL FUNKCAUCHY
        REAL FUNKCAUCHY
        COMMON/BLKX0INI/X0INI
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        AMP=0.
        EEAMP=0.
        X0=0.
        EEX0=0.
        SIGMA=0.
        EESIGMA=0.
C
        IF(NP.LT.3)THEN
          WRITE(*,101)'FATAL ERROR: in subroutine CAUCHYFIT.'
          WRITE(*,101)' No. of points for fit < 3'
          STOP
        END IF
C
C como primera estimacion del centro de la funcion calculamos el ajuste
C a una parabola y determinamos el maximo/minimo
        CALL POLFIT(XF,YF,YF,NP,3,0,A,CHISQR)
        X0INI=-REAL(A(2)/(2*A(3)))
        X0=0.              !amoeba supondra que el maximo esta alrededor de x=0
C
C Como estimacion de la amplitud tomamos 3 veces el valor de la dispersion
C de los datos alrededor de la media
        MEAN=0.D0
        DO I=1,NP
          MEAN=MEAN+DBLE(YF(I))
        END DO
        MEAN=MEAN/DBLE(NP)
        DISPER=0.D0
        DO I=1,NP
          DISPER=DISPER+(DBLE(YF(I))-MEAN)*(DBLE(YF(I))-MEAN)
        END DO
        DISPER=DSQRT(DISPER/DBLE(NP-1))
        AMP=3.*REAL(DISPER)
C si la parabola tiene un minimo (derivada segunda positiva), la 
C amplitud es negativa
        IF(A(3).GT.0.D0) AMP=-AMP
C
C Como estimacion del valor de SIGMA tomamos un cuarto del recorrido en la
C variable X
        XMIN=XF(1)
        XMAX=XMIN
        DO I=2,NP
          IF(XF(I).LT.XMIN) XMIN=XF(I)
          IF(XF(I).GT.XMAX) XMAX=XF(I)
        END DO
        SIGMA=(XMAX-XMIN)/4.
C------------------------------------------------------------------------------
C Con estas estimaciones ya podemos calcular los parametros de entrada de
C la subrutina DOWNHILL
        NDIM=3
        NDIM=3
        XX0(1)=AMP
        XX0(2)=X0
        XX0(3)=SIGMA
        IF(AMP.EQ.0.)THEN
          DXX0(1)=0.1
        ELSE
          DXX0(1)=0.1*AMP
        END IF
        DXX0(2)=1.
        IF(SIGMA.EQ.0.)THEN
          DXX0(3)=0.1
        ELSE
          DXX0(3)=0.1*SIGMA
        END IF
C
        CALL DOWNHILL(NDIM,XX0,DXX0,FUNKCAUCHY,1.0,0.5,2.0,YRMSTOL,
     +   XX,DXX,NEVAL)
C
        AMP=XX(1)
        X0=XX(2)+X0INI
        SIGMA=XX(3)
        EEAMP=DXX(1)
        EEX0=DXX(2)
        EESIGMA=DXX(3)
C------------------------------------------------------------------------------
101     FORMAT(A)
        END
C
C******************************************************************************
C
C XX(1)=AMP, XX(2)=X0, XX(3)=SIGMA
        REAL FUNCTION FUNKCAUCHY(XX)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        REAL XX(3)
C
        INTEGER I
        INTEGER NP
        REAL XF(NCMAX),YF(NCMAX)
        COMMON/BLKFITG1/NP
        COMMON/BLKFITG2/XF,YF
        COMMON/BLKX0INI/X0INI
        REAL FF,X0INI
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
C Introducimos X0INI para que el ajuste se realice alrededor de x=0. De esta
C forma el ajuste es mejor
        FUNKCAUCHY=0.
        DO I=1,NP
          FF=XX(1)/(XX(3)*XX(3)+(XF(I)-XX(2)-X0INI)*(XF(I)-XX(2)-X0INI))
          FUNKCAUCHY=FUNKCAUCHY+(YF(I)-FF)*(YF(I)-FF)
        END DO
        END
