C------------------------------------------------------------------------------
C Version 28-February-1997                                     File: gauscfit.f
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
C SUBROUTINE GAUSCFIT(X0,SIGMA,AMP,Y0,EEX0,EESIGMA,EEAMP,EEY0,YRMSTOL)
C
C Input: YRMSTOL
C Input (COMMON): NP,XF,YF
C Output: X0,SIGMA,AMP,Y0,EEX0,EESIGMA,EEAM,EEY0
C
C Fit numerically a gaussian + constant (using DOWNHILL):
C Y=Y0+AMP*EXP[-((X-X0)^2/(2*SIGMA^2))]
C
C REAL X0 -> center of the fitted gaussian
C REAL SIGMA -> sigma value of the fitted gaussian
C REAL AMP -> maximum of the fitted gaussian
C REAL Y0 -> constant
C REAL EEX0,EESIGMA,EEAMP,EEY0 -> errors in X0,SIGMA,AMP,Y0 (rms from DOWNHILL)
C REAL YRMSTOL -> stopping criterion for DOWNHILL
C
Comment
C------------------------------------------------------------------------------
        SUBROUTINE GAUSCFIT(X0,SIGMA,AMP,Y0,EEX0,EESIGMA,EEAMP,EEY0,
     +   YRMSTOL)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INTEGER NP
        REAL XF(NCMAX),YF(NCMAX)
        COMMON/BLKFITG1/NP
        COMMON/BLKFITG2/XF,YF
        REAL X0,SIGMA,AMP,Y0
        REAL EEX0,EESIGMA,EEAMP,EEY0
        REAL YRMSTOL
        REAL X0INI
C
        INTEGER I
        REAL XMIN,XMAX
        REAL A(3),CHISQR
        DOUBLE PRECISION MEAN,DISPER
C
        INTEGER NDIM,NEVAL
        REAL XX0(4),DXX0(4)
        REAL XX(4),DXX(4)
        EXTERNAL FUNKGC
        REAL FUNKGC
        COMMON/BLKX0INI/X0INI
C
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        AMP=0.
        EEAMP=0.
        X0=0.
        EEX0=0.
        SIGMA=0.
        EESIGMA=0.
        Y0=0.
        EEY0=0.
C
        IF(NP.LT.4)THEN
          WRITE(*,101)'FATAL ERROR: in subroutine GAUSCFIT.'
          WRITE(*,101)' No. of points for fit < 4'
          STOP
        END IF
        IF(NP.GT.NCMAX)THEN
          WRITE(*,101)'FATAL ERROR: in subroutine GAUSCFIT.'
          WRITE(*,100)'NP, NCMAX: '
          WRITE(*,*) NP,NCMAX
          WRITE(*,101)' No. of points for fit too large'
          STOP
        END IF
C
C como primera estimacion del centro de la gaussiana calculamos el ajuste
C a una parabola y determinamos el maximo/minimo
        CALL POLFIT(XF,YF,YF,NP,3,0,A,CHISQR)
        X0INI=-REAL(A(2)/(2*A(3)))
        X0=0.1       !DOWNHILL supondra que la gaussiana esta alrededor de cero
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
C
C Como estimacion de la cte. tomamos el valor medio de los datos
        Y0=REAL(MEAN)
ccc        type*,amp,x0,sigma,y0
C------------------------------------------------------------------------------
C Con estas estimaciones ya podemos calcular los parametros de entrada de
C la subrutina DOWNHILL
        NDIM=4
        XX0(1)=AMP
        XX0(2)=X0
        XX0(3)=SIGMA
        XX0(4)=Y0
C Hay que evitar que los DXX sean nulos
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
        DXX0(4)=DXX0(1)
C
        CALL DOWNHILL(NDIM,XX0,DXX0,FUNKGC,1.0,0.5,2.0,YRMSTOL,
     +   XX,DXX,NEVAL)
C
        AMP=XX(1)
        X0=XX(2)+X0INI
        SIGMA=XX(3)
        Y0=XX(4)
        EEAMP=DXX(1)
        EEX0=DXX(2)
        EESIGMA=DXX(3)
        EEY0=DXX(4)
C------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
        END
C
C******************************************************************************
C
C XX(1)=AMP, XX(2)=X0, XX(3)=SIGMA, XX(4)=Y0
        REAL FUNCTION FUNKGC(XX)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        REAL XX(4)
C
        INTEGER I
        INTEGER NP
        REAL XF(NCMAX),YF(NCMAX)
        COMMON/BLKFITG1/NP
        COMMON/BLKFITG2/XF,YF
        COMMON/BLKX0INI/X0INI
        DOUBLE PRECISION DFUNK
        REAL FF,X0INI
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
C Introducimos un offset de X0INI para que el ajuste se realiza alrededor
C de x=0 (funciona mucho mejor)
        DFUNK=0.D0
        DO I=1,NP
          FF=XX(4)+XX(1)*
     +     EXP(-(((XF(I)-XX(2)-X0INI)*(XF(I)-XX(2)-X0INI))/
     +     (2.*XX(3)*XX(3))))
          DFUNK=DFUNK+DBLE((YF(I)-FF)*(YF(I)-FF))
        END DO
        FUNKGC=REAL(DFUNK)
        END
