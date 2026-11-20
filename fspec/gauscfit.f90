!------------------------------------------------------------------------------
! Version 20-November-2026                                     File: gauscfit.f
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
! SUBROUTINE GAUSCFIT(X0,SIGMA,AMP,Y0,EEX0,EESIGMA,EEAMP,EEY0,YRMSTOL)
!
! Input: YRMSTOL
! Input (COMMON): NP,XF,YF
! Output: X0,SIGMA,AMP,Y0,EEX0,EESIGMA,EEAM,EEY0
!
! Fit numerically a gaussian + constant (using DOWNHILL):
! Y=Y0+AMP*EXP[-((X-X0)^2/(2*SIGMA^2))]
!
! REAL X0 -> center of the fitted gaussian
! REAL SIGMA -> sigma value of the fitted gaussian
! REAL AMP -> maximum of the fitted gaussian
! REAL Y0 -> constant
! REAL EEX0,EESIGMA,EEAMP,EEY0 -> errors in X0,SIGMA,AMP,Y0 (rms from DOWNHILL)
! REAL YRMSTOL -> stopping criterion for DOWNHILL
!
!Comment
!------------------------------------------------------------------------------
        SUBROUTINE GAUSCFIT(X0,SIGMA,AMP,Y0,EEX0,EESIGMA,EEAMP,EEY0,&
     &   YRMSTOL)
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
!
        INTEGER I
        REAL XMIN,XMAX
        REAL A(3),CHISQR
        DOUBLE PRECISION MEAN,DISPER
!
        INTEGER NDIM,NEVAL
        REAL XX0(4),DXX0(4)
        REAL XX(4),DXX(4)
        EXTERNAL FUNKGC
        REAL FUNKGC
        COMMON/BLKX0INI/X0INI
!
!------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
!
        AMP=0.
        EEAMP=0.
        X0=0.
        EEX0=0.
        SIGMA=0.
        EESIGMA=0.
        Y0=0.
        EEY0=0.
!
        IF(NP.LT.4)THEN
          WRITE(*,101)'FATAL ERROR: in subroutine GAUSCFIT.'
          WRITE(*,101)' No. of points for fit < 4'
          STOP
        END IF
        IF(NP.GT.NCMAX)THEN
          WRITE(*,101)'FATAL ERROR: in subroutine GAUSCFIT.'
          WRITE(*,100,ADVANCE='NO')'NP, NCMAX: '
          WRITE(*,*) NP,NCMAX
          WRITE(*,101)' No. of points for fit too large'
          STOP
        END IF
!
! como primera estimacion del centro de la gaussiana calculamos el ajuste
! a una parabola y determinamos el maximo/minimo
        CALL POLFIT(XF,YF,YF,NP,3,0,A,CHISQR)
        X0INI=-REAL(A(2)/(2*A(3)))
        X0=0.1       !DOWNHILL supondra que la gaussiana esta alrededor de cero
!
! Como estimacion de la amplitud tomamos 3 veces el valor de la dispersion
! de los datos alrededor de la media
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
! si la parabola tiene un minimo (derivada segunda positiva), la 
! amplitud es negativa
        IF(A(3).GT.0.D0) AMP=-AMP
!
! Como estimacion del valor de SIGMA tomamos un cuarto del recorrido en la
! variable X
        XMIN=XF(1)
        XMAX=XMIN
        DO I=2,NP
          IF(XF(I).LT.XMIN) XMIN=XF(I)
          IF(XF(I).GT.XMAX) XMAX=XF(I)
        END DO
        SIGMA=(XMAX-XMIN)/4.
!
! Como estimacion de la cte. tomamos el valor medio de los datos
        Y0=REAL(MEAN)
!cc        type*,amp,x0,sigma,y0
!------------------------------------------------------------------------------
! Con estas estimaciones ya podemos calcular los parametros de entrada de
! la subrutina DOWNHILL
        NDIM=4
        XX0(1)=AMP
        XX0(2)=X0
        XX0(3)=SIGMA
        XX0(4)=Y0
! Hay que evitar que los DXX sean nulos
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
!
        CALL DOWNHILL(NDIM,XX0,DXX0,FUNKGC,1.0,0.5,2.0,YRMSTOL,&
     &   XX,DXX,NEVAL)
!
        AMP=XX(1)
        X0=XX(2)+X0INI
        SIGMA=XX(3)
        Y0=XX(4)
        EEAMP=DXX(1)
        EEX0=DXX(2)
        EESIGMA=DXX(3)
        EEY0=DXX(4)
!------------------------------------------------------------------------------
100     FORMAT(A)
101     FORMAT(A)
        END
!
!******************************************************************************
!
! XX(1)=AMP, XX(2)=X0, XX(3)=SIGMA, XX(4)=Y0
        REAL FUNCTION FUNKGC(XX)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        REAL XX(4)
!
        INTEGER I
        INTEGER NP
        REAL XF(NCMAX),YF(NCMAX)
        COMMON/BLKFITG1/NP
        COMMON/BLKFITG2/XF,YF
        COMMON/BLKX0INI/X0INI
        DOUBLE PRECISION DFUNK
        REAL FF,X0INI
!------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
!
! Introducimos un offset de X0INI para que el ajuste se realiza alrededor
! de x=0 (funciona mucho mejor)
        DFUNK=0.D0
        DO I=1,NP
          FF=XX(4)+XX(1)*&
     &     EXP(-(((XF(I)-XX(2)-X0INI)*(XF(I)-XX(2)-X0INI))/&
     &     (2.*XX(3)*XX(3))))
          DFUNK=DFUNK+DBLE((YF(I)-FF)*(YF(I)-FF))
        END DO
        FUNKGC=REAL(DFUNK)
        END
