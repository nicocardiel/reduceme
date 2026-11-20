!------------------------------------------------------------------------------
! Version 20-November-2026                                 File: gauss2afit.f90
!------------------------------------------------------------------------------
!Comment
!
! SUBROUTINE GAUSS2AFIT(NPFIT,XFIT,YFIT,EYFIT,DELTAX,X0,SIGMA,AMP,
!                       EX0,ESIGMA,EAMP,EEX0,EESIGMA,EEAMP,
!                       YRMSTOL,NSIMUL)
!
! Input: NPFIT,XFIT,YFIT,EYFIT,DELTAX,YRMSTOL,NSIMUL
! Output: X0,SIGMA,AMP,EX0,ESIGMA,EAMP,EEX0,EESIGMA,EEAMP
!
! Numerical fit of 2 gaussians with the same width and area (using DOWNHILL),
! with a fixed separation given by DELTAX.
! Y=AMP*EXP[-((X-X0)^2/(2*SIGMA^2))]+AMP*EXP[-((X-X0-DELTAX)^2/(2*SIGMA^2))]
!
! INTEGER NPFIT -> number of points to be fitted
! REAL XFIT,YFIT,EYFIT -> x, y and error
! REAL DELTAX -> separation between gaussians
! REAL X0 -> center of the fitted gaussian
! REAL SIGMA -> sigma value of the fitted gaussian
! REAL AMP -> maximum of the fitted gaussian
! REAL EX0,ESIGMA,EAMP -> errors in X0,SIGMA,AMP (due to EYFIT --simulations--)
! REAL EEX0,EESIGMA,EEAMP -> errors in X0,SIGMA,AMP (rms from DOWNHILL)
! REAL YRMSTOL -> stopping criterion for DOWNHILL
! INTEGER NSIMUL -> number of simulations to compute errors
!
!Comment
!------------------------------------------------------------------------------
        SUBROUTINE GAUSS2AFIT(NPFIT,XFIT,YFIT,EYFIT,DELTAX,X0,SIGMA,AMP,&
     &   EX0,ESIGMA,EAMP,EEX0,EESIGMA,EEAMP,YRMSTOL,NSIMUL) 
        USE REDLIB_M
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INTEGER NPFIT
        REAL XFIT(NPFIT),YFIT(NPFIT),EYFIT(NPFIT)
        REAL DELTAX
        REAL X0,SIGMA,AMP,EX0,ESIGMA,EAMP,EEX0,EESIGMA,EEAMP
        REAL YRMSTOL
        INTEGER NSIMUL
!
        INTEGER NSIMULMAX
        PARAMETER(NSIMULMAX=1000)                !numero maximo de simulaciones
        REAL PI2
        PARAMETER(PI2=6.283185307)               !2 x pi
!
        INTEGER I,J,ISIMUL
        INTEGER NP
        INTEGER NSEED
        REAL X0INI
        REAL XMIN,XMAX
        REAL XF(NCMAX),YF(NCMAX),ERRYF
        REAL A(3),CHISQR
        REAL X0_SIMUL(NSIMULMAX)
        REAL SIGMA_SIMUL(NSIMULMAX)
        REAL AMP_SIMUL(NSIMULMAX)
        REAL DELTAX_
        REAL RANRED,R1,R2
        DOUBLE PRECISION MEAN,DISPER
!
        INTEGER NDIM,NEVAL
        REAL XX0(3),DXX0(3)      !valores iniciales y desplazamientos de prueba
        REAL XX(3),DXX(3)
        EXTERNAL FUNKGSS_CASEA
        REAL FUNKGSS_CASEA
!
        COMMON/BLKFITG1_CASEA/NP
        COMMON/BLKFITG2_CASEA/XF,YF
        COMMON/BLKX0INI_CASEA/X0INI
        COMMON/BLKDELTAX_CASEA/DELTAX_
!------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
!------------------------------------------------------------------------------
        DELTAX_=DELTAX !duplicate variable for common
!------------------------------------------------------------------------------
        AMP=0.
        EAMP=0.
        EEAMP=0.
        X0=0.
        EX0=0.
        EEX0=0.
        SIGMA=0.
        ESIGMA=0.
        EESIGMA=0.
!
        IF(NPFIT.LT.5)THEN
          WRITE(*,101)'ERROR: in subroutine GAUSSFIT:'
          WRITE(*,101)' No. of points for fit < 3'
          WRITE(*,100,ADVANCE='NO')'Press <CR> to continue...'
          RETURN
        ELSEIF(NPFIT.GT.NCMAX)THEN
          WRITE(*,101)'ERROR: in subroutine GAUSSFIT:'
          WRITE(*,101)' No. of points for fit > NCMAX'
          WRITE(*,100,ADVANCE='NO')'Press <CR> to continue...'
          RETURN
        END IF
!
        IF(NSIMUL.GT.NSIMULMAX)THEN
          WRITE(*,101)'ERROR: in subroutine GAUSSFIT:'
          WRITE(*,101)' No. of simulations > NSIMULMAX'
          WRITE(*,100,ADVANCE='NO')'Press <CR> to continue...'
          RETURN
        END IF
!
        NP=NPFIT
!------------------------------------------------------------------------------
        DO I=1,NP
          XF(I)=XFIT(I)
          YF(I)=YFIT(I)
        END DO
!------------------------------------------------------------------------------
! como primera estimacion del centro de la gaussiana calculamos el ajuste
! a una parabola y determinamos el maximo/minimo
        CALL POLFIT(XF,YF,YF,NP,3,0,A,CHISQR)
        X0INI=-REAL(A(2)/(2*A(3)))
        X0=0.0           !DOWNHILL supondra que el maximo esta alrededor de x=0
! como estimacion de la amplitud tomamos 3 veces el valor de la dispersion
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
! Como estimacion del valor de SIGMA tomamos un octavo del recorrido en la
! variable X
        XMIN=XF(1)
        XMAX=XMIN
        DO I=2,NP
          IF(XF(I).LT.XMIN) XMIN=XF(I)
          IF(XF(I).GT.XMAX) XMAX=XF(I)
        END DO
        SIGMA=(XMAX-XMIN)/8.
!------------------------------------------------------------------------------
! Con estas estimaciones ya podemos calcular los parametros de entrada de
! la subrutina DOWNHILL
        NDIM=3
        XX0(1)=AMP
        XX0(2)=X0-DELTAX/2.0
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
!
        CALL DOWNHILL(NDIM,XX0,DXX0,FUNKGSS_CASEA,1.0,0.5,2.0,YRMSTOL,&
     &   XX,DXX,NEVAL)
!
        AMP=XX(1)
        X0=XX(2)+X0INI
        SIGMA=XX(3)
        EEAMP=DXX(1)
        EEX0=DXX(2)
        EESIGMA=DXX(3)
!------------------------------------------------------------------------------
        IF(NSIMUL.LT.2)RETURN       !si no hay que calcular errores, regresamos
!------------------------------------------------------------------------------
! guardamos los valores finales como valores de prueba para el calculo de
! errores
        DO J=1,3
          XX0(J)=XX(J)
          IF(DXX(J).GT.0.0)DXX0(J)=DXX(J) !de lo contrario usa el valor inicial
        END DO
!------------------------------------------------------------------------------
        NSEED=-1
        DO ISIMUL=1,NSIMUL
          DO I=1,NP
            R1=RANRED(NSEED)
            R2=RANRED(NSEED)
            ERRYF=1.41421356*EYFIT(I)*SQRT(-1.*LOG(1.-R1))*COS(PI2*R2)
            YF(I)=YFIT(I)+ERRYF
          END DO
          CALL DOWNHILL(NDIM,XX0,DXX0,FUNKGSS_CASEA,1.0,0.5,2.0,YRMSTOL,&
     &     XX,DXX,NEVAL)
          AMP_SIMUL(ISIMUL)=XX(1)
          X0_SIMUL(ISIMUL)=XX(2)+X0INI
          SIGMA_SIMUL(ISIMUL)=XX(3)
        END DO
!
! EAMP: error en AMP
        MEAN=0.D0
        DO ISIMUL=1,NSIMUL
          MEAN=MEAN+DBLE(AMP_SIMUL(ISIMUL))
        END DO
        MEAN=MEAN/DBLE(NSIMUL)
        DISPER=0.D0
        DO ISIMUL=1,NSIMUL
          DISPER=DISPER+(DBLE(AMP_SIMUL(ISIMUL))-MEAN)*&
     &     (DBLE(AMP_SIMUL(ISIMUL))-MEAN)
        END DO
        DISPER=DSQRT(DISPER/DBLE(NSIMUL-1))
        EAMP=REAL(DISPER)
!
! EX0: error en X0
        MEAN=0.D0
        DO ISIMUL=1,NSIMUL
          MEAN=MEAN+DBLE(X0_SIMUL(ISIMUL))
        END DO
        MEAN=MEAN/DBLE(NSIMUL)
        DISPER=0.D0
        DO ISIMUL=1,NSIMUL
          DISPER=DISPER+(DBLE(X0_SIMUL(ISIMUL))-MEAN)*&
     &     (DBLE(X0_SIMUL(ISIMUL))-MEAN)
        END DO
        DISPER=DSQRT(DISPER/DBLE(NSIMUL-1))
        EX0=REAL(DISPER)
!
! ESIGMA: error en SIGMA
        MEAN=0.D0
        DO ISIMUL=1,NSIMUL
          MEAN=MEAN+DBLE(SIGMA_SIMUL(ISIMUL))
        END DO
        MEAN=MEAN/DBLE(NSIMUL)
        DISPER=0.D0
        DO ISIMUL=1,NSIMUL
          DISPER=DISPER+(DBLE(SIGMA_SIMUL(ISIMUL))-MEAN)*&
     &     (DBLE(SIGMA_SIMUL(ISIMUL))-MEAN)
        END DO
        DISPER=DSQRT(DISPER/DBLE(NSIMUL-1))
        ESIGMA=REAL(DISPER)
!------------------------------------------------------------------------------
100     FORMAT(A)
101     FORMAT(A)
        END
!
!******************************************************************************
!
! XX(1)=AMP, XX(2)=X0, XX(3)=SIGMA
        REAL FUNCTION FUNKGSS_CASEA(XX)
        USE REDLIB_M
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        REAL XX(3)
!
        INTEGER I
        INTEGER NF
        REAL XF(NCMAX),YF(NCMAX)
        REAL DELTAX
        COMMON/BLKFITG1_CASEA/NF
        COMMON/BLKFITG2_CASEA/XF,YF
        COMMON/BLKX0INI_CASEA/X0INI
        COMMON/BLKDELTAX_CASEA/DELTAX
        REAL FF,X0INI
!------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
! Introducimos X0INI para que el ajuste se realice alrededor de x=0. De esta
! forma el ajuste es mejor
        FUNKGSS_CASEA=0.
        DO I=1,NF
          FF=XX(1)*EXP(-(((XF(I)-XX(2)-X0INI)*(XF(I)-XX(2)-X0INI))/&
     &     (2.*XX(3)*XX(3))))
          FF=FF+XX(1)*EXP(-(((XF(I)-XX(2)-DELTAX-X0INI)*&
     &                       (XF(I)-XX(2)-DELTAX-X0INI))/&
     &     (2.*XX(3)*XX(3))))
          FUNKGSS_CASEA=FUNKGSS_CASEA+(YF(I)-FF)*(YF(I)-FF)
        END DO
        END
