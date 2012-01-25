C------------------------------------------------------------------------------
C Version 07-September-2007                                  File: gauss2bfit.f
C------------------------------------------------------------------------------
Comment
C
C SUBROUTINE GAUSS2BFIT(NPFIT,XFIT,YFIT,EYFIT,DELTAX,X0,SIGMA,AMP1,AMP2,
C                       EX0,ESIGMA,EAMP1,EAMP2,EEX0,EESIGMA,EEAMP1,EEAMP2,
C                       YRMSTOL,NSIMUL)
C
C Input: NPFIT,XFIT,YFIT,EYFIT,DELTAX,YRMSTOL,NSIMUL
C Output: X0,SIGMA,AMP1,AMP2,EX0,ESIGMA,EAMP1,EAMP2,EEX0,EESIGMA,EEAMP1,EEAMP2
C
C Numerical fit of 2 gaussians with the same width and different area 
C (using DOWNHILL), with a fixed separation given by DELTAX.
C Y=AMP1*EXP[-((X-X0)^2/(2*SIGMA^2))]+AMP2*EXP[-((X-X0-DELTAX)^2/(2*SIGMA^2))]
C
C INTEGER NPFIT -> number of points to be fitted
C REAL XFIT,YFIT,EYFIT -> x, y and error
C REAL DELTAX -> separation between gaussians
C REAL X0 -> center of the fitted gaussian
C REAL SIGMA -> sigma value of the fitted gaussian
C REAL AMP1 -> maximum of the fitted gaussian #1
C REAL AMP2 -> maximum of the fitted gaussian #2
C REAL EX0,ESIGMA,EAMPn-> errors in X0,SIGMA,AMPn (due to EYFIT --simulations--)
C REAL EEX0,EESIGMA,EEAMPn-> errors in X0,SIGMA,AMPn (rms from DOWNHILL)
C REAL YRMSTOL -> stopping criterion for DOWNHILL
C INTEGER NSIMUL -> number of simulations to compute errors
C
Comment
C------------------------------------------------------------------------------
        SUBROUTINE GAUSS2BFIT(NPFIT,XFIT,YFIT,EYFIT,DELTAX,
     +   X0,SIGMA,AMP1,AMP2,
     +   EX0,ESIGMA,EAMP1,EAMP2,
     +   EEX0,EESIGMA,EEAMP1,EEAMP2,YRMSTOL,NSIMUL) 
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INTEGER NPFIT
        REAL XFIT(NPFIT),YFIT(NPFIT),EYFIT(NPFIT)
        REAL DELTAX
        REAL X0,SIGMA,AMP1,AMP2
        REAL EX0,ESIGMA,EAMP1,EAMP2
        REAL EEX0,EESIGMA,EEAMP1,EEAMP2
        REAL YRMSTOL
        INTEGER NSIMUL
C
        INTEGER NSIMULMAX
        PARAMETER(NSIMULMAX=1000)                !numero maximo de simulaciones
        REAL PI2
        PARAMETER(PI2=6.283185307)               !2 x pi
C
        INTEGER I,J,ISIMUL
        INTEGER NP
        INTEGER NSEED
        REAL X0INI
        REAL XMIN,XMAX
        REAL XF(NCMAX),YF(NCMAX),ERRYF
        REAL A(3),CHISQR
        REAL X0_SIMUL(NSIMULMAX)
        REAL SIGMA_SIMUL(NSIMULMAX)
        REAL AMP1_SIMUL(NSIMULMAX)
        REAL AMP2_SIMUL(NSIMULMAX)
        REAL DELTAX_
        REAL RANRED,R1,R2
        DOUBLE PRECISION MEAN,DISPER
C
        INTEGER NDIM,NEVAL
        REAL XX0(4),DXX0(4)      !valores iniciales y desplazamientos de prueba
        REAL XX(4),DXX(4)
        EXTERNAL FUNKGSS_CASEB
        REAL FUNKGSS_CASEB
C
        COMMON/BLKFITG1_CASEB/NP
        COMMON/BLKFITG2_CASEB/XF,YF
        COMMON/BLKX0INI_CASEB/X0INI
        COMMON/BLKDELTAX_CASEB/DELTAX_
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C------------------------------------------------------------------------------
        DELTAX_=DELTAX !duplicate variable for common
C------------------------------------------------------------------------------
        AMP1=0.
        AMP2=0.
        EAMP1=0.
        EAMP2=0.
        EEAMP1=0.
        EEAMP2=0.
        X0=0.
        EX0=0.
        EEX0=0.
        SIGMA=0.
        ESIGMA=0.
        EESIGMA=0.
C
        IF(NPFIT.LT.5)THEN
          WRITE(*,101)'ERROR: in subroutine GAUSS2BFIT:'
          WRITE(*,101)' No. of points for fit < 5'
          WRITE(*,100)'Press <CR> to continue...'
          RETURN
        ELSEIF(NPFIT.GT.NCMAX)THEN
          WRITE(*,101)'ERROR: in subroutine GAUSS2BFIT:'
          WRITE(*,101)' No. of points for fit > NCMAX'
          WRITE(*,100)'Press <CR> to continue...'
          RETURN
        END IF
C
        IF(NSIMUL.GT.NSIMULMAX)THEN
          WRITE(*,101)'ERROR: in subroutine GAUSS2BFIT:'
          WRITE(*,101)' No. of simulations > NSIMULMAX'
          WRITE(*,100)'Press <CR> to continue...'
          RETURN
        END IF
C
        NP=NPFIT
C------------------------------------------------------------------------------
        DO I=1,NP
          XF(I)=XFIT(I)
          YF(I)=YFIT(I)
        END DO
C------------------------------------------------------------------------------
C como primera estimacion del centro de la gaussiana calculamos el ajuste
C a una parabola y determinamos el maximo/minimo
        CALL POLFIT(XF,YF,YF,NP,3,0,A,CHISQR)
        X0INI=-REAL(A(2)/(2*A(3)))
        X0=0.0           !DOWNHILL supondra que el maximo esta alrededor de x=0
C como estimacion de la amplitud tomamos 3 veces el valor de la dispersion
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
        AMP1=3.*REAL(DISPER)
        AMP2=AMP1 !hacemos las dos amplitudes iguales inicialmente
C si la parabola tiene un minimo (derivada segunda positiva), la 
C amplitud es negativa
        IF(A(3).GT.0.D0)THEN
          AMP1=-AMP1
          AMP2=-AMP2
        END IF
C
C Como estimacion del valor de SIGMA tomamos un octavo del recorrido en la
C variable X
        XMIN=XF(1)
        XMAX=XMIN
        DO I=2,NP
          IF(XF(I).LT.XMIN) XMIN=XF(I)
          IF(XF(I).GT.XMAX) XMAX=XF(I)
        END DO
        SIGMA=(XMAX-XMIN)/8.
C------------------------------------------------------------------------------
C Con estas estimaciones ya podemos calcular los parametros de entrada de
C la subrutina DOWNHILL
        NDIM=4
        XX0(1)=AMP1
        XX0(2)=AMP2
        XX0(3)=X0-DELTAX/2.0
        XX0(4)=SIGMA
        IF(AMP1.EQ.0.)THEN
          DXX0(1)=0.1
        ELSE
          DXX0(1)=0.1*AMP1
        END IF
        IF(AMP2.EQ.0.)THEN
          DXX0(2)=0.1
        ELSE
          DXX0(2)=0.1*AMP2
        END IF
        DXX0(3)=1.
        IF(SIGMA.EQ.0.)THEN
          DXX0(4)=0.1
        ELSE
          DXX0(4)=0.1*SIGMA
        END IF
C
        CALL DOWNHILL(NDIM,XX0,DXX0,FUNKGSS_CASEB,1.0,0.5,2.0,YRMSTOL,
     +   XX,DXX,NEVAL)
C
        AMP1=XX(1)
        AMP2=XX(2)
        X0=XX(3)+X0INI
        SIGMA=XX(4)
        EEAMP1=DXX(1)
        EEAMP2=DXX(2)
        EEX0=DXX(3)
        EESIGMA=DXX(4)
C------------------------------------------------------------------------------
        IF(NSIMUL.LT.2)RETURN       !si no hay que calcular errores, regresamos
C------------------------------------------------------------------------------
C guardamos los valores finales como valores de prueba para el calculo de
C errores
        DO J=1,4
          XX0(J)=XX(J)
          IF(DXX(J).GT.0.0)DXX0(J)=DXX(J) !de lo contrario usa el valor inicial
        END DO
C------------------------------------------------------------------------------
        NSEED=-1
        DO ISIMUL=1,NSIMUL
          DO I=1,NP
            R1=RANRED(NSEED)
            R2=RANRED(NSEED)
            ERRYF=1.41421356*EYFIT(I)*SQRT(-1.*LOG(1.-R1))*COS(PI2*R2)
            YF(I)=YFIT(I)+ERRYF
          END DO
          CALL DOWNHILL(NDIM,XX0,DXX0,FUNKGSS_CASEB,1.0,0.5,2.0,YRMSTOL,
     +     XX,DXX,NEVAL)
          AMP1_SIMUL(ISIMUL)=XX(1)
          AMP2_SIMUL(ISIMUL)=XX(2)
          X0_SIMUL(ISIMUL)=XX(3)+X0INI
          SIGMA_SIMUL(ISIMUL)=XX(4)
        END DO
C
C EAMP: error en AMPn
        MEAN=0.D0
        DO ISIMUL=1,NSIMUL
          MEAN=MEAN+DBLE(AMP1_SIMUL(ISIMUL))
        END DO
        MEAN=MEAN/DBLE(NSIMUL)
        DISPER=0.D0
        DO ISIMUL=1,NSIMUL
          DISPER=DISPER+(DBLE(AMP1_SIMUL(ISIMUL))-MEAN)*
     +     (DBLE(AMP1_SIMUL(ISIMUL))-MEAN)
        END DO
        DISPER=DSQRT(DISPER/DBLE(NSIMUL-1))
        EAMP1=REAL(DISPER)
C
        MEAN=0.D0
        DO ISIMUL=1,NSIMUL
          MEAN=MEAN+DBLE(AMP2_SIMUL(ISIMUL))
        END DO
        MEAN=MEAN/DBLE(NSIMUL)
        DISPER=0.D0
        DO ISIMUL=1,NSIMUL
          DISPER=DISPER+(DBLE(AMP2_SIMUL(ISIMUL))-MEAN)*
     +     (DBLE(AMP2_SIMUL(ISIMUL))-MEAN)
        END DO
        DISPER=DSQRT(DISPER/DBLE(NSIMUL-1))
        EAMP2=REAL(DISPER)
C
C EX0: error en X0
        MEAN=0.D0
        DO ISIMUL=1,NSIMUL
          MEAN=MEAN+DBLE(X0_SIMUL(ISIMUL))
        END DO
        MEAN=MEAN/DBLE(NSIMUL)
        DISPER=0.D0
        DO ISIMUL=1,NSIMUL
          DISPER=DISPER+(DBLE(X0_SIMUL(ISIMUL))-MEAN)*
     +     (DBLE(X0_SIMUL(ISIMUL))-MEAN)
        END DO
        DISPER=DSQRT(DISPER/DBLE(NSIMUL-1))
        EX0=REAL(DISPER)
C
C ESIGMA: error en SIGMA
        MEAN=0.D0
        DO ISIMUL=1,NSIMUL
          MEAN=MEAN+DBLE(SIGMA_SIMUL(ISIMUL))
        END DO
        MEAN=MEAN/DBLE(NSIMUL)
        DISPER=0.D0
        DO ISIMUL=1,NSIMUL
          DISPER=DISPER+(DBLE(SIGMA_SIMUL(ISIMUL))-MEAN)*
     +     (DBLE(SIGMA_SIMUL(ISIMUL))-MEAN)
        END DO
        DISPER=DSQRT(DISPER/DBLE(NSIMUL-1))
        ESIGMA=REAL(DISPER)
C------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
        END
C
C******************************************************************************
C
C XX(1)=AMP1, XX(2)=AMP2, XX(3)=X0, XX(4)=SIGMA
        REAL FUNCTION FUNKGSS_CASEB(XX)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        REAL XX(4)
C
        INTEGER I
        INTEGER NF
        REAL XF(NCMAX),YF(NCMAX)
        REAL DELTAX
        COMMON/BLKFITG1_CASEB/NF
        COMMON/BLKFITG2_CASEB/XF,YF
        COMMON/BLKX0INI_CASEB/X0INI
        COMMON/BLKDELTAX_CASEB/DELTAX
        REAL FF,X0INI
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C Introducimos X0INI para que el ajuste se realice alrededor de x=0. De esta
C forma el ajuste es mejor
        FUNKGSS_CASEB=0.
        DO I=1,NF
          FF=XX(1)*EXP(-(((XF(I)-XX(3)-X0INI)*(XF(I)-XX(3)-X0INI))/
     +     (2.*XX(4)*XX(4))))
          FF=FF+XX(2)*EXP(-(((XF(I)-XX(3)-DELTAX-X0INI)*
     +                       (XF(I)-XX(3)-DELTAX-X0INI))/
     +     (2.*XX(4)*XX(4))))
          FUNKGSS_CASEB=FUNKGSS_CASEB+(YF(I)-FF)*(YF(I)-FF)
        END DO
        END
