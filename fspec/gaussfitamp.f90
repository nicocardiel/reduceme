!------------------------------------------------------------------------------
! Version 20-November-2026                                File: gaussfitamp.f90
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
! SUBROUTINE GAUSSFITAMP(NPFIT,XFIT,YFIT,EYFIT,X0,SIGMA,EX0,ESIGMA,MINSIGMA,
!                        AMP,EAMP,EEAMP,NSIMUL)
!
! Input: NPFIT,XFIT,YFIT,EYFIT,X0,SIGMA,EX0,ESIGMA,MINSIGMA,NSIMUL
! Output: AMP,EAMP,EEAMP
!
! Fit of AMP of a gaussian with X0 and SIGMA fixed:
! Y=AMP*EXP[-((X-X0)^2/(2*SIGMA^2))]
!
! INTEGER NPFIT -> number of points to be fitted
! REAL XFIT,YFIT,EYFIT -> x, y and error
! REAL X0, EX0 -> center of the gaussian and its error
! REAL SIGMA, ESIGMA -> sigma value of the gaussian and its error
! REAL MINSIGMA -> minimum SIGMA allowed in simulations (typically MINSIGMA
!                  must be the spectral resolution)
! REAL AMP -> maximum of the fitted gaussian
! REAL EAMP -> error in SIGMA (due to EYFIT)
! REAL EEAMP -> error in SIGMA (due to EX0 and ESIGMA ---simulations---)
! INTEGER NSIMUL -> number of simulations to compute EEAMP
!
!Comment
!------------------------------------------------------------------------------
        SUBROUTINE GAUSSFITAMP(NPFIT,XFIT,YFIT,EYFIT,X0,SIGMA,&
     &   EX0,ESIGMA,MINSIGMA,AMP,EAMP,EEAMP,NSIMUL)
        IMPLICIT NONE
        INTEGER NPFIT
        REAL XFIT(NPFIT),YFIT(NPFIT),EYFIT(NPFIT)
        REAL X0,SIGMA,EX0,ESIGMA,MINSIGMA
        REAL AMP,EAMP,EEAMP
        INTEGER NSIMUL
!
        INTEGER NSIMULMAX
        PARAMETER(NSIMULMAX=1000)                !numero maximo de simulaciones
        REAL PI2
        PARAMETER(PI2=6.283185307)               !2 x pi
!
        INTEGER I,J,ISIMUL
        INTEGER NEXIT
        INTEGER NSEED
        REAL X0_,SIGMA_
        REAL ERR_X0,ERR_SIGMA
        REAL AMP_SIMUL(NSIMULMAX)
        REAL R1,R2
        DOUBLE PRECISION DSUM1,DSUM2,DERIVATIVE
        DOUBLE PRECISION MEAN,DISPER
        LOGICAL LOOP
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
        AMP=0.
        EAMP=0.
        EEAMP=0.
!
        IF(NPFIT.LT.3)THEN
          WRITE(*,101)'ERROR: in subroutine GAUSSFITAMP:'
          WRITE(*,101)' No. of points for fit < 3'
          WRITE(*,100,ADVANCE='NO')'Press <CR> to continue...'
          RETURN
        END IF
!
        IF(NSIMUL.GT.NSIMULMAX)THEN
          WRITE(*,101)'ERROR: in subroutine GAUSSFITAMP:'
          WRITE(*,101)' No. of simulations > NSIMULMAX'
          WRITE(*,100,ADVANCE='NO')'Press <CR> to continue...'
          RETURN
        END IF
!------------------------------------------------------------------------------
        DSUM1=0.D0
        DSUM2=0.D0
        DO I=1,NPFIT
          DSUM1=DSUM1+DBLE(&
     &     YFIT(I)*EXP(-(XFIT(I)-X0)*(XFIT(I)-X0)/(2.*SIGMA*SIGMA))/&
     &     (EYFIT(I)*EYFIT(I))&
     &                    )
          DSUM2=DSUM2+DBLE(&
     &     EXP(-(XFIT(I)-X0)*(XFIT(I)-X0)/(SIGMA*SIGMA))/&
     &     (EYFIT(I)*EYFIT(I))&
     &                    )
        END DO
        AMP=REAL(DSUM1/DSUM2)
        DSUM1=0.D0
        DO J=1,NPFIT
          DERIVATIVE=DBLE(&
     &     EXP(-(XFIT(J)-X0)*(XFIT(J)-X0)/(2.*SIGMA*SIGMA))/&
     &     (EYFIT(J)*EYFIT(J))&
     &                    )/DSUM2
           DSUM1=DSUM1+DERIVATIVE*DERIVATIVE*EYFIT(J)*EYFIT(J)
        END DO
        EAMP=REAL(DSQRT(DSUM1))
!------------------------------------------------------------------------------
        IF(NSIMUL.LT.2)RETURN     !si no hay que hacer simulaciones, regresamos
!------------------------------------------------------------------------------

        NSEED=-1
        DO ISIMUL=1,NSIMUL
          R1=RANRED(NSEED)
          R2=RANRED(NSEED)
          ERR_X0=1.41421356*EX0*SQRT(-1.*LOG(1.-R1))*COS(PI2*R2)
          X0_=X0+ERR_X0
          LOOP=.TRUE.
          NEXIT=0
          DO WHILE(LOOP)
            R1=RANRED(NSEED)
            R2=RANRED(NSEED)
            ERR_SIGMA=1.41421356*ESIGMA*SQRT(-1.*LOG(1.-R1))*COS(PI2*R2)
            SIGMA_=SIGMA+ERR_SIGMA
            LOOP=(SIGMA_.LT.MINSIGMA)
            IF(LOOP)THEN
              NEXIT=NEXIT+1
              IF(NEXIT.GT.100*NSIMUL)THEN
                WRITE(*,101) 'FATAL ERROR: in subroutine GAUSSFITAMP'
                WRITE(*,100,ADVANCE='NO') 'SIGMA, ERRSIGMA, MINSIGMA: '
                WRITE(*,*) SIGMA,ERR_SIGMA,MINSIGMA
                WRITE(*,101) 'Check previous values!'
                STOP
              END IF
            END IF
          END DO
          DSUM1=0.D0
          DSUM2=0.D0
          DO I=1,NPFIT
            DSUM1=DSUM1+DBLE(&
     &       YFIT(I)*EXP(-(XFIT(I)-X0_)*(XFIT(I)-X0_)/&
     &       (2.*SIGMA_*SIGMA_))/(EYFIT(I)*EYFIT(I))&
     &                      )
            DSUM2=DSUM2+DBLE(&
     &       EXP(-(XFIT(I)-X0_)*(XFIT(I)-X0_)/(SIGMA_*SIGMA_))/&
     &       (EYFIT(I)*EYFIT(I))&
     &                      )
          END DO
          AMP_SIMUL(ISIMUL)=REAL(DSUM1/DSUM2)
        END DO
!
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
        EEAMP=REAL(DISPER)
!------------------------------------------------------------------------------
100     FORMAT(A)
101     FORMAT(A)
        END
