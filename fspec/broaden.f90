!------------------------------------------------------------------------------
! Version 20-November-2026                                    File: broaden.f90
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
! SUBROUTINE BROADEN(S1,S2,NCHAN,STWV,DISP,SIGMA,LERR)
!
! Input: S1,NCHAN,STWV,DISP,SIGMA,LERR
! Output: S2
!
! Broadens a single spectrum by convolving with a gaussian (variable along the
! wavelength scale).
!
! REAL    S1(NCMAX) -> input spectrum
! REAL    S2(NCMAX) -> output spectrum
! INTEGER NCHAN -> no. of channels
! REAL    STWV -> central wavelength of the first channel
! REAL    DISP -> dispersion (Angstrom/pixel)
! REAL    SIGMA(NCMAX) -> sigma value of the gaussian to be applied in each
!                         channel
! LOGICAL LERR -> if .TRUE. S1 corresponds to an error spectrum
!
!Comment
!------------------------------------------------------------------------------
        SUBROUTINE BROADEN(S1,S2,NCHAN,STWV,DISP,SIGMA,LERR)
!
        IMPLICIT NONE
        INTEGER NCHAN
        REAL S1(NCHAN),S2(NCHAN)
        REAL STWV,DISP
        REAL SIGMA(NCHAN)
        LOGICAL LERR
!
        REAL PI
        PARAMETER (PI=3.141592654)
        REAL C
        PARAMETER (C=299792.458)
!
        INTEGER I,K,K1,K2,N,INC
        REAL W0,W,W1,W2
        REAL SUM,FACTOR1,FACTOR2
        REAL PIXSIGMA
        REAL FACTOR01,FACTOR02,PIXSIGMA0
!------------------------------------------------------------------------------
! calculamos todas las constantes fuera de los bucles para aumentar velocidad
!------------------------------------------------------------------------------
        IF(LERR)THEN
          DO I=1,NCHAN
            IF(SIGMA(I).GT.0.0)THEN
              FACTOR01=-C*C/(2.*SIGMA(I)*SIGMA(I))
              FACTOR02=C/(SQRT(2.*PI)*SIGMA(I))
              PIXSIGMA0=SIGMA(I)/(C*DISP)
              W0=REAL(I-1)*DISP+STWV                          !longitud de onda
              SUM=0.0
              FACTOR1=FACTOR01/(W0*W0)
              FACTOR2=FACTOR02/W0
              PIXSIGMA=PIXSIGMA0*W0     !numero de pixels que equivalen a sigma
              INC=NINT(6.*PIXSIGMA)               !numero de pixels a cada lado
              K1=I-INC                               !limite inferior: -6 sigma
              IF(K1.LT.1) K1=1                       !ojo con el borde inferior
              K2=I+INC                               !limite superior: +6 sigma
              IF(K2.GT.NCHAN) K2=NCHAN               !ojo con el borde superior
              N=NINT(10./PIXSIGMA) !exigimos 10 intervalos para muestrear sigma
              IF(MOD(N,2).NE.0) N=N+1                 !N debe ser un numero par
              IF(N.LT.10) N=10          !como minimo usamos 10 intervalos/pixel
              DO K=K1,K2  !sumamos entre +-6 sigma (salvo en los bordes, claro)
                W=REAL(K-1)*DISP+STWV      !l.d.o. del centro del pixel K-esimo
                W1=W-DISP/2.                 !l.d.o. inferior del pixel K-esimo
                W2=W+DISP/2.                 !l.d.o. superior del pixel K-esimo
                SUM=SUM+S1(K)*S1(K)*FINTGAUSSE(W1,W2,N,W0,FACTOR1)
              END DO
              SUM=SQRT(SUM)
              S2(I)=FACTOR2*SUM
            ELSE
              S2(I)=S1(I)
            END IF
          END DO
        ELSE
          DO I=1,NCHAN
            IF(SIGMA(I).GT.0.0)THEN
              FACTOR01=-C*C/(2.*SIGMA(I)*SIGMA(I))
              FACTOR02=C/(SQRT(2.*PI)*SIGMA(I))
              PIXSIGMA0=SIGMA(I)/(C*DISP)
              W0=REAL(I-1)*DISP+STWV                          !longitud de onda
              SUM=0.0
              FACTOR1=FACTOR01/(W0*W0)
              FACTOR2=FACTOR02/W0
              PIXSIGMA=PIXSIGMA0*W0     !numero de pixels que equivalen a sigma
              INC=NINT(6.*PIXSIGMA)               !numero de pixels a cada lado
              K1=I-INC                               !limite inferior: -6 sigma
              IF(K1.LT.1) K1=1                       !ojo con el borde inferior
              K2=I+INC                               !limite superior: +6 sigma
              IF(K2.GT.NCHAN) K2=NCHAN               !ojo con el borde superior
              N=NINT(10./PIXSIGMA) !exigimos 10 intervalos para muestrear sigma
              IF(MOD(N,2).NE.0) N=N+1                 !N debe ser un numero par
              IF(N.LT.10) N=10          !como minimo usamos 10 intervalos/pixel
              DO K=K1,K2  !sumamos entre +-6 sigma (salvo en los bordes, claro)
                W=REAL(K-1)*DISP+STWV      !l.d.o. del centro del pixel K-esimo
                W1=W-DISP/2.                 !l.d.o. inferior del pixel K-esimo
                W2=W+DISP/2.                 !l.d.o. superior del pixel K-esimo
                SUM=SUM+S1(K)*FINTGAUSS(W1,W2,N,W0,FACTOR1)
              END DO
              S2(I)=FACTOR2*SUM
            ELSE
              S2(I)=S1(I)
            END IF
          END DO
        END IF
!
        END
