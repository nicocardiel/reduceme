C------------------------------------------------------------------------------
C Version 5-May-1999                                            File: broaden.f
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
C SUBROUTINE BROADEN(S1,S2,NCHAN,STWV,DISP,SIGMA,LERR)
C
C Input: S1,NCHAN,STWV,DISP,SIGMA,LERR
C Output: S2
C
C Broadens a single spectrum by convolving with a gaussian (variable along the
C wavelength scale).
C
C REAL    S1(NCMAX) -> input spectrum
C REAL    S2(NCMAX) -> output spectrum
C INTEGER NCHAN -> no. of channels
C REAL    STWV -> central wavelength of the first channel
C REAL    DISP -> dispersion (Angstrom/pixel)
C REAL    SIGMA(NCMAX) -> sigma value of the gaussian to be applied in each
C                         channel
C LOGICAL LERR -> if .TRUE. S1 corresponds to an error spectrum
C
Comment
C------------------------------------------------------------------------------
        SUBROUTINE BROADEN(S1,S2,NCHAN,STWV,DISP,SIGMA,LERR)
C
        IMPLICIT NONE
        INTEGER NCHAN
        REAL S1(NCHAN),S2(NCHAN)
        REAL STWV,DISP
        REAL SIGMA(NCHAN)
        LOGICAL LERR
C
        REAL PI
        PARAMETER (PI=3.141592654)
        REAL C
        PARAMETER (C=299792.458)
C
        INTEGER I,K,K1,K2,N,INC
        REAL W0,W,W1,W2
        REAL SUM,FACTOR1,FACTOR2
        REAL PIXSIGMA
        REAL FINTGAUSS
        REAL FINTGAUSSE
        REAL FACTOR01,FACTOR02,PIXSIGMA0
C------------------------------------------------------------------------------
C calculamos todas las constantes fuera de los bucles para aumentar velocidad
C------------------------------------------------------------------------------
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
C
        END
