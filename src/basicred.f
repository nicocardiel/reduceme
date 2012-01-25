C------------------------------------------------------------------------------
C Version 28-November-1996                                     file: basicred.f
C------------------------------------------------------------------------------
C Copyright N. Cardiel & J. Gorgas, Departamento de Astrofisica
C Universidad Complutense de Madrid, 28040-Madrid, Spain
C E-mail: ncl@astrax.fis.ucm.es or fjg@astrax.fis.ucm.es
C------------------------------------------------------------------------------
C This program is free software; you can redistribute it and/or modify it
C under the terms of the GNU General Public License as published by the Free
C Software Foundation; either version 2 of the License, or (at your option) any
C later version. See the file gnu-public-license.txt for details.
C------------------------------------------------------------------------------
Comment
C
C Program: basicred
C Classification: arithmetic & manipulations
C Description: Determine the BIAS level and output the useful region of a 
C frame after dividing it by the required flatfields. This program also
C generates the error frames (gain and readout noise of the employed detector
C must be known).
C
Comment
C
C Realiza estadistica en region de BIAS, calcula el valor medio, mediana,...
C lo sustrae a la imagen y recorta esta al taman~o util. Finalmente puede
C sustraerse una imagen con la estructura del BIAS, una imagen con la
C estructura del DARK y dividir
C por una imagen de flatfield de alta frecuencia (flat de cupula o de lampara)
C y por una imagen de flatfield de baja frecuencia (flat de cielo). 
C El programa permite realizar asimismo el correspondiente calculo con las
C imagenes de errores. Por ultimo, permite salvar las imagenes invirtiendo
C la direccion espectral (por si acaso los espectros estan al reves).
C
        PROGRAM BASICRED
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        REAL READF
C
        INTEGER I,J,JJ,L,K
        INTEGER NS1,NS2,NC1,NC2,NSF,NCF
        INTEGER NPIX
        REAL A(NCMAX,NSMAX),F(NCMAX),ERR(NCMAX,NSMAX),EF(NCMAX)
        REAL S(NCMAX)
        REAL PIXEL(NCMAX*NSMAX)
        REAL DARK,DARKERR,DARKERR2,STD2
        REAL RNOISE,RNOISE2,GAIN
        REAL AIRMASS2,TIMEXPOS2
        REAL EXPOSURE
        REAL FBIAS,BIASMIN,BIASMAX
        REAL FMEDIAN,FMEDIAN1,FMEANSIGMA,FMEAN2
        DOUBLE PRECISION MEAN,STD
        CHARACTER*1 CFLAT,CIMAD,CERR,CDEIMA,CSURE,CFFEIMA,CREV,CBIAS
        CHARACTER*50 CDUMMY
        CHARACTER*75 INFILE,DARKFILE,FLATFILE,OUTFILE,ERRFILE
C NOTA: si en el fichero 'redlib.inc' se ha cambiado las dimensiones de
C OBJECT, FITSFILE y COMMENT, tambien hay que cambiarlas aqui en las
C variables que siguen:
        CHARACTER*255 OBJECT2,FITSFILE2,COMMENT2
        LOGICAL AL(NCMAX,NSMAX),LDATANEG,LDATANEGBIS
C------------------------------------------------------------------------------
        THISPROGRAM='basicred'
        CALL WELCOME('28-November-1996')
        CALL SHOWHLP('explanation')
C
        WRITE(*,100)'Create error image (y/n) '
        CERR(1:1)=READC('n','yn')
C
        WRITE(*,100)'Input file name'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        OBJECT2=OBJECT                 !guardamos las variables de este fichero
        FITSFILE2=FITSFILE
        COMMENT2=COMMENT
        AIRMASS2=AIRMASS
        TIMEXPOS2=TIMEXPOS
        DO I=1,NSCAN
          READ(20) (A(J,I),J=1,NCHAN)
        END DO
        CLOSE(20)
C------------------------------------------------------------------------------
C la matriz AL(,) determina que pixels se utilan para calcular el valor de BIAS
        DO I=1,NSCAN
          DO J=1,NCHAN
            AL(J,I)=.FALSE.
          END DO
        END DO
9       WRITE(*,101)'Introduce region to calculate BIAS:'
        CALL SHOWHLP('bias region')
10      WRITE(*,100)'1st & last scan (0,0=EXIT) '
        CALL READ2I('0,0',NS1,NS2)
        IF((NS1.EQ.0).AND.(NS2.EQ.0)) GOTO 25
        IF((NS1.LT.1).OR.(NS2.GT.NSCAN).OR.(NS1.GT.NS2))THEN
          WRITE(*,101)'ERROR: invalid entry. Try again.'
          GOTO 10
        END IF
20      WRITE(*,100)'1st & last channel.............. '
        CALL READ2I('@',NC1,NC2)
        IF((NC1.LT.1).OR.(NC2.GT.NCHAN).OR.(NC1.GT.NC2))THEN
          WRITE(*,101)'ERROR: invalid entry. Try again.'
          GOTO 20
        END IF
        DO I=NS1,NS2
          DO J=NC1,NC2
            AL(J,I)=.TRUE.
          END DO
        END DO
        GOTO 9
C
25      NPIX=0
        MEAN=0.D0
        BIASMIN=1.E30
        BIASMAX=-1.E30
        DO I=1,NSCAN
          DO J=1,NCHAN
            IF(AL(J,I))THEN
              MEAN=MEAN+DBLE(A(J,I))
              NPIX=NPIX+1
              IF(A(J,I).LT.BIASMIN) BIASMIN=A(J,I)
              IF(A(J,I).GT.BIASMAX) BIASMAX=A(J,I)
            END IF
          END DO
        END DO
C
        IF(NPIX.EQ.0)THEN
          WRITE(*,101)'WARNING: No. of pixels to calculate BIAS = 0.'
          MEAN=0.D0
          STD=0.D0
          FMEANSIGMA=0.
          FMEDIAN=0.
        ELSE
          MEAN=MEAN/DBLE(NPIX)
          STD=0.D0
          IF(NPIX.GT.1)THEN
            DO I=1,NSCAN
              DO J=1,NCHAN
                IF(AL(J,I))THEN
                  STD=STD+(DBLE(A(J,I))-MEAN)*(DBLE(A(J,I))-MEAN)
                END IF
              END DO
            END DO
            STD=DSQRT(STD/DBLE(NPIX-1))
          END IF
        END IF
C
        IF(NPIX.GT.1)THEN
          K=0
          DO I=1,NSCAN
            DO J=1,NCHAN
              IF(AL(J,I))THEN
                K=K+1
                PIXEL(K)=A(J,I)
              END IF
            END DO
          END DO
          FMEANSIGMA=FMEAN2(NPIX,PIXEL,3.0)  !media eliminando puntos a 3 sigma
          FMEDIAN=FMEDIAN1(NPIX,PIXEL)                                 !mediana
        END IF
C
        WRITE(*,*)
        WRITE(*,110)'>>> No. of pixels to measure bias.: ',NPIX
        WRITE(*,100)'>>> Measured mean bias value......: '
        WRITE(*,*)REAL(MEAN)
        WRITE(*,100)'>>> Mean (removing points > 3 sig): '
        WRITE(*,*)FMEANSIGMA
        WRITE(*,100)'>>> Measured median...............: '
        WRITE(*,*)FMEDIAN
        WRITE(*,100)'>>> Standard deviation............: '
        WRITE(*,*)REAL(STD)
        WRITE(*,100)'>>> Minimum & times stand.dev.....: '
        IF(NPIX.GE.1)THEN
          WRITE(*,*) BIASMIN,(REAL(MEAN)-BIASMIN)/REAL(STD)
        ELSE
          WRITE(*,*)0.,0.
        END IF
        WRITE(*,100)'>>> Maximum & times stand.dev.....: '
        IF(NPIX.GE.1)THEN
          WRITE(*,*) BIASMAX,(BIASMAX-REAL(MEAN))/REAL(STD)
        ELSE
          WRITE(*,*)0.,0.
        END IF
        WRITE(*,100)'>>> Bias error....................: '
        IF(NPIX.GE.1)THEN
          STD=STD/DBLE(NPIX)
        ELSE
          WRITE(*,*)0.,0.
        END IF
        WRITE(*,*)REAL(STD)
        WRITE(*,*)
        WRITE(*,101)'* Select BIAS value to be employed:'
        WRITE(*,101)'1 - Mean'
        WRITE(*,101)'2 - Mean excluding pixels > 3 sigma'
        WRITE(*,101)'3 - Median'
        WRITE(*,100)'Option '
        CBIAS(1:1)=READC('3','123')
        IF(CBIAS.EQ.'1')THEN
          FBIAS=REAL(MEAN)
        ELSEIF(CBIAS.EQ.'2')THEN
          FBIAS=FMEANSIGMA
        ELSE
          FBIAS=FMEDIAN
        END IF
C------------------------------------------------------------------------------
        WRITE(*,*)
        WRITE(*,101)'* NOTE: if necessary, a dark current image'//
     +   ' can be can be also'
        WRITE(*,101)'  subtracted later in this program.'
        WRITE(*,100)'Constant dark current value (in ADU/sec) '
        DARK=READF('0.0')
        IF(CERR.EQ.'y')THEN
          WRITE(*,100)'Error in constant dark current (in ADU/sec) '
          DARKERR=READF('0.0')
        END IF
        WRITE(CDUMMY,*) TIMEXPOS
        WRITE(*,100)'Exposure time (sec) '
        EXPOSURE=READF(CDUMMY)
        DARK=DARK*EXPOSURE
        DARKERR=DARKERR*EXPOSURE
        WRITE(CDUMMY,*) DARK
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,101)'>>> Efective Dark Current.......(ADUs): '//
     +   CDUMMY(1:L)
        WRITE(CDUMMY,*) DARKERR
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,101)'>>> Efective Dark Current Error (ADUs): '//
     +   CDUMMY(1:L)
        WRITE(*,*)
C------------------------------------------------------------------------------
        DO I=1,NSCAN
          DO J=1,NCHAN
            A(J,I)=A(J,I)-FBIAS-DARK
          END DO
        END DO
C
        LDATANEG=.FALSE.
        IF(CERR.EQ.'y')THEN
          WRITE(*,*)
          WRITE(*,100)'Read-out Noise (counts=Diginal Units, ADU)'
          RNOISE=READF('@')
          RNOISE2=RNOISE*RNOISE
          WRITE(*,100)'Gain (electrons/ADU)                      '
          GAIN=READF('@')
C error fotonico, ruido de lectura y sustraccion de bias+dark
          STD2=REAL(STD*STD)
          DARKERR2=DARKERR*DARKERR
          DO I=1,NSCAN                     !sumamos los errores cuadraticamente
            DO J=1,NCHAN
              IF(A(J,I).GE.0.)THEN
                ERR(J,I)=A(J,I)/GAIN                !A(J,I) sustraido BIAS+DARK
              ELSE             !si la senhal es negativa, el error se hace cero
                ERR(J,I)=0.                                     !por hacer algo
                LDATANEG=.TRUE.
              END IF
              ERR(J,I)=ERR(J,I)+RNOISE2+STD2+DARKERR2
            END DO
          END DO
          IF(LDATANEG)THEN
            WRITE(*,101)'WARNING#1: negative data values have '//
     +       'been found.'
          END IF
        END IF
C------------------------------------------------------------------------------
        CALL SHOWHLP('useful region')
        WRITE(*,101)'Introduce frame region to be saved:'
30      WRITE(*,100)'1st & last scan...'
        CALL READ2I('@',NS1,NS2)
        IF((NS1.LT.1).OR.(NS2.GT.NSCAN).OR.(NS1.GT.NS2))THEN
          WRITE(*,101)'ERROR: invalid entry. Try again.'
          GOTO 30
        END IF
        NSF=NS2-NS1+1
40      WRITE(*,100)'1st & last channel'
        CALL READ2I('@',NC1,NC2)
        IF((NC1.LT.1).OR.(NC2.GT.NCHAN).OR.(NC1.GT.NC2))THEN
          WRITE(*,101)'ERROR: invalid entry. Try again.'
          GOTO 40
        END IF
        NCF=NC2-NC1+1
C..............................................................................
        IF(LDATANEG)THEN
          LDATANEGBIS=.FALSE.
          DO I=NS1,NS2
            DO J=NC1,NC2
              IF(A(J,I).LT.0.) LDATANEGBIS=.TRUE.
            END DO
          END DO
          IF(LDATANEGBIS)THEN
            WRITE(*,101)'WARNING#2: negative data values have '//
     +       'been found.'
          ELSE
            WRITE(*,101)'WARNING#2: negative data values have '//
     +       'NOT been found.'
          END IF
        END IF
C------------------------------------------------------------------------------
        WRITE(*,100)'Subtract image of BIAS (y/n) '
        CIMAD(1:1)=READC('n','yn')
        IF(CIMAD.EQ.'n') GOTO 44
        WRITE(*,100)'Bias image file name......'
        DARKFILE=INFILEX(22,'@',NSF,NCF,STWV,DISP,21,.FALSE.)!............match
        IF(CERR.EQ.'y')THEN
          WRITE(*,100)'Are you using a Bias image error file (y/n) '
          CDEIMA(1:1)=READC('n','yn')
          IF(CDEIMA.EQ.'y')THEN
            WRITE(*,100)'Are you completely sure (y/n) '
            CSURE(1:1)=READC('y','yn')
            IF(CSURE.EQ.'y')THEN
              WRITE(*,100)'Bias image error file name'
              CALL GUESSEF(DARKFILE,ERRFILE)
              ERRFILE=INFILEX(23,ERRFILE,NSF,NCF,STWV,DISP,21,.TRUE.)!....match
            END IF
          END IF
        END IF
C sustraemos imagen de BIAS
        DO I=NS1,NS2
          IF(CERR.EQ.'y')THEN              !error de sustraer la imagen de BIAS
            IF(CDEIMA.EQ.'y')THEN
              IF(CSURE.EQ.'y')THEN
                READ(23) (EF(J),J=1,NCF)
                DO J=NC1,NC2
                  JJ=J-NC1+1
                  ERR(J,I)=ERR(J,I)+EF(JJ)*EF(JJ)
                END DO
              END IF
            END IF
          END IF
          READ(22) (F(J),J=1,NCF)
          DO J=NC1,NC2
            JJ=J-NC1+1
            A(J,I)=A(J,I)-F(JJ)
          END DO
        END DO
        CLOSE(22)
        IF(CERR.EQ.'y')THEN
          IF(CDEIMA.EQ.'y')THEN
            IF(CSURE.EQ.'y')THEN
              CLOSE(23)
            END IF
          END IF
        END IF
C..............................................................................
        IF(LDATANEG)THEN
          LDATANEGBIS=.FALSE.
          DO I=NS1,NS2
            DO J=NC1,NC2
              IF(A(J,I).LT.0.) LDATANEGBIS=.TRUE.
            END DO
          END DO
          IF(LDATANEGBIS)THEN
            WRITE(*,101)'WARNING#3: negative data values have '//
     +       'been found.'
          ELSE
            WRITE(*,101)'WARNING#3: negative data values have '//
     +       'NOT been found.'
          END IF
        END IF
C------------------------------------------------------------------------------
44      WRITE(*,100)'Subtract image of DARK (y/n) '
        CIMAD(1:1)=READC('n','yn')
        IF(CIMAD.EQ.'n') GOTO 45
        WRITE(*,100)'Dark image file name (ADUs/sec)'
        DARKFILE=INFILEX(22,'@',NSF,NCF,STWV,DISP,21,.FALSE.)!............match
        IF(CERR.EQ.'y')THEN
          WRITE(*,100)'Are you using a Dark image error file (y/n) '
          CDEIMA(1:1)=READC('n','yn')
          IF(CDEIMA.EQ.'y')THEN
            WRITE(*,100)'Are you completely sure (y/n) '
            CSURE(1:1)=READC('y','yn')
            IF(CSURE.EQ.'y')THEN
              WRITE(*,100)'Dark image error file name'
              CALL GUESSEF(DARKFILE,ERRFILE)
              ERRFILE=INFILEX(23,ERRFILE,NSF,NCF,STWV,DISP,21,.TRUE.)!....match
            END IF
          END IF
        END IF
C sustraemos imagen de DARK
        DO I=NS1,NS2
          IF(CERR.EQ.'y')THEN              !error de sustraer la imagen de BIAS
            IF(CDEIMA.EQ.'y')THEN
              IF(CSURE.EQ.'y')THEN
                READ(23) (EF(J),J=1,NCF)
                DO J=NC1,NC2
                  JJ=J-NC1+1
                  ERR(J,I)=ERR(J,I)+EF(JJ)*EF(JJ)*EXPOSURE*EXPOSURE
                END DO
              END IF
            END IF
          END IF
          READ(22) (F(J),J=1,NCF)
          DO J=NC1,NC2
            JJ=J-NC1+1
            A(J,I)=A(J,I)-F(JJ)*EXPOSURE
          END DO
        END DO
        CLOSE(22)
        IF(CERR.EQ.'y')THEN
          IF(CDEIMA.EQ.'y')THEN
            IF(CSURE.EQ.'y')THEN
              CLOSE(23)
            END IF
          END IF
        END IF
C..............................................................................
        IF(LDATANEG)THEN
          LDATANEGBIS=.FALSE.
          DO I=NS1,NS2
            DO J=NC1,NC2
              IF(A(J,I).LT.0.) LDATANEGBIS=.TRUE.
            END DO
          END DO
          IF(LDATANEGBIS)THEN
            WRITE(*,101)'WARNING#4: negative data values have '//
     +       'been found.'
          ELSE
            WRITE(*,101)'WARNING#4: negative data values have '//
     +       'NOT been found.'
          END IF
        END IF
C------------------------------------------------------------------------------
45      CALL SHOWHLP('high-frequency flatfield')
        WRITE(*,100)'Divide by a high-frequency flat field image (y/n) '
        CFLAT(1:1)=READC('y','yn')
        IF(CFLAT.EQ.'n') GOTO 46
        WRITE(*,100)'High-frequency flat field file name......'
        FLATFILE=INFILEX(25,'@',NSF,NCF,STWV,DISP,21,.FALSE.)!............match
        IF(CERR.EQ.'y')THEN
          WRITE(*,100)'Are you using a flat field error file name '
          CFFEIMA(1:1)=READC('y','yn')
          IF(CFFEIMA.EQ.'y')THEN
            WRITE(*,100)'Flat field error file name'
            CALL GUESSEF(FLATFILE,ERRFILE)
            ERRFILE=INFILEX(26,ERRFILE,NSF,NCF,STWV,DISP,21,.TRUE.)!......match
          ELSE
            DO J=1,NCF  !hay que usarlo, aunque sea cero, para escalar el error
              EF(J)=0.
            END DO
          END IF
        END IF
C dividimos por imagen de FLATFIELD
        DO I=NS1,NS2
          READ(25) (F(J),J=1,NCF)
          IF(CERR.EQ.'y')THEN                !error de dividir por el flatfield
            IF(CFFEIMA.EQ.'y')THEN
              READ(26) (EF(J),J=1,NCF)
            END IF
            DO J=NC1,NC2
              JJ=J-NC1+1
              ERR(J,I)=A(J,I)*A(J,I)*EF(JJ)*EF(JJ)+F(JJ)*F(JJ)*ERR(J,I)
              ERR(J,I)=ERR(J,I)/(F(JJ)*F(JJ)*F(JJ)*F(JJ))
            END DO
          END IF
          DO J=NC1,NC2
            JJ=J-NC1+1
            A(J,I)=A(J,I)/F(JJ)
          END DO
        END DO
        CLOSE(25)
        IF(CERR.EQ.'y')THEN
          IF(CFFEIMA.EQ.'y')THEN
            CLOSE(26)
          END IF
        END IF
C------------------------------------------------------------------------------
46      CALL SHOWHLP('low-frequency flatfield')
        WRITE(*,100)'Divide by a low-frequency flat field image (y/n) '
        CFLAT(1:1)=READC('y','yn')
        IF(CFLAT.EQ.'n') GOTO 50
        WRITE(*,100)'Low-frequency flat field file name.......'
        FLATFILE=INFILEX(25,'@',NSF,NCF,STWV,DISP,21,.FALSE.)!............match
        IF(CERR.EQ.'y')THEN
          WRITE(*,100)'Are you using a flat field error file name '
          CFFEIMA(1:1)=READC('y','yn')
          IF(CFFEIMA.EQ.'y')THEN
            WRITE(*,100)'Flat field error file name'
            CALL GUESSEF(FLATFILE,ERRFILE)
            ERRFILE=INFILEX(26,ERRFILE,NSF,NCF,STWV,DISP,21,.TRUE.)!......match
          ELSE
            DO J=1,NCF  !hay que usarlo, aunque sea cero, para escalar el error
              EF(J)=0.
            END DO
          END IF
        END IF
C dividimos por imagen de FLATFIELD
        DO I=NS1,NS2
          READ(25) (F(J),J=1,NCF)
          IF(CERR.EQ.'y')THEN                !error de dividir por el flatfield
            IF(CFFEIMA.EQ.'y')THEN
              READ(26) (EF(J),J=1,NCF)
            END IF
            DO J=NC1,NC2
              JJ=J-NC1+1
              ERR(J,I)=A(J,I)*A(J,I)*EF(JJ)*EF(JJ)+F(JJ)*F(JJ)*ERR(J,I)
              ERR(J,I)=ERR(J,I)/(F(JJ)*F(JJ)*F(JJ)*F(JJ))
            END DO
          END IF
          DO J=NC1,NC2
            JJ=J-NC1+1
            A(J,I)=A(J,I)/F(JJ)
          END DO
        END DO
        CLOSE(25)
        IF(CERR.EQ.'y')THEN
          IF(CFFEIMA.EQ.'y')THEN
            CLOSE(26)
          END IF
        END IF
C------------------------------------------------------------------------------
50      WRITE(*,*)
        WRITE(*,101)'* Bias, Dark & Flatfield files SHOULD NOT '//
     +   'be reversed'
        CALL SHOWHLP('reverse spectra')
        WRITE(*,100)'Reverse spectra in the wavelength direction (y/n) '
        CREV(1:1)=READC('n','yn')
C
        WRITE(*,100)'Output file name......'
        OBJECT=OBJECT2            !guardamos las variables del fichero original
        FITSFILE=FITSFILE2
        COMMENT=COMMENT2
        AIRMASS=AIRMASS2
        TIMEXPOS=TIMEXPOS2
        OUTFILE=OUTFILEX(30,'@',NSF,NCF,STWV,DISP,1,.FALSE.)
        IF(CREV.EQ.'n')THEN
          DO I=NS1,NS2
            DO J=NC1,NC2
              S(J)=A(J,I)
            END DO
            WRITE(30) (S(J),J=NC1,NC2)
          END DO
        ELSE                                    !salvamos invirtiendo espectros
          DO I=NS1,NS2
            DO J=NC1,NC2
              S(J)=A(NC2+NC1-J,I)
            END DO
            WRITE(30) (S(J),J=NC1,NC2)
          END DO
        END IF
        CLOSE(30)
C
        IF(CERR.EQ.'y')THEN
          DO I=NS1,NS2                           !pasamos de varianza a errores
            DO J=NC1,NC2
              ERR(J,I)=SQRT(ERR(J,I))
            END DO
          END DO
          WRITE(*,100)'Output error file name '
          CALL GUESSEF(OUTFILE,ERRFILE)
          ERRFILE=OUTFILEX(31,ERRFILE,NSF,NCF,STWV,DISP,1,.TRUE.)
          IF(CREV.EQ.'n')THEN
            DO I=NS1,NS2
              DO J=NC1,NC2
                S(J)=ERR(J,I)
              END DO
              WRITE(31) (S(J),J=NC1,NC2)
            END DO
          ELSE                                  !salvamos invirtiendo espectros
            DO I=NS1,NS2
              DO J=NC1,NC2
                S(J)=ERR(NC2+NC1-J,I)
              END DO
              WRITE(31) (S(J),J=NC1,NC2)
            END DO
          END IF
          CLOSE(31)
        END IF
C
        STOP
100     FORMAT(A,$)
101     FORMAT(A)
110     FORMAT(A,I6)
        END
C
C******************************************************************************
C Retorna el valor medio de la matriz X(N), eliminando los puntos que se
C alejen de la media mas de TIMES veces sigma. La subrutina permite recuperar
C puntos que han sido eliminados en iteraciones anteriores.
        REAL FUNCTION FMEAN2(N,X,TIMES)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INTEGER N
        REAL X(N)
        REAL TIMES
C
        INTEGER NMAX
        PARAMETER(NMAX=NSMAX*NCMAX)
C
        INTEGER I,NN
        REAL FSIGMA
        DOUBLE PRECISION SUM,SIGMA
        LOGICAL IFX(NMAX),IFXX(NMAX)
        LOGICAL LREPEAT
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        IF(N.EQ.0) STOP 'FATAL ERROR in FMEAN2: N=0.'
        IF(N.GT.NMAX) STOP 'FATAL ERROR in FMEAN2: N too large.'
C
        DO I=1,N
          IFX(I)=.TRUE.
        END DO
C
10      NN=0
        SUM=0.D0
        DO I=1,N
          IF(IFX(I))THEN
            NN=NN+1
            SUM=SUM+DBLE(X(I))
          END IF
        END DO
        FMEAN2=REAL(SUM)/REAL(NN)
        IF(N.EQ.1) RETURN
C
        SIGMA=0.D0
        IF(NN.GT.1)THEN
          DO I=1,N
            IF(IFX(I)) SIGMA=SIGMA+DBLE(X(I)-FMEAN2)*DBLE(X(I)-FMEAN2)
          END DO
          SIGMA=DSQRT(SIGMA/DBLE(NN-1))
        END IF
C
        FSIGMA=REAL(SIGMA)
        DO I=1,N
          IFXX(I)=(ABS(X(I)-FMEAN2).LE.TIMES*FSIGMA)
        END DO
C
        LREPEAT=.FALSE.
        DO I=1,N
          IF(IFX(I).NEQV.IFXX(I)) LREPEAT=.TRUE.
        END DO
        IF(.NOT.LREPEAT) RETURN
C
        DO I=1,N
          IFX(I)=IFXX(I)
        END DO
        GOTO 10
C
        END
