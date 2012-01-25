C------------------------------------------------------------------------------
C Version 8-July-1998                                           file: wcnoarc.f
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
C Program: wcnoarc
C Classification: wavelengths
C Description: Performs the wavelength calibration of an image by using the 
C polynomial fits obtained with fitlin for other frames.
C
Comment
C
C Calibra en l.d.o. un espectro para el cual no tenemos arco. Requiere
C el polinomio de calibracion en l.d.o. de otra imagen (o un promedio),
C un espectro de referencia (a velocidad radial = 0 km/sec), el espectro
C que queremos calibrar, y la velocidad radial de la estrella que queremos
C calibrar en el instante de observacion (geocentrica).
C
C------------------------------------------------------------------------------
C Este programa reune en uno solo las tareas realizadas por:
C shpol
C rebincw (sin tener presente correccion C)
C rvshift
C corrfft
C------------------------------------------------------------------------------
C Para poder calibrar imagenes sin arco a partir de imagenes con arco se puede
C utilizar este programa cuyos detalles esenciales se dan a continuacion.
C
C Partimos de la idea de que los polinomios de calibracion en ldo son
C aproximadamente los mismos para todas las imagenes. Las diferencias entre los
C polinomios se deberan exclusivamente a un desplazamiento (offset) debido a las
C flexiones del telescopio/instrumentos segun la orientacion del telescopio en
C el momento de la observacion. Esto quiere decir que podemos utilizar un
C polinomio cualquiera (obtenido a partir de un arco en el que no hayamos tenido
C que eliminar muchas lineas por ejemplo), e introducir offsets de prueba hasta
C dar con el offset adecuado. El procedimiento a seguir es el siguiente:
C (1) introducimos un offset de prueba
C (2) corregimos los coeficientes del polinomio de dicho offset (shpol)
C (3) calibramos en ldo el espectro problema (rebincw)
C (4) pasamos el espectro problema a velocidad radial cero (rvshift)
C (5) hacemos correlacion cruzada del espectro problema a velocidad radial cero
C  con el de una estrella calibrada con arco (tambien puesta a velocidad radial
C  cero para evitar problemas de distorsion en el espectro debido a velocidades
C  radiales) (corrfft)
C (6) obtenemos un offset en la correlacion
C Si representamos graficamente el offset obtenido en la correlacion cruzada (6)
C en funcion del offset introducido en el polinomio (1), obtenemos una funcion
C que podemos ajustar a un polinomio de grado bajo (generalmente 2). Obteniendo
C la raiz del polinomio (el valor de x que nos da y=0), tenemos un offset que
C introducido en (1) nos calibra en ldo el espectro problema.
C NOTA: en un principio se penso que el la altura del maximo de la correlacion
C seria mayor para el valor del offset inicial que produce offset_corrfft=0.
C Sin embargo se ha visto de forma practica que esto no es asi. De todas
C maneras el programa dibuja tambien la altura del maximo en la correlacion
C cruzada en funcion del offset de entrada inicial.
C------------------------------------------------------------------------------
C
        PROGRAM WCNOARC
C
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER READILIM
        REAL READF
C
        INTEGER NMAX
        PARAMETER (NMAX=8192)      !Si se cambia, cambiar tambien en subrutinas
        INTEGER MAXNINTERV
        PARAMETER (MAXNINTERV=101)               !numero maximo de intervalos+1
        INTEGER MAXNTEMP
        PARAMETER (MAXNTEMP=99)            !numero maximo de espectros template
C
        INTEGER I,J,L
        INTEGER K,IDUM
        INTEGER NDEG              !grado del polinomio de calibracion en l.d.o.
        INTEGER NTERMS                                                 !=NDEG+1
        INTEGER NINTERV    !numero de intervalos para buscar el offset adecuado
        INTEGER NC1,NC2            !primer y ultimo canal a utilizar en CORRFFT
        INTEGER NCHANEFF                                    !NCHANEFF=NC2-NC1+1
        INTEGER K1,K2,K3,K4      !frecuencias para filtrar espectros en CORRFFT
        INTEGER NSIZE                 !tamanho efectivo del espectro en CORRFFT
        INTEGER NSIDE       !No. points at each side of maximum to fit gaussian
        INTEGER NTERMS0           !grado+1 del polinomio final XOFF_FFT vs XOFF
        INTEGER NTEMPT,NTEMP           !Numero de espectros Template a utilizar
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        REAL AA(0:19)      !coeficientes del polinomio de calibracion en l.d.o.
        REAL BB(0:19)              !coeficientes del polinomio traslado en XOFF
        REAL A(20)  !igual que BB pero el indice comienza en 1 en lugar de en 0
        REAL AF(20)                       !coeficientes ajuste XOFF_FFT vs XOFF
        REAL XOFF(MAXNINTERV)       !offset de prueba introducido para calibrar
        REAL XOFF_FFT(MAXNINTERV)       !offset final en la correlacion cruzada
        REAL MAX_FFT(MAXNINTERV)    !valor del maximo en la correlacion cruzada
        REAL SPTEMP(NMAX)     !espectro calibrado y a velocidad radial 0 km/sec
        REAL SP(NCMAX)                !espectro que queremos calibrar en l.d.o.
        REAL SPLDO(NCMAX)                !espectro problema calibrado en l.d.o.
        REAL SPLDORV(NCMAX)   !espectro problema cal. en l.d.o. y a Rv=0 km/sec
        REAL SPLDORVLONG(NMAX)      !igual que SPLDORV pero dimensionado a NMAX
                                    !lo cual hace falta
        REAL RADVEL      !velocidad radial de la estrella que queremos calibrar
        REAL MAXOFF                  !offset maximo para buscar offset adecuado
        REAL FL  !fraccion de espectro sobre la que se aplica la campana coseno
        REAL STWV0,DISP0                                      !variables basura
        REAL KFILTER(NMAX)                 !filtro en el espacio de frecuencias
        REAL COSBELL(NCMAX)         !campana de coseno para correlacion cruzada
        REAL XCORR(NMAX),FCORR(NMAX)                    !funcion de correlacion
        REAL XMIN,XMAX,YMIN,YMAX,DX,DY
        REAL CHISQR,XXX,POLXXX,XX0,XX1,FUN1,FUN2
        REAL XOFF_FINAL(MAXNTEMP)          !offset calculado para cada template
        REAL MEAN,DISPERS                     !media y dispersion de XOFF_FINAL
        CHARACTER*1 CSHOW,CSHOW2,COUT
        CHARACTER*1 CSILENT                                       !run silently
        CHARACTER*75 INFILE0,INFILE1,INFILE2,FILEXOFF
        CHARACTER*80 CDUMMY
        LOGICAL LCOLOR(MAX_ID_RED)
        LOGICAL LSHOW,LSHOW2
C
        COMMON/BLKGEN0/STWV,DISP
        COMMON/BLKGEN1/NCHAN
        COMMON/BLKGEN2/NTERMS
        COMMON/BLKGEN3/A
        COMMON/BLKGEN4/CSILENT
        COMMON/BLKSHPOL1/NDEG
        COMMON/BLKSHPOL2/AA,BB
        COMMON/BLKREBINCW1/SP
        COMMON/BLKREBINCW2/SPLDO
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
C------------------------------------------------------------------------------
C Avoid compilation warnings
        RADVEL=0.0
        NINTERV=0
        MAXOFF=0.0
C------------------------------------------------------------------------------
        THISPROGRAM='wcnoarc'
        CALL WELCOME('8-July-1998')
C------------------------------------------------------------------------------
        WRITE(*,100)'Run silently (y/n) '
        CSILENT(1:1)=READC('y','yn')
C------------------------------------------------------------------------------
C numero de estrellas template a utilizar
        WRITE(*,100)'No. of template spectra to be used '
        NTEMPT=READILIM('1',1,MAXNTEMP)
        NTEMP=0                                       !el proximo es el primero
C------------------------------------------------------------------------------
C abrimos salida grafica
        CALL PIDEGTER(NTERM,IDN,LCOLOR)
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGSCH(1.2)
        END DO
C
        WRITE(*,100)'Show all plots (y/n) '
        CSHOW(1:1)=READC('n','yn')
        LSHOW=(CSHOW.EQ.'y')
        IF(LSHOW)THEN
          CSHOW2='y'
        ELSE
          WRITE(*,100)'Show final correlation plot (y/n) '
          CSHOW2(1:1)=READC('n','yn')
        END IF
        LSHOW2=(CSHOW2.EQ.'y')
C------------------------------------------------------------------------------
        WRITE(*,*)
        WRITE(*,100)'File name (file with generic polynomial'//
     +   ' coefficients)'
        INFILE0=INFILEX(20,'@',0,0,0.,0.,3,.FALSE.)
        K=-1
10      CONTINUE
        READ(20,*,END=11) IDUM,AA(K+1)
        K=K+1
        WRITE(*,*)K,AA(K)
        GOTO 10
11      CONTINUE
        CLOSE(20)
        NDEG=K
        NTERMS=NDEG+1
C------------------------------------------------------------------------------
14      NTEMP=NTEMP+1
        WRITE(*,*)
        WRITE(*,100)'Calibrated template spectrum (Rv=0 km/sec)'
        INFILE1=INFILEX(14,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        IF(NSCAN.GT.1)THEN
          WRITE(*,101)'ERROR: this file contains more than 1 spectrum.'
          CLOSE(14)
          STOP
        END IF
        READ(14) (SPTEMP(J),J=1,NCHAN)
        CLOSE(14)
        WRITE(*,*)
        WRITE(*,101)'* NOTE: these values of STWV and DISP will be '//
     +   'employed in the wavelength '
        WRITE(*,101)'        calibration of the uncalibrated spectrum.'
        IF(NTEMP.GT.1) GOTO 50
C------------------------------------------------------------------------------
        WRITE(*,*)
        WRITE(*,100)'Spectrum to be calibrated'
        INFILE2=INFILEX(15,'@',NSCAN,NCHAN,STWV0,DISP0,1,.FALSE.)
        IF(NSCAN.GT.1)THEN
          WRITE(*,101)'ERROR: this file contains more than 1 spectrum.'
          CLOSE(15)
          STOP
        END IF
        READ(15) (SP(J),J=1,NCHAN)
        CLOSE(15)
C------------------------------------------------------------------------------
        WRITE(*,100)'Radial velocity of last spectrum'
        RADVEL=READF('@')
C------------------------------------------------------------------------------
15      WRITE(*,100)'Maximum offset around initial position (channels) '
        MAXOFF=READF('10.0')
        IF(MAXOFF.LT.0)THEN
          WRITE(*,101)'ERROR: this value must be greater or equal '//
     +     'than zero. Try again.'
          GOTO 15
        END IF
        WRITE(*,100)'Number of intervals '
        NINTERV=READILIM('10',1,MAXNINTERV-1)
C------------------------------------------------------------------------------
C Definimos los parametros necesarios para realizar la correlacion cruzada
C de manera automatica
C..............................................................................
C Definimos region a utilizar de los espectros
20      WRITE(CDUMMY,'(A2,I10)')'1,',NCHAN
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,100)'1st & last channel to be employed '
        CALL READ2I(CDUMMY(1:L),NC1,NC2)
        IF((NC1.LT.1).OR.(NC2.GT.NCHAN).OR.(NC1.GT.NC2))THEN
          WRITE(*,101)'ERROR: invalid numbers. Try again.'
          GOTO 20
        END IF
        NCHANEFF=NC2-NC1+1
C..............................................................................
C Potencia de 2 superior a NCHANEFF
        CALL FFT2POWER(NCHANEFF,NSIZE)
C..............................................................................
C Definimos frecuencias de filtrado
        K1=10
        K2=20
        K3=50
        K4=100
        CALL FFTKFILTER(NSIZE,KFILTER,K1,K2,K3,K4)
C..............................................................................
C Definimos la campana de coseno
        FL=0.10
        CALL FFTCOSBELL(NCHANEFF,COSBELL,FL)
C..............................................................................
C Numero de puntos a cada lado del maximo para encontrar el offset
        WRITE(*,100)'No. of points at each side of maximum to fit '//
     +   'polynomial '
        NSIDE=READILIM('3',1,NSIZE/2)
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C Bucle principal de busqueda del offset adecuado. Recorremos desde -MAXOFF
C hasta +MAXOFF, diviendo dicho intervalo en NINTERV regiones, es decir,
C calculando todo el proceso en NINTERV+1 puntos
50      IF(NCHAN.LT.NSIZE)THEN
          DO J=NCHAN+1,NSIZE
            SPTEMP(J)=0.
          END DO
        END IF
        CALL FFTPREP(NSIZE,SPTEMP,NC1,NC2,.TRUE.,COSBELL,.TRUE.,
     +   KFILTER,LSHOW,INFILE1)
        CALL STOPPLOT(LSHOW)
C
        DO I=1,NINTERV+1
          XOFF(I)=-MAXOFF+REAL(I-1)/REAL(NINTERV)*2.0*MAXOFF
          CALL SUB_SHPOL(XOFF(I))      !generamos un nuevo polinomio trasladado
          DO K=0,NDEG
            A(K+1)=BB(K)             !cambiamos de variables por compatibilidad
          END DO
          CALL SUB_REBINCW         !calibramos en l.d.o. con el nuevo polinomio
          CALL RVREBIN(-RADVEL,NCHAN,SPLDO,SPLDORV,STWV,DISP)       !corregimos
                                                  !espectro de velocidad radial
          DO J=1,NCHAN
            SPLDORVLONG(J)=SPLDORV(J)
          END DO
          IF(NCHAN.LT.NSIZE)THEN
            DO J=NCHAN+1,NSIZE
              SPLDORVLONG(J)=0.
            END DO
          END IF
          CALL FFTPREP(NSIZE,SPLDORVLONG,NC1,NC2,.TRUE.,COSBELL,.TRUE.,
     +     KFILTER,LSHOW,INFILE2)
          CALL STOPPLOT(LSHOW)
          CALL FFTCORREL(NSIZE,SPTEMP,SPLDORVLONG,XCORR,FCORR)
          IF(LSHOW2)THEN
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              CALL PGPAGE
              CALL PGIDEN_RED
              CALL AUTOPLOT(NSIZE,XCORR,FCORR,1,NSIZE,
     +         'X-offset','correlation value',
     +          'Correlation of the two data sets',
     +         .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +         0,.FALSE.,'BCNTS','BCNTS',
     +         101,5,
     +         0.,1.,0.5,1.0)
              CALL PGMTEXT('T',1.5,0.0,0.0,INFILE1)
              CALL PGMTEXT('T',1.5,1.0,1.0,INFILE2)
            END DO
          END IF
          CALL FFTCORRZOOM(XCORR,FCORR,NSIZE,NSIZE-(NC2-NC1+1),
     +     INFILE1,INFILE2,NSIDE,LSHOW2,XOFF_FFT(I),MAX_FFT(I))
          CALL STOPPLOT(LSHOW2)
        END DO
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C Dibujamos el resultado
        XMIN=-MAXOFF
        XMAX=MAXOFF
        DX=XMAX-XMIN
        XMIN=XMIN-DX/50.
        XMAX=XMAX+DX/50.
C
        CALL FINDMM(NINTERV+1,MAX_FFT,YMIN,YMAX)
        DY=YMAX-YMIN
        YMIN=YMIN-DY/50.
        YMAX=YMAX+DY/50.
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,0)
          CALL PGPOINT(NINTERV+1,XOFF,MAX_FFT,3)
          CALL PGLINE(NINTERV+1,XOFF,MAX_FFT)
          CALL PGLABEL('XOFF','MAX_FFT',CHAR(32))
          CALL PGMTEXT('T',1.,1.,1.,INFILE1)
          CALL PGMTEXT('T',1.,0.,0.,INFILE2)
        END DO
C
        WRITE(*,100)'Press <RETURN> to continue...'
        READ(*,*)
C
        CALL FINDMM(NINTERV+1,XOFF_FFT,YMIN,YMAX)
        DY=YMAX-YMIN
        YMIN=YMIN-DY/50.
        YMAX=YMAX+DY/50.
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,0)
          CALL PGPOINT(NINTERV+1,XOFF,XOFF_FFT,3)
          CALL PGLINE(NINTERV+1,XOFF,XOFF_FFT)
          CALL PGLABEL('XOFF','XOFF_FFT',CHAR(32))
          CALL PGMTEXT('T',1.,1.,1.,INFILE1)
          CALL PGMTEXT('T',1.,0.,0.,INFILE2)
        END DO
        WRITE(*,100)'Polynomial degree to fit data '
        NTERMS0=READILIM('2',0,19)
        NTERMS0=NTERMS0+1
        CALL POLFIT(XOFF,XOFF_FFT,XOFF_FFT,NINTERV+1,NTERMS0,0,AF,
     +   CHISQR)
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM)) CALL PGSCI(2)
          DO K=1,201                          !dibujamos el polinomio calculado
            XXX=XMIN+REAL(K-1)/200.*(XMAX-XMIN)
            POLXXX=AF(NTERMS0)
            DO J=NTERMS0-1,1,-1
              POLXXX=POLXXX*XXX+AF(J)
            END DO
            IF(K.EQ.1)THEN
              CALL PGMOVE(XXX,POLXXX)
            ELSE
              CALL PGDRAW(XXX,POLXXX)
            END IF
          END DO
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
        END DO
C Calculamos el cero por Newton-Raphson
        XX0=0.                                 !como valor inicial supongo cero
70      CONTINUE
        FUN1=AF(NTERMS0)                        !calculamos el polinomio en XX0
        DO J=NTERMS0-1,1,-1
          FUN1=FUN1*XX0+AF(J)
        END DO
        FUN2=REAL(NTERMS0-1)*AF(NTERMS0)          !calculamos la derivada en XX0
        DO J=NTERMS0-1,2,-1
          FUN2=FUN2*XX0+REAL(J-1)*AF(J)
        END DO
        XX1=XX0-FUN1/FUN2
        IF(ABS(XX1-XX0).GT.1.0E-4)THEN
          XX0=XX1
          GOTO 70
        END IF
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM)) CALL PGSCI(3)
          CALL PGSLS(2)
          CALL PGMOVE(XX1,YMIN)
          CALL PGDRAW(XX1,YMAX)
          CALL PGSLS(1)
          CALL PGMOVE(XMIN,0.)
          CALL PGDRAW(XMAX,0.)
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
        END DO
C------------------------------------------------------------------------------
        WRITE(*,101)'Fit result:'
        DO J=1,NTERMS0
          WRITE(*,*)J,AF(J)
        END DO
        WRITE(*,100)'Final XOFF = '
        WRITE(*,*)XX1
        XOFF_FINAL(NTEMP)=XX1
        IF(NTEMP.LT.NTEMPT) GOTO 14
C------------------------------------------------------------------------------
        CALL PGEND
C------------------------------------------------------------------------------
        IF(NTEMPT.EQ.1)THEN
          MEAN=XOFF_FINAL(1)
          GOTO 90
        END IF
        WRITE(*,*)
        MEAN=0.
        DO K=1,NTEMPT
          WRITE(*,*)K,XOFF_FINAL(K)
          MEAN=MEAN+XOFF_FINAL(K)
        END DO
        MEAN=MEAN/REAL(NTEMPT)
        DISPERS=0.
        DO K=1,NTEMPT
          DISPERS=DISPERS+(MEAN-XOFF_FINAL(K))*(MEAN-XOFF_FINAL(K))
        END DO
        DISPERS=SQRT(DISPERS/REAL(NTEMPT-1))
        WRITE(*,*)
        WRITE(*,100)'>>> Mean......: '
        WRITE(*,*)MEAN
        WRITE(*,100)'>>> Dispersion: '
        WRITE(*,*)DISPERS
C------------------------------------------------------------------------------
90      WRITE(*,100)'Write final XOFF into file (y/n) '
        COUT(1:1)=READC('n','yn')
        IF(COUT.EQ.'y')THEN
          WRITE(*,100)'Output file name for XOFF'
          FILEXOFF=OUTFILEX(37,'@',0,0,0.,0.,3,.FALSE.)
          WRITE(37,*)MEAN
          CLOSE(37)
        END IF
C------------------------------------------------------------------------------
        STOP
100     FORMAT(A,$)
101     FORMAT(A)
        END
C
C******************************************************************************
C******************************************************************************
C Calibracion en longitud de onda del espectro SP(). La salida se produce en
C SPLDO(). Esta subrutina esta adaptada del programa rebincw.f
        SUBROUTINE SUB_REBINCW
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
C
        INTEGER I,J,IB0,J1,J2,JMAX
        INTEGER NTERMS
        REAL SP(NCMAX)
        REAL SPLDO(NCMAX)
        REAL A(20)
        REAL WAV,SUM,XA,XB,X1,X1M,X2M
        REAL F1,F2,A1,A2,A0,PIXELM,ECH
        REAL ST,ST2,MAXERROR
        DOUBLE PRECISION PIXEL,PIXEL0,Q,P,ERW,XACC
        DOUBLE PRECISION PIXLAM
        CHARACTER*1 CSILENT
C
        COMMON/BLKGEN0/STWV,DISP
        COMMON/BLKGEN1/NCHAN
        COMMON/BLKGEN2/NTERMS
        COMMON/BLKGEN3/A
        COMMON/BLKGEN4/CSILENT
        COMMON/BLKREBINCW1/SP
        COMMON/BLKREBINCW2/SPLDO
        COMMON/LOCREBINCW1/JMAX
        COMMON/LOCREBINCW2/XACC
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        DO I=1,NCHAN                       !inicializamos el espectro de salida
          SPLDO(I)=0.0
        END DO
C------------------------------------------------------------------------------
        XACC=1.D-10
        JMAX=200
C
C   FROM VISTA:
C       To rebin the old spectrum, we calculate the extreme lambdas of each
C       bin in the new spectrum. The location of these lambdas in pixel space
C       of the original spectrum is calculated using the inverse dispersion.
C       Whole pixels lying inside the new bin are directly added. Any fractional
C       pixel contribution to the bin is estimated by integration of the
C       parabola whose area matches the intensity of the pixel and its two
C       neighbors (Gonzalez-integration).
C       NOTE: this is different from Sympson-integration which assumes the
C       function is only sampled (not binned and sampled). The Sympson-parabola
C       C0+C1*X+C2*X*X passing through (-1,Ym),(0,Y0),(1,Yp), has coefficients
C       C0=Y0, C1=(Yp-Ym)/2, and C2=(Yp+Ym)/2-Y0, while the Gonzalez-parabola
C       A0+A1*X+A2*X*X, has A0=C0-(C2/12), A1=C1, and A2=C2. This scheme takes
C       care of the fact that the original spectrum is also binned data.
C
        DO I=1,NCHAN
          IB0=0 
C-- Find the limits of the new-pixel in the old pixel scale:
          IF (I.GT.1) THEN
            PIXEL=PIXEL0
          ELSE
            WAV = STWV + (FLOAT(I)-1.5)*DISP
            PIXEL = PIXLAM(WAV,0.D0,I)
          END IF
          X1 = SNGL(PIXEL)
          PIXEL0=PIXEL
          WAV = STWV + (FLOAT(I)-0.5)*DISP
          PIXEL = PIXLAM(WAV,PIXEL0,I)
          Q=PIXEL
          P=DBLE(A(NTERMS))
          DO J=NTERMS-1,1,-1
            P=P*Q+DBLE(A(J))
          END DO
          ERW=DABS(DBLE(WAV)-P)
          IF(ERW.GT.5.D-12) THEN
            WRITE(*,*)'Warning: Error in wavelength assignations'
            WRITE(*,'(A,I5,A,$)')'Pixel ',I,'    Diff='
            WRITE(*,*)ERW
            STOP
          END IF
          PIXEL0=PIXEL
C
          SUM=0.0
          XA=X1
          XB=SNGL(PIXEL)  
          X1M=X1
          PIXELM=SNGL(PIXEL)
          X2M = MAX(0.5,MIN(NCHAN+0.5,MAX(PIXELM,X1M)))
          X1M = MAX(0.5,MIN(NCHAN+0.5,MIN(PIXELM,X1M)))
          IF (X1M.EQ.X2M) THEN                  !The new pixel maps outside the
            SPLDO(I) = 0.0E0                    !range of the old spectrum.
            IB0=I
            GO TO 200
          END IF
          J1 = MIN0(NCHAN-1,MAX0(NINT(X1M),2))    ! This takes care of the
          J2 = MAX0(2,MIN0(NINT(X2M),NCHAN-1))    ! interpolation near edges.
          F1 = X1M - FLOAT(J1)
          F2 = X2M - FLOAT(J2)
          A1 = (SP(J1+1)-SP(J1-1))/4.
          A2 = (SP(J1+1)+SP(J1-1))/6. - SP(J1)/3.
          IF (J2 .EQ. J1) THEN
C           A0 = SP(J1) + A2/4.0                        ! Binned parabola.
            A0 = SP(J1) - A2/4.
            SUM = (F2-F1)*(A0+(F1+F2)*A1+A2*(F1*F1+F2*F1+F2*F2))
          ELSE
C-- Left-pixel fractional contribution to flux:
            SUM = (0.5-F1)*(SP(J1)+(F1+0.5)*(A1+A2*F1))
C-- Whole pixel contribution to flux:
            IF (J1.LE.J2-2) THEN                        ! This check was put
              DO J=J1+1,J2-1,1                  ! here for portability.
                SUM = SUM + SP(J)
              END DO
            END IF
C-- Right-pixel fractional contribution to flux:
            A1 = (SP(J2+1)-SP(J2-1))/4.
            A2 = (SP(J2+1)+SP(J2-1))/6. - SP(J2)/3.
            SUM = SUM + (0.5+F2)*(SP(J2)+(F2-0.5)*(A1+A2*F2))
          END IF
          SPLDO(I)=SUM
C         SPA(I) = SUM/(X2-X1)          ! Mean flux inside the new bin.
200       CONTINUE
          IF(CSILENT.EQ.'n')THEN
            IF(IB0.NE.0) WRITE(*,'(3(A,I5))')'Pixel ',I,'  -------> '//
     +       'Blank'
          END IF
        END DO
C------------------------------------------------------------------------------
C Conservacion del numero de cuentas
        IF(CSILENT.EQ.'n')THEN
          WRITE(*,*)
          WRITE(*,100)'>>> Checking conservation of number of counts...'
        END IF
        MAXERROR=0.
        ST=0.
        ST2=0.
        DO I=1,NCHAN
          ST=ST+SP(I)
          ST2=ST2+SPLDO(I)
        END DO
        IF(ST.EQ.0.)THEN
          WRITE(*,101)'ERROR: No. of counts = 0.'
        ELSE
          ECH=ABS(ST-ST2)/ST*100
          IF(ECH.GT.0.5) THEN
            IF(ABS(ST-ST2)/ST*100. .GT. MAXERROR) 
     +       MAXERROR=ABS(ST-ST2)/ST*100.
          END IF
        END IF
        IF(CSILENT.EQ.'n')THEN
          WRITE(*,101)'   ...OK!'
          WRITE(*,100)'>>> Maximum error (%): '
          WRITE(*,*)MAXERROR
        END IF
C------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
        END
C
C******************************************************************************
C******************************************************************************
C Funcion complementaria de SUB_REBINCW
        DOUBLE PRECISION FUNCTION PIXLAM(WAV,X0,II)
        IMPLICIT NONE
        INTEGER II,J,JMAX,NTERMS
        REAL WAV,A(20)
        DOUBLE PRECISION XP,X1,X2,PRO,PRO2,DPRO,X0,XACC,RTNEWT,F,DF
        DOUBLE PRECISION DX,DX0
        COMMON/BLKGEN2/NTERMS
        COMMON/BLKGEN3/A
        COMMON/LOCREBINCW1/JMAX
        COMMON/LOCREBINCW2/XACC
C------------------------------------------------------------------------------
        pixlam=0.d0    !ncl: evita un Warning que no tiene sentido
        XP=X0
10      CALL FUNCD(XP,PRO,DPRO,WAV)
        IF(PRO.GT.0.D00) THEN
          XP=XP-1.D00
          GO TO 10
        ELSE
          X1=XP
        END IF
        X2=XP
20      X2=X2+2.D00
        CALL FUNCD(X2,PRO2,DPRO,WAV)
        IF(PRO2.LE.0.D00) GO TO 20
C
C       WRITE(*,*)X1,X2,PRO,PRO2
        RTNEWT=.5D00*(X1+X2)
        DO J=1,JMAX
          CALL FUNCD(RTNEWT,F,DF,WAV)
          DX=F/DF
          IF(J.GT.100) WRITE(*,*)'DX',DX,DX0,DX+DX0
          IF(DABS(DX+DX0).LT.(XACC*1.D-10)) THEN
            WRITE(*,*)II
            WRITE(*,*)'DX',DX,DX0,DX+DX0
            DX=DX/2.D00
            WRITE(*,*)'OOOHO'
          END IF
          DX0=DX
C         WRITE(*,*)RTNEWT,F,DX
          RTNEWT=RTNEWT-DX
          IF((X1-RTNEWT)*(RTNEWT-X2).LT.0.)THEN
            WRITE(*,*)PRO,PRO2
            WRITE(*,*)X1,X2,RTNEWT,WAV
            STOP 'Jumped out of brackets'
          END IF
          IF(ABS(DX).LT.XACC) THEN
            PIXLAM=RTNEWT
            RETURN
          END IF
        END DO
        WRITE(*,*)X1,X2,RTNEWT,WAV
        STOP 'PIXLAM exceeding maximum iterations'
        END
C
C******************************************************************************
C******************************************************************************
C Funcion complementaria de SUB_REBINCW
        SUBROUTINE FUNCD(X,F,DF,WAV)
        IMPLICIT NONE
        INTEGER J,NTERMS,JMAX
        REAL WAV,A(20)
        DOUBLE PRECISION X,F,DF,P,DP,XACC
        COMMON/BLKGEN2/NTERMS
        COMMON/BLKGEN3/A
        COMMON/LOCREBINCW1/JMAX
        COMMON/LOCREBINCW2/XACC
C------------------------------------------------------------------------------
        P=DBLE(A(NTERMS))
        DP=0.D00
        DO J=NTERMS-1,1,-1
          DP=DP*X+P
          P=P*X+DBLE(A(J))
        END DO
        F=P-DBLE(WAV)
        DF=DP
        RETURN
        END
C
C******************************************************************************
C******************************************************************************
C Traslada un polinomio un cierto offset XOFF en el eje x. Los coeficientes 
C del polinomio de entrada son A() y los de salida B().
        SUBROUTINE SUB_SHPOL(XOFF)
        IMPLICIT NONE
        REAL XOFF
C
        INTEGER I,J
        INTEGER NDEG              !grado del polinomio de calibracion en l.d.o.
        REAL AA(0:19)      !coeficientes del polinomio de calibracion en l.d.o.
        REAL BB(0:19)              !coeficientes del polinomio traslado en XOFF
        REAL COMB             !funcion numero combinatorio (definida mas abajo)
C
        COMMON/BLKSHPOL1/NDEG
        COMMON/BLKSHPOL2/AA,BB
C------------------------------------------------------------------------------
        DO I=0,NDEG
          BB(I)=0.
        END DO
        DO I=0,NDEG
          DO J=0,I
            BB(I-J)=BB(I-J)+AA(I)*COMB(I,J)*(XOFF**REAL(J))
          END DO
        END DO
        END
C
C******************************************************************************
C******************************************************************************
C Calcula el numero combinatorio N sobre M
        REAL FUNCTION COMB(N,M)
        IMPLICIT NONE
        INTEGER N,M
        REAL FACT                       !funcion factorial (definida mas abajo)
C------------------------------------------------------------------------------
        COMB=FACT(N)/(FACT(M)*FACT(N-M))
        END
C
C******************************************************************************
C******************************************************************************
C Calcula el factorial de K
        REAL FUNCTION FACT(K)
        IMPLICIT NONE
        INTEGER K
        INTEGER J
C------------------------------------------------------------------------------
        FACT=1.
        IF(K.LT.0)THEN
          STOP 'Factorial of negative number.'
        ELSEIF(K.GT.0)THEN
          DO J=1,K
            FACT=FACT*REAL(J)
          END DO
        END IF
        END
C
C******************************************************************************
C******************************************************************************
C Determina si seguimos o no dibujando
        SUBROUTINE STOPPLOT(LPLOT)
        IMPLICIT NONE
        INCLUDE 'futils.inc'
C
        LOGICAL LPLOT
C
        CHARACTER*1 CMORE
C------------------------------------------------------------------------------
        IF(LPLOT)THEN
          WRITE(*,100)'More plots (y/n) '
          CMORE(1:1)=READC('y','yn')
          LPLOT=(CMORE.EQ.'y')
        END IF
C
100     FORMAT(A,$)
        END
