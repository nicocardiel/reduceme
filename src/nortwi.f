C------------------------------------------------------------------------------
C Version 7-July-1998                                            file: nortwi.f
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
C Program: nortwi
C Classification: distortion
C Description: Wavelength normalization of twilight flatfields calculating 
C C-distortion with cross correlation.
C
Comment
C
C Normaliza un flatfield de twilight, calculando el offset (respecto a un
C espectro promedio) de cada espectro de la imagen antes de dividir por un
C espectro promedio.
C
C Modus Operandi: se define una region para extraer el espectro promedio
C (preferiblemente cerca del centro de la imagen y sin utilizar demasiados
C scans). Luego definimos el intervalo en l.d.o. (channels) a utilizar para
C realizar la correlacion cruzada. No conviene utilizar todo el espectro
C porque sino el programa puede ser muy lento. Luego, el programa calcula
C el offset (mediante la correlacion cruzada) entre el espectro promedio y
C cada uno de los espectros individuales de la imagen. El resultado permite
C mover el espectro promedio para que al dividirlo, no tengamos ningun
C problema de distorsion C.
C El programa tambien permite utilizar la distorsion calculada para refinar
C el calculo, empleando un nuevo espectro promedio calculado sumandos mas
C SCANS (corrigiendo apropiadamente de distorsion), y reptiendo el proceso.
C Asimismo el programa puede generar como salida una imagen igual a la
C original pero corregida de la distorsion calculada.
C
        PROGRAM NORTWI
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER READILIM
C
        INTEGER NMAX
        PARAMETER (NMAX=8192)      !Si se cambia, cambiar tambien en subrutinas
        REAL PI
        PARAMETER (PI=3.141593)
C
        INTEGER I,J,L
        INTEGER K1,K2,K3,K4
        INTEGER NS1,NS2,NS0
        INTEGER NC1,NC2
        INTEGER NCHANEFF
        INTEGER NSIZE,K,NSIDE
        INTEGER NTERMS,NPT
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        REAL A(NCMAX,NSMAX),B(NCMAX,NSMAX)
        REAL X(NMAX)
        REAL MEANSP(NMAX),TESTSP(NMAX)
        REAL MEANY(NSMAX)
        REAL FL,COSBELL(NCMAX),KFILTER(NMAX)
        REAL XOFF_FFT(NSMAX),MAX_FFT(NSMAX)
        REAL XOFF_POL(NSMAX)
        REAL XMIN,XMAX,YMIN,YMAX,DX,DY
        REAL XMINSC,XMAXSC,YMINSC,YMAXSC
        REAL XPOL(NSMAX),YPOL(NSMAX),CHISQR,XX,YY,AA(20)
        REAL SPP1(NCMAX),SPP2(NCMAX)
        REAL XCORR(NMAX),FCORR(NMAX)
        CHARACTER*1 CSHOW,CFITOK,CFITALL,COPC,CSURE
        CHARACTER*75 INFILE,OUTFILE
        CHARACTER*80 CDUMMY
        LOGICAL LLOOP,LSHOW
        LOGICAL IFSCAN_MEAN(NSMAX)
        LOGICAL IFSCAN_UTIL(NSMAX)
        LOGICAL IFSCAN_POLI(NSMAX)
        LOGICAL LCOLOR(MAX_ID_RED)
C
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
C------------------------------------------------------------------------------
        THISPROGRAM='nortwi'
        CALL WELCOME('7-July-1998')
C abrimos salida grafica
        CALL PIDEGTER(NTERM,IDN,LCOLOR)
        WRITE(*,100)'Show all plots (y/n) '
        CSHOW(1:1)=READC('n','yn')
        LSHOW=(CSHOW.EQ.'y')
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGSCH(1.2)
        END DO
C------------------------------------------------------------------------------
C leemos imagen
        WRITE(*,100)'Input file name'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        IF(NSCAN.EQ.1)THEN
          WRITE(*,101)'ERROR: this file only contains 1 spectrum.'
          CLOSE(20)
          STOP
        END IF
        DO I=1,NSCAN
          READ(20) (A(J,I),J=1,NCHAN)
        END DO
        CLOSE(20)
C
        DO J=1,NMAX
          X(J)=REAL(J)
        END DO
C
        IF((NCMAX.GT.NMAX).OR.(NSMAX.GT.NMAX))THEN
          WRITE(*,101)'FATAL ERROR in program nortwi: redim NMAX'
          STOP
        END IF
C------------------------------------------------------------------------------
C calculamos el corte espacial promedio y lo dibujamos
        DO I=1,NSCAN
          MEANY(I)=0.
          DO J=1,NCHAN
            MEANY(I)=MEANY(I)+A(J,I)
          END DO
          MEANY(I)=MEANY(I)/REAL(NCHAN)
        END DO
        XMIN=1.
        XMAX=REAL(NSCAN)
        DX=XMAX-XMIN
        XMIN=XMIN-DX/30.
        XMAX=XMAX+DX/30.
        CALL FINDMM(NSCAN,MEANY,YMIN,YMAX)
        DY=YMAX-YMIN
        YMIN=YMIN-DY/30.
        YMAX=YMAX+DY/30.
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGSLW(3)
          CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,0)
          CALL PGLABEL('Scan','Averaged no. of Counts',
     +     'File: '//INFILE)
          CALL PGSLW(1)
          CALL PGIDEN_RED
          CALL PGBIN(NSCAN,X,MEANY,.TRUE.)
        END DO
        XMINSC=XMIN
        XMAXSC=XMAX
        YMINSC=YMIN
        YMAXSC=YMAX
C------------------------------------------------------------------------------
C determinamos regiones a utilizar para calcular MEAN spectrum
        LLOOP=.TRUE.
        DO I=1,NSCAN
          IFSCAN_MEAN(I)=.FALSE.
        END DO
        DO WHILE(LLOOP)
          WRITE(*,100)'Scan region to calculate mean spectrum '//
     +     '(0,0=EXIT) '
          CALL READ2I('0,0',NS1,NS2)
          IF((NS1.EQ.0).AND.(NS2.EQ.0))THEN
            NS0=0
            DO I=1,NSCAN
              IF(IFSCAN_MEAN(I)) NS0=NS0+1
            END DO
            IF(NS0.LT.1)THEN
              WRITE(*,101)'WARNING: no. of scans = 0. Try again.'
            ELSE
              LLOOP=.FALSE.
            END IF
          ELSEIF((NS1.LT.1).OR.(NS2.GT.NSCAN).OR.(NS1.GT.NS2))THEN
            WRITE(*,101)'WARNING: numbers out of range. Try again.'
          ELSE
            DO I=NS1,NS2
              IFSCAN_MEAN(I)=.TRUE.
            END DO
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              IF(LCOLOR(ITERM)) CALL PGSCI(3)
              CALL PGMOVE(REAL(NS1),YMIN)
              CALL PGDRAW(REAL(NS1),YMAX)
              CALL PGMOVE(REAL(NS2),YMIN)
              CALL PGDRAW(REAL(NS2),YMAX)
              CALL PGRECT(REAL(NS1),REAL(NS2),YMIN,YMIN+(YMAX-YMIN)*.01)
              CALL PGRECT(REAL(NS1),REAL(NS2),YMAX,YMAX-(YMAX-YMIN)*.01)
              IF(LCOLOR(ITERM)) CALL PGSCI(1)
            END DO
          END IF
        END DO
C------------------------------------------------------------------------------
C calculamos espectro promedio
        DO J=1,NCHAN
          MEANSP(J)=0.
        END DO
        NS0=0
        DO I=1,NSCAN
          IF(IFSCAN_MEAN(I))THEN
            NS0=NS0+1
            DO J=1,NCHAN
              MEANSP(J)=MEANSP(J)+A(J,I)
            END DO
          END IF
        END DO
        DO J=1,NCHAN
          MEANSP(J)=MEANSP(J)/REAL(NS0)
        END DO
C------------------------------------------------------------------------------
C dibujamos espectro promedio
        XMIN=1.
        XMAX=REAL(NCHAN)
        DX=XMAX-XMIN
        XMIN=XMIN-DX/30.
        XMAX=XMAX+DX/30.
        CALL FINDMM(NCHAN,MEANSP,YMIN,YMAX)
        DY=YMAX-YMIN
        YMIN=YMIN-DY/30.
        YMAX=YMAX+DY/30.
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGSLW(3)
          CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,0)
          CALL PGLABEL('Channel','Averaged no. of Counts',
     +     'File: '//INFILE)
          CALL PGSLW(1)
          CALL PGIDEN_RED
          CALL PGBIN(NCHAN,X,MEANSP,.TRUE.)
        END DO
C------------------------------------------------------------------------------
        WRITE(*,*)
        WRITE(*,101)'* Define channel region to calculate cross '//
     +   'correlation:'
        LLOOP=.TRUE.
        DO WHILE(LLOOP)
          WRITE(*,100)'Channel region'
          CALL READ2I('@',NC1,NC2)
          IF((NC1.LT.1).OR.(NC2.GT.NCHAN).OR.(NC1.GT.NC2))THEN
            WRITE(*,101)'ERROR: numbers out of range. Try again.'
          ELSE
            LLOOP=.FALSE.
          END IF
        END DO
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM)) CALL PGSCI(3)
          CALL PGMOVE(REAL(NC1),YMIN)
          CALL PGDRAW(REAL(NC1),YMAX)
          CALL PGMOVE(REAL(NC2),YMIN)
          CALL PGDRAW(REAL(NC2),YMAX)
          CALL PGRECT(REAL(NC1),REAL(NC2),YMIN,YMIN+(YMAX-YMIN)*0.01)
          CALL PGRECT(REAL(NC1),REAL(NC2),YMAX,YMAX-(YMAX-YMIN)*0.01)
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
        END DO
C------------------------------------------------------------------------------
C determinamos regiones a utilizar para calcular el offset mediante 
C correlacion cruzada
        IF(LSHOW)THEN
          WRITE(*,100)'Press <CR> to continue...'
          READ(*,*)
        END IF
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGPAGE
          CALL PGWINDOW(XMINSC,XMAXSC,YMINSC,YMAXSC)
          CALL PGSLW(3)
          CALL PGBOX('BCTNS',0.0,0,'BTNS',0.0,0)
          CALL PGBOX(' ',0.0,0,'C',0.0,0)
          CALL PGLABEL('Scan','Averaged no. of Counts','File: '//INFILE)
          CALL PGSLW(1)
          CALL PGIDEN_RED
          CALL PGBIN(NSCAN,X,MEANY,.TRUE.)
        END DO
C
        WRITE(*,*)
        WRITE(*,101)'* Define scan region to calculate cross '//
     +   'correlation:'
        LLOOP=.TRUE.
        DO I=1,NSCAN
          IFSCAN_UTIL(I)=.FALSE.
        END DO
        DO WHILE(LLOOP)
          WRITE(*,100)'Scan region (0,0=EXIT) '
          CALL READ2I('0,0',NS1,NS2)
          IF((NS1.EQ.0).AND.(NS2.EQ.0))THEN
            NS0=0
            DO I=1,NSCAN
              IF(IFSCAN_UTIL(I)) NS0=NS0+1
            END DO
            IF(NS0.LT.1)THEN
              WRITE(*,101)'WARNING: no. of scans = 0. Try again.'
            ELSE
              LLOOP=.FALSE.
            END IF
          ELSEIF((NS1.LT.1).OR.(NS2.GT.NSCAN).OR.(NS1.GT.NS2))THEN
            WRITE(*,101)'WARNING: numbers out of range. Try again.'
          ELSE
            DO I=NS1,NS2
              IFSCAN_UTIL(I)=.TRUE.
            END DO
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              IF(LCOLOR(ITERM)) CALL PGSCI(3)
              CALL PGMOVE(REAL(NS1),YMINSC)
              CALL PGDRAW(REAL(NS1),YMAXSC)
              CALL PGMOVE(REAL(NS2),YMINSC)
              CALL PGDRAW(REAL(NS2),YMAXSC)
              CALL PGRECT(REAL(NS1),REAL(NS2),
     +         YMINSC,YMINSC+(YMAXSC-YMINSC)*0.01)
              CALL PGRECT(REAL(NS1),REAL(NS2),
     +         YMAXSC,YMAXSC-(YMAXSC-YMINSC)*0.01)
              IF(LCOLOR(ITERM)) CALL PGSCI(1)
            END DO
          END IF
        END DO
C------------------------------------------------------------------------------
C Definimos los parametros necesarios para realizar la correlacion cruzada
C de manera automatica
C..............................................................................
C region a utilizar de los espectros (ya definida)
        NCHANEFF=NC2-NC1+1
C..............................................................................
C potencia de 2 superior a NCHANEFF
        CALL FFT2POWER(NCHANEFF,NSIZE)
C..............................................................................
C definimos frecuencias de filtrado
        K1=10
        K2=20
        K3=50
        K4=100
        CALL FFTKFILTER(NSIZE,KFILTER,K1,K2,K3,K4)
C..............................................................................
C definimos la campana de coseno
        FL=0.10
        CALL FFTCOSBELL(NCHANEFF,COSBELL,FL)
C..............................................................................
C Numero de puntos a cada lado del maximo para encontrar el offset
        WRITE(*,*)
        WRITE(*,100)'No. of points at each side of maximum to fit '//
     +   'polynomial '
        NSIDE=READILIM('3',1,NSIZE/2)
C------------------------------------------------------------------------------
C bucle principal: correlacion cruzada
40      IF(NCHAN.LT.NSIZE)THEN
          DO J=NCHAN+1,NSIZE
            MEANSP(J)=0.
          END DO
        END IF
        CALL FFTPREP(NSIZE,MEANSP,NC1,NC2,.TRUE.,COSBELL,.TRUE.,
     +   KFILTER,LSHOW,'mean spectrum')
        CALL STOPPLOT(LSHOW)
        WRITE(*,100)'Working...:#####'
        DO I=1,NSCAN
          X(I)=REAL(I)
          IF(IFSCAN_UTIL(I))THEN
            WRITE(*,'(A,I5,$)')'\\b\\b\\b\\b\\b',I
            DO J=1,NCHAN
              TESTSP(J)=A(J,I)
            END DO
            IF(NCHAN.LT.NSIZE)THEN
              DO J=NCHAN+1,NSIZE
                TESTSP(J)=0.
              END DO
            END IF
            WRITE(CDUMMY,'(A,I10)')'Scan #',I
            CALL RMBLANK(CDUMMY,CDUMMY,L)
            CALL FFTPREP(NSIZE,TESTSP,NC1,NC2,.TRUE.,COSBELL,.TRUE.,
     +       KFILTER,LSHOW,CDUMMY(1:L))
            CALL STOPPLOT(LSHOW)
            CALL FFTCORREL(NSIZE,MEANSP,TESTSP,XCORR,FCORR)
            IF(LSHOW)THEN
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                CALL PGPAGE
                CALL PGIDEN_RED
                CALL AUTOPLOT(NSIZE,XCORR,FCORR,1,NSIZE,
     +           'X-offset','correlation value',
     +            'Correlation of the two data sets',
     +           .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +           0,.FALSE.,'BCNTS','BCNTS',
     +           101,5,
     +           0.,1.,0.5,1.0)
                CALL PGMTEXT('T',1.5,0.0,0.0,'mean spectrum')
                CALL PGMTEXT('T',1.5,1.0,1.0,CDUMMY(1:L))
              END DO
            END IF
            CALL FFTCORRZOOM(XCORR,FCORR,NSIZE,NSIZE-(NC2-NC1+1),
     +       'mean spectrum',CDUMMY(1:L),NSIDE,LSHOW,
     +       XOFF_FFT(I),MAX_FFT(I))
            CALL STOPPLOT(LSHOW)
          ELSE
            XOFF_FFT(I)=0.
            MAX_FFT(I)=0.
          END IF
        END DO
        WRITE(*,*)
C
50      K=0
        DO I=1,NSCAN
          IF(IFSCAN_UTIL(I))THEN
            K=K+1
            XPOL(K)=X(I)
            YPOL(K)=XOFF_FFT(I)
          END IF
        END DO
        NPT=K
C------------------------------------------------------------------------------
C dibujamos los offsets calculados respecto al espectro promedio
        CALL FINDMM(NPT,YPOL,YMIN,YMAX)
        DY=YMAX-YMIN
        YMIN=YMIN-DY/10.
        YMAX=YMAX+DY/10.
C
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGPAGE
          CALL PGENV(XMINSC,XMAXSC,YMINSC,YMAXSC,0,-2)
          CALL PGSLW(3)
          CALL PGBOX('BCTNS',0.0,0,'BTNS',0.0,0)
          CALL PGBOX(' ',0.0,0,'C',0.0,0)
          CALL PGLABEL('Scan','Averaged no. of Counts','File: '//INFILE)
          CALL PGSLW(1)
          CALL PGIDEN_RED
          CALL PGBIN(NSCAN,X,MEANY,.TRUE.)
          CALL PGWINDOW(XMINSC,XMAXSC,YMIN,YMAX)
          IF(LCOLOR(ITERM)) CALL PGSCI(2)
          CALL PGSLW(3)
          CALL PGBOX(' ',0.0,0,'CTMS',0.0,0)
          CALL PGMTEXT('R',2.5,0.5,0.5,'Deviation (channels)')
          CALL PGSLW(1)
          CALL PGPOINT(NPT,XPOL,YPOL,17)
        END DO
C------------------------------------------------------------------------------
C ajustamos polinomio a los offsets calculados y lo dibujamos
        WRITE(*,101)'* Fitting polynomial to the deviations:'
        WRITE(*,110)'>>> Initial no. of points to be fitted: ',NPT
        WRITE(*,100)'Are you fitting all the points (y/n) '
        CFITALL(1:1)=READC('y','yn')
        IF(CFITALL.EQ.'n')THEN
          WRITE(*,101)'* Define scan region to be employed:'
          LLOOP=.TRUE.
          DO I=1,NSCAN
            IFSCAN_POLI(I)=.FALSE.
          END DO
          DO WHILE(LLOOP)
            WRITE(*,100)'Scan region (0,0=EXIT) '
            CALL READ2I('0,0',NS1,NS2)
            IF((NS1.EQ.0).AND.(NS2.EQ.0))THEN
              NS0=0
              DO I=1,NSCAN
                IF(IFSCAN_POLI(I)) NS0=NS0+1
              END DO
              IF(NS0.LT.1)THEN
                WRITE(*,101)'WARNING: no. of scans = 0. Try again.'
              ELSE
                LLOOP=.FALSE.
              END IF
            ELSEIF((NS1.LT.1).OR.(NS2.GT.NSCAN).OR.(NS1.GT.NS2))THEN
              WRITE(*,101)'WARNING: numbers out of range. Try again.'
            ELSE
              DO I=NS1,NS2
                IFSCAN_POLI(I)=.TRUE.
              END DO
            END IF
          END DO
          K=0
          DO I=1,NSCAN
            IF(IFSCAN_UTIL(I).AND.IFSCAN_POLI(I))THEN
              K=K+1
              XPOL(K)=X(I)
              YPOL(K)=XOFF_FFT(I)
            END IF
          END DO
          NPT=K
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(4)
            CALL PGPOINT(NPT,XPOL,YPOL,17)
          END DO
          WRITE(*,110)'>>> Final no. of points to be fitted: ',NPT
        END IF
C
        WRITE(*,100)'Polynomial degree '
        NTERMS=READILIM('2',0,19)
        NTERMS=NTERMS+1
        WRITE(*,100)'Fitting polynomial...'
        CALL POLFIT(XPOL,YPOL,YPOL,NPT,NTERMS,0,AA,CHISQR)
        WRITE(*,101)'...OK!'
        DO I=1,NSCAN
          XX=REAL(I)
          YY=AA(NTERMS)
          DO K=NTERMS-1,1,-1
            YY=YY*XX+AA(K)
          END DO
          XOFF_POL(I)=YY
        END DO
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM)) CALL PGSCI(5)
          CALL PGLINE(NSCAN,X,XOFF_POL)
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
        END DO
        WRITE(*,100)'Is this fit OK (y/n) '
        CFITOK(1:1)=READC('y','yn')
        IF(CFITOK.EQ.'n') GOTO 50
C------------------------------------------------------------------------------
60      WRITE(*,*)
        WRITE(*,101)'(1) Iterate with new mean spectrum '//
     +   '(distortion corrected)'
        WRITE(*,101)'(2) Divide initial frame (using distortion) '//
     +   'by mean spectrum and save'
        WRITE(*,101)'(3) Shift individual spectrum and save whole frame'
        WRITE(*,101)'(0) QUIT'
        WRITE(*,100)'Option (1/2/3/0)'
        COPC(1:1)=READC('@','0123')
C..............................................................................
        IF(COPC.EQ.'0')THEN
          WRITE(*,100)'Do you really want to quit (y/n) '
          CSURE(1:1)=READC('n','yn')
          IF(CSURE.EQ.'y')THEN
            CALL PGEND
            STOP
          END IF
          GOTO 60
C..............................................................................
        ELSEIF(COPC.EQ.'1')THEN
          GOTO 70
C..............................................................................
        ELSEIF(COPC.EQ.'2')THEN
          GOTO 70
C..............................................................................
        ELSEIF(COPC.EQ.'3')THEN
          GOTO 90
C..............................................................................
        ELSE
          GOTO 60
        END IF
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C recalculamos un nuevo espectro promedio, pudiendo usar ahora un intervalo
C en SCANS mucho mayor dado que antes de sumar los espectros y calcular
C el promedio, corregimos estos del offset calculado anteriormente
70      WRITE(*,*)
        WRITE(*,101)'* Define scan region to calculate NEW '//
     +   'mean spectrum (distortion corrected):'
        LLOOP=.TRUE.
        DO I=1,NSCAN
          IFSCAN_MEAN(I)=.FALSE.
        END DO
        DO WHILE(LLOOP)
          WRITE(*,100)'Scan region (0,0=EXIT) '
          CALL READ2I('0,0',NS1,NS2)
          IF((NS1.EQ.0).AND.(NS2.EQ.0))THEN
            NS0=0
            DO I=1,NSCAN
              IF(IFSCAN_MEAN(I)) NS0=NS0+1
            END DO
            IF(NS0.LT.1)THEN
              WRITE(*,101)'WARNING: no. of scans = 0. Try again.'
            ELSE
              LLOOP=.FALSE.
            END IF
          ELSEIF((NS1.LT.1).OR.(NS2.GT.NSCAN).OR.(NS1.GT.NS2))THEN
            WRITE(*,101)'WARNING: numbers out of range. Try again.'
          ELSE
            DO I=NS1,NS2
              IFSCAN_MEAN(I)=.TRUE.
            END DO
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              IF(LCOLOR(ITERM)) CALL PGSCI(3)
              CALL PGMOVE(REAL(NS1),YMIN)
              CALL PGDRAW(REAL(NS1),YMAX)
              CALL PGMOVE(REAL(NS2),YMIN)
              CALL PGDRAW(REAL(NS2),YMAX)
              CALL PGRECT(REAL(NS1),REAL(NS2),YMIN,YMIN+(YMAX-YMIN)*.01)
              CALL PGRECT(REAL(NS1),REAL(NS2),YMAX,YMAX-(YMAX-YMIN)*.01)
              IF(LCOLOR(ITERM)) CALL PGSCI(1)
            END DO
          END IF
        END DO
C calculamos el espectro promedio
        DO J=1,NCHAN
          MEANSP(J)=0.
        END DO
        NS0=0
        WRITE(*,100)'Obtaining new mean spectrum...:#####'
        DO I=1,NSCAN
          IF(IFSCAN_MEAN(I))THEN
            WRITE(*,'(A,I5,$)')'\\b\\b\\b\\b\\b',I
            NS0=NS0+1
            DO J=1,NCHAN
              SPP1(J)=A(J,I)
            END DO
            CALL CHREBIN(XOFF_POL(I),NCHAN,SPP1,SPP2)
            DO J=1,NCHAN
              MEANSP(J)=MEANSP(J)+SPP2(J)
            END DO
          END IF
        END DO
        WRITE(*,*)
        DO J=1,NCHAN
          MEANSP(J)=MEANSP(J)/REAL(NS0)
        END DO
        IF(COPC.EQ.'2') GOTO 80
        GOTO 40
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
80      WRITE(*,*)
        WRITE(*,100)'Dividing initial frame by mean spectrum...:#####'
        DO I=1,NSCAN
          WRITE(*,'(A,I5,$)')'\\b\\b\\b\\b\\b',I
          CALL CHREBIN(-XOFF_POL(I),NCHAN,MEANSP,SPP2)
          DO J=1,NCHAN
            IF(SPP2(J).GT.0.0)THEN
              B(J,I)=A(J,I)/SPP2(J)
            ELSE
              B(J,I)=1.
            END IF
          END DO
        END DO
        WRITE(*,*)
        WRITE(*,*)
        WRITE(*,100)'Output file name'
        OUTFILE=OUTFILEX(30,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        WRITE(*,100)'Saving file...'
        DO I=1,NSCAN
          WRITE(30) (B(J,I),J=1,NCHAN)
        END DO
        CLOSE(30)
        WRITE(*,101)'...OK! File saved and closed.'
        GOTO 60
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
90      WRITE(*,*)
        WRITE(*,100)'Shifting individual spectrum...:#####'
        DO I=1,NSCAN
          WRITE(*,'(A,I5,$)')'\\b\\b\\b\\b\\b',I
          DO J=1,NCHAN
            SPP1(J)=A(J,I)
          END DO
          CALL CHREBIN(XOFF_POL(I),NCHAN,SPP1,SPP2)
          DO J=1,NCHAN
            B(J,I)=SPP2(J)
          END DO
        END DO
        WRITE(*,*)
        WRITE(*,100)'Output file name'
        OUTFILE=OUTFILEX(30,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        WRITE(*,100)'Saving file...'
        DO I=1,NSCAN
          WRITE(30) (B(J,I),J=1,NCHAN)
        END DO
        CLOSE(30)
        WRITE(*,101)'...OK! File saved and closed.'
        GOTO 60
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
110     FORMAT(A,I6)
        END
C
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
