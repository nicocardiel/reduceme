C------------------------------------------------------------------------------
C Version 3-July-1998                                           file: corrfft.f
C @ ncl & fjg
C------------------------------------------------------------------------------
Comment
C
C Program: corrfft
C Classification: arithmetic & manipulations
C Description: Cross-correlates spectra using FFT.
C
Comment
C
C Realiza correlacion cruzada de dos senhales.
C NOTA: el programa ha sido ampliado para utilizar varios espectros REFERENCE
C y determinar asi una estimacion del error del offset.
C
C Este programa toma dos espectros y realiza una correlacion de ambos, para
C asi determinar el posible offset entre ellos. Los espectros son leidos,
C se indica la region de ambos espectros que va a utilizarse, y luego
C se tiene varias opciones. En primer lugar, es posible filtrar los espectros
C utilizando la transformada de Fourier. Para ello se utiliza un filtro de 
C frecuencias lineal definido en 4 puntos (K1,K2,K3 y K4), de modo que hasta
C K1 el filtro es cero, entre K1 y K2 aumenta linealmente de cero a uno,
C entre K2 y K3 es uno, entre K3 y K4 disminuye linealmente de uno a cero,
C y a partir de K4 es cero. Posteriormente (con los espectros filtrados
C o sin filtrar), se puede ajustar un polinomio para eliminar la forma
C a gran escala (si hemos filtrado con un valor adecuado de K1 el ajuste
C al polinomio no es necesario). Las frecuencias por debajo de K1 eliminan
C la forma a gran escala, mientras que las frecuencias por encima de K3
C deben eliminar el ruido.
C
C La rutina que calcula FFT (y por tanto tambien la que calcula correlacion
C cruzada) necesita que el numero de puntos en el espectro sea un potencia de
C 2**N (con N un numero entero). El programa determina el valor de N mas proximo
C al tamanho de nuestros espectros, y rellena con ceros hasta completar el valor
C 2**N. Es importante que el tamanho util de nuestros espectros no sea
C exactamente 2**N para que asi el programa pueda insertar dichos ceros (zero
C padding), dado que los ceros evitan posibles problemas de borde al realizar la
C correlacion. De hecho, el numero de ceros introducidos determina el offset
C maximo para el cual es fiable la correlacion realizada.  Asi, si el tamanho
C util de los espectros a correlacionar es exactamente 2**N, el programa sugiere
C utilizar 2**(N+1) para introducir tantos ceros como dicho tamanho util (es el
C caso mas exagerado).
C
C Por ultimo, el programa permite corregir el segundo espectro del offset
C encontrado, pudiendose salvar dicho espectro.
C
        PROGRAM CORRFFT
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER READILIM
        REAL READF
C
        INTEGER NMAX
        PARAMETER (NMAX=8192)      !Si se cambia, cambiar tambien en subrutinas
        INTEGER MAXNTERMS
        PARAMETER (MAXNTERMS=20)
        INTEGER NTEMPMAX
        PARAMETER (NTEMPMAX=20)              !No. maximo de espectros REFERENCE
        REAL PI
        PARAMETER (PI=3.141593)
C
        INTEGER I,J,L
        INTEGER NTEMP,NTEMPTOT
        INTEGER N,NZPAD
        INTEGER NS0,NC1,NC2
        INTEGER NSCAN2,NCHAN2
        INTEGER NDEG         !numero de terminos en el ajuste de los polinomios
        INTEGER K1,K2,K3,K4
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        REAL STWV2,DISP2
        REAL X(NMAX)
        REAL SORIG1(NCMAX),SORIG2(NCMAX),SORIGFIN(NCMAX)
        REAL SP1(NCMAX),SP2(NCMAX)
        REAL SPR1(NMAX),SPR2(NMAX),SPI1(NMAX),SPI2(NMAX)
        REAL SP1R(NMAX),SP2R(NMAX)
        REAL KFILTER(NMAX)
        REAL A1(MAXNTERMS),A2(MAXNTERMS)
        REAL POL1(NCMAX),POL2(NCMAX)
        REAL CHISQR,POL
        REAL COSBELL(NCMAX),FL
        REAL DATA1(NMAX),DATA2(NMAX)
        REAL FCORR(NMAX),XCORR(NMAX)
        REAL MODULO1(NMAX),MODULO2(NMAX),FASE1(NMAX),FASE2(NMAX)
        REAL POWERLOG1(NMAX),POWERLOG2(NMAX)
        REAL XMIN,XMAX,YMIN,YMAX
        REAL X0,Y0
        REAL X0TAB(NTEMPMAX),X0MEAN,X0SIGMA
        DOUBLE PRECISION MEAN
        CHARACTER*1 CCOLOR,CFILT,CRFILT,CSHIFT,CSHIFTERR,COUT
        CHARACTER*1 CCONT,CREP,CSHOW
        CHARACTER*120 CDUMMY
        CHARACTER*75 INFILE1,INFILE2,OUTFILE,FILEXOFF,ERRFILE
        LOGICAL LPL1
        LOGICAL LCOLOR(MAX_ID_RED)
C
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
C------------------------------------------------------------------------------
        THISPROGRAM='corrfft'
        CALL WELCOME('3-July-1998')
C valores por defecto para filtrar
        K1=10
        K2=20
        K3=50
        K4=100
C
        DO J=1,NMAX
          X(J)=REAL(J)
        END DO
C abrimos salida grafica
        CALL PIDEGTER(NTERM,IDN,LCOLOR)
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGSCH(1.2)
        END DO
        IF(LCOLOR(1))THEN
          CCOLOR='y'
        ELSE
          CCOLOR='n'
        END IF
        WRITE(*,100)'Show all plots (y/n) '
        CSHOW(1:1)=READC('y','yn')
C------------------------------------------------------------------------------
        WRITE(*,100)'Number of REFERENCE spectra '
        NTEMPTOT=READILIM('1',1,NTEMPMAX)
        NTEMP=0
C entramos espectros
10      NTEMP=NTEMP+1
        WRITE(*,100)'REFERENCE spectrum file name'
        INFILE1=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        IF(NSCAN.GT.1)THEN
          WRITE(*,101)'WARNING: this file contains more than 1 '//
     +     'spectrum.'
          WRITE(*,100)'Scan number to be employed '
          NS0=READILIM('@',1,NSCAN)
        ELSE
          NS0=1
        END IF
        DO I=1,NS0
          READ(20) (SORIG1(J),J=1,NCHAN)
        END DO
        CLOSE(20)
        DO J=1,NCHAN
          SP1(J)=SORIG1(J)
        END DO
        MEAN=0.D0
        DO J=1,NCHAN
          MEAN=MEAN+SP1(J)
        END DO
        MEAN=MEAN/DBLE(NCHAN)
        IF(MEAN.LT.0.0) MEAN=-MEAN
        DO J=1,NCHAN
          SP1(J)=SP1(J)/REAL(MEAN)
        END DO
C
        IF(NTEMP.GT.1) GOTO 12          !si no es la primera template, saltamos
C
        WRITE(*,100)'PROBLEM   spectrum file name'
        INFILE2=INFILEX(20,'@',NSCAN2,NCHAN2,STWV2,DISP2,1,.FALSE.)
        IF((NSCAN.NE.NSCAN2).OR.(NCHAN.NE.NCHAN2).OR.(STWV.NE.STWV2)
     +   .OR.(DISP.NE.DISP2))THEN
          WRITE(*,101)'ERROR: header information in PROBLEM is '//
     +     'different.'
          WRITE(*,100)'Do you want to continue anyway (y/n) '
          CCONT(1:1)=READC('n','yn')
          IF(CCONT.EQ.'n')THEN
            CLOSE(20)
            STOP
          END IF
        END IF
        IF(NSCAN.GT.1)THEN
          WRITE(*,101)'WARNING: this file contains more than '//
     +     '1 spectrum.'
          WRITE(*,100)'Scan number to be employed '
          NS0=READILIM('@',1,NSCAN)
        ELSE
          NS0=1
        END IF
        DO I=1,NS0
          READ(20) (SORIG2(J),J=1,NCHAN)
        END DO
        CLOSE(20)
C inicializamos limites a utilizar para la primera estrella REFERENCE
        NC1=1
        NC2=NCHAN
C
12      DO J=1,NCHAN
          SP2(J)=SORIG2(J)
        END DO
        MEAN=0.D0
        DO J=1,NCHAN
          MEAN=MEAN+SP2(J)
        END DO
        IF(MEAN.LT.0.0) MEAN=-MEAN
        MEAN=MEAN/DBLE(NCHAN)
        DO J=1,NCHAN
          SP2(J)=SP2(J)/REAL(MEAN)
        END DO
C si tenemos color, dibujamos los dos espectros en la misma grafica
        IF(CSHOW.EQ.'y')THEN
          IF(CCOLOR.EQ.'y')THEN
            YMAX=SP1(1)
            DO J=2,NCHAN
              IF(SP1(J).GT.YMAX) YMAX=SP1(J)
            END DO
            LPL1=.TRUE.
            DO J=1,NCHAN
              IF(SP2(J).GT.YMAX)THEN
                LPL1=.FALSE.
                YMAX=SP2(J)
              END IF
            END DO
            IF(LPL1)THEN
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                CALL PGIDEN_RED
                CALL AUTOPLOT(NCHAN,X,SP1,1,NCHAN,
     +           'channel','counts',' ',
     +           .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +           0,.FALSE.,'BCNTS','BCNTS',
     +           101,3,
     +           0.,1.,0.0,1.0)
                CALL AUTOPLOT(NCHAN,X,SP2,1,NCHAN,
     +           'channel','counts',' ',
     +           .FALSE.,XMIN,XMAX,YMIN,YMAX,0.00,
     +           0,.FALSE.,' ',' ',
     +           101,2,
     +           0.,1.,0.0,1.0)
                CALL PGSCI(3)
                CALL PGMTEXT('T',2.,0.,0.,'REFERENCE: '//INFILE1)
                CALL PGSCI(2)
                CALL PGMTEXT('T',2.,1.,1.,'PROBLEM: '//INFILE2)
                CALL PGSCI(1)
              END DO
            ELSE
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                CALL PGIDEN_RED
                CALL AUTOPLOT(NCHAN,X,SP2,1,NCHAN,
     +           'channel','counts',' ',
     +           .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +           0,.FALSE.,'BCNTS','BCNTS',
     +           101,3,
     +           0.,1.,0.0,1.0)
                CALL AUTOPLOT(NCHAN,X,SP1,1,NCHAN,
     +           'channel','counts',' ',
     +           .FALSE.,XMIN,XMAX,YMIN,YMAX,0.00,
     +           0,.FALSE.,' ',' ',
     +           101,2,
     +           0.,1.,0.0,1.0)
                CALL PGSCI(3)
                CALL PGMTEXT('T',2.,0.,0.,'REFERENCE: '//INFILE1)
                CALL PGSCI(2)
                CALL PGMTEXT('T',2.,1.,1.,'PROBLEM: '//INFILE2)
                CALL PGSCI(1)
              END DO
            END IF
          ELSE
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              CALL PGIDEN_RED
              CALL AUTOPLOT(NCHAN,X,SP1,1,NCHAN,
     +         'channel','counts','REFERENCE: '//INFILE1,
     +         .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +         0,.FALSE.,'BCNTS','BCNTS',
     +         101,3,
     +         0.,1.,0.5,1.0)
              CALL AUTOPLOT(NCHAN,X,SP2,1,NCHAN,
     +         'channel','counts','PROBLEM: '//INFILE2,
     +         .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +         0,.FALSE.,'BCNTS','BCNTS',
     +         101,2,
     +         0.,1.,0.0,0.5)
            END DO
          END IF
        END IF
C------------------------------------------------------------------------------
C definimos region a utilizar de los espectros
20      WRITE(CDUMMY,'(I10,A1,I10)')NC1,',',NC2
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,100)'1st & last channel to be employed '
        CALL READ2I(CDUMMY(1:L),NC1,NC2)
        IF((NC1.LT.1).OR.(NC2.GT.NCHAN).OR.(NC1.GT.NC2))THEN
          WRITE(*,101)'ERROR: invalid numbers. Try again.'
          GOTO 20
        END IF
        DO J=NC1,NC2
          SP1(J-NC1+1)=SP1(J)
          SP2(J-NC1+1)=SP2(J)
        END DO
        NCHAN=NC2-NC1+1
        IF(CSHOW.EQ.'y')THEN
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGPAGE
            CALL PGIDEN_RED
            CALL AUTOPLOT(NCHAN,X,SP1,1,NCHAN,
     +       'channel','counts','REFERENCE: '//INFILE1,
     +       .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +       0,.FALSE.,'BCNTS','BCNTS',
     +       101,3,
     +       0.,1.,0.5,1.0)
            CALL AUTOPLOT(NCHAN,X,SP2,1,NCHAN,
     +       'channel','counts','PROBLEM: '//INFILE2,
     +       .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +       0,.FALSE.,'BCNTS','BCNTS',
     +       101,2,
     +       0.,1.,0.0,0.5)
          END DO
        END IF
C------------------------------------------------------------------------------
C calculamos una campana de coseno para suavizar los bordes (ver Brault &
C White, A&A, 13, 169)
        FL=0.10
        CALL FFTCOSBELL(NCHAN,COSBELL,FL)
C dibujamos la campana de coseno
        IF(CSHOW.EQ.'y')THEN
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGPAGE
            CALL AUTOPLOT(NCHAN,X,COSBELL,1,NCHAN,
     +                    'channel','function value','cosine bell',
     +                    .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +                    0,.FALSE.,'BCNTS','BCNTS',
     +                    101,5,
     +                    0.,1.,0.0,1.0)
            CALL PGIDEN_RED
          END DO
          WRITE(*,100)'Press <CR> to continue...'
          READ(*,*)
        END IF
C aplicamos la campana de coseno
        WRITE(*,100)'Applying cosine bell...'
        DO J=1,NCHAN
          SP1(J)=SP1(J)*COSBELL(J)
          SP2(J)=SP2(J)*COSBELL(J)
        END DO
        IF(CSHOW.EQ.'y')THEN
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGPAGE
            CALL PGIDEN_RED
            CALL AUTOPLOT(NCHAN,X,SP1,1,NCHAN,
     +       'channel','counts','REFERENCE: '//INFILE1,
     +       .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +       0,.FALSE.,'BCNTS','BCNTS',
     +       101,3,
     +       0.,1.,0.5,1.0)
            CALL AUTOPLOT(NCHAN,X,SP2,1,NCHAN,
     +       'channel','counts','PROBLEM: '//INFILE2,
     +       .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +       0,.FALSE.,'BCNTS','BCNTS',
     +       101,2,
     +       0.,1.,0.0,0.5)
          END DO
        END IF
        WRITE(*,*)
C------------------------------------------------------------------------------
        CALL FFT2POWER(NCHAN,N)
        NZPAD=N-NCHAN
C------------------------------------------------------------------------------
        WRITE(*,*)
        WRITE(*,100)'Are you filtering the spectra (y/n) '
        CFILT(1:1)=READC('y','yn')
        IF(CFILT.EQ.'y')THEN
28        DO J=1,NCHAN
            SPR1(J)=SP1(J)
            SPR2(J)=SP2(J)
            SPI1(J)=0.
            SPI2(J)=0.
          END DO
          IF(N.GT.NCHAN)THEN
            WRITE(*,100)'Applying zero padding...'
            DO J=1,NZPAD
              SPR1(NCHAN+J)=0.
              SPR2(NCHAN+J)=0.
              SPI1(NCHAN+J)=0.
              SPI2(NCHAN+J)=0.
            END DO
            WRITE(*,*)
            IF(CSHOW.EQ.'y')THEN
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                CALL PGPAGE
                CALL PGIDEN_RED
                CALL AUTOPLOT(N,X,SPR1,1,N,
     +           'channel','counts','REFERENCE: '//INFILE1,
     +           .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +           0,.FALSE.,'BCNTS','BCNTS',
     +           101,3,
     +           0.,1.,0.5,1.0)
                CALL AUTOPLOT(N,X,SPR2,1,N,
     +           'channel','counts','PROBLEM: '//INFILE2,
     +           .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +           0,.FALSE.,'BCNTS','BCNTS',
     +           101,2,
     +           0.,1.,0.0,0.5)
              END DO
              WRITE(*,100)'Press <CR>...'
              READ(*,*)
            END IF
          END IF
          CALL CFFT(N,SPR1,SPI1,1)
          CALL CFFT(N,SPR2,SPI2,1)
          DO J=1,N
            MODULO1(J)=SQRT(SPR1(J)*SPR1(J)+SPI1(J)*SPI1(J))
            FASE1(J)=ATAN2(SPI1(J),SPR1(J))
            MODULO2(J)=SQRT(SPR2(J)*SPR2(J)+SPI2(J)*SPI2(J))
            FASE2(J)=ATAN2(SPI2(J),SPR2(J))
          END DO
          IF(CSHOW.EQ.'y')THEN
            DO J=1,N
              POWERLOG1(J)=ALOG10(MODULO1(J))
              POWERLOG2(J)=ALOG10(MODULO2(J))
            END DO
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              CALL PGPAGE
              CALL PGIDEN_RED
              CALL AUTOPLOT(N,X,POWERLOG1,1,N,
     +         'K axis: discrete transform domain',
     +          'Power Spectrum: log\\d10\\u(P)',
     +          'frequency (in units of the Nyquist frequency)',
     +         .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +         0,.FALSE.,'BNTS','BNTS',
     +         101,3,
     +         0.,1.,0.5,1.0)
              CALL PGMTEXT('T',2.0,1.0,1.0,'REFERENCE: '//INFILE1)
              CALL PGBOX(' ',0.0,0,'C',0.0,0)
              CALL PGWINDOW(-0.10,2.10,YMIN,YMAX)
              CALL PGBOX('CMTS',0.0,0,' ',0.0,0)
              CALL PGWINDOW(XMIN,XMAX,YMIN,YMAX)
              CALL AUTOPLOT(N,X,POWERLOG2,1,N,
     +         'K axis: discrete transform domain',
     +          'Power Spectrum: log\\d10\\u(P)',
     +          'frequency (in units of the Nyquist frequency)',
     +         .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +         0,.FALSE.,'BNTS','BNTS',
     +         101,2,
     +         0.,1.,0.0,0.5)
              CALL PGMTEXT('T',2.0,1.0,1.0,'PROBLEM: '//INFILE2)
              CALL PGBOX(' ',0.0,0,'C',0.0,0)
              CALL PGWINDOW(-0.10,2.10,YMIN,YMAX)
              CALL PGBOX('CMTS',0.0,0,' ',0.0,0)
              CALL PGWINDOW(XMIN,XMAX,YMIN,YMAX)
            END DO
          END IF
          CALL FFTKFILTER(N,KFILTER,K1,K2,K3,K4)
          IF(CSHOW.EQ.'y')THEN
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              CALL PGSCI(5)
              CALL AUTOPLOT(N,X,KFILTER,1,N,
     +         ' ',' ',' ',
     +         .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +         0,.FALSE.,' ','CMTS',
     +         101,5,
     +         0.,1.,0.5,1.0)
              CALL PGMTEXT('R',3.0,0.5,0.5,'Filter')
              CALL AUTOPLOT(N,X,KFILTER,1,N,
     +         ' ',' ',' ',
     +         .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +         0,.FALSE.,' ','CMTS',
     +         101,5,
     +         0.,1.,0.0,0.5)
              CALL PGMTEXT('R',3.0,0.5,0.5,'Filter')
              CALL PGSCI(1)
            END DO
            WRITE(*,100)'Press <CR>...'
            READ(*,*)
          END IF
C aplicamos el filtro
          DO J=1,N
            MODULO1(J)=MODULO1(J)*KFILTER(J)
            MODULO2(J)=MODULO2(J)*KFILTER(J)
          END DO
C aplicamos los valores filtrados
          DO J=1,N
            SPR1(J)=MODULO1(J)*COS(FASE1(J))
            SPI1(J)=MODULO1(J)*SIN(FASE1(J))
            SPR2(J)=MODULO2(J)*COS(FASE2(J))
            SPI2(J)=MODULO2(J)*SIN(FASE2(J))
          END DO
C FFT inversa
          CALL CFFT(N,SPR1,SPI1,-1)
          CALL CFFT(N,SPR2,SPI2,-1)
          DO J=1,N
            SP1R(J)=SPR1(J)
            SP2R(J)=SPR2(J)
          END DO
          IF(CSHOW.EQ.'y')THEN
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              CALL PGPAGE
              CALL PGIDEN_RED
              CALL AUTOPLOT(NCHAN,X,SP1R,1,NCHAN,
     +         'channel','counts','REFERENCE: '//INFILE1,
     +         .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +         0,.FALSE.,'BCNTS','BCNTS',
     +         101,3,
     +         0.,1.,0.5,1.0)
              CALL AUTOPLOT(NCHAN,X,SP2R,1,NCHAN,
     +         'channel','counts','PROBLEM: '//INFILE2,
     +         .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +         0,.FALSE.,'BCNTS','BCNTS',
     +         101,2,
     +         0.,1.,0.0,0.5)
            END DO
          END IF
          WRITE(*,100)'Repeat filtering (y/n) '
          CRFILT(1:1)=READC('n','yn')
          IF(CRFILT.EQ.'y') GOTO 28
        ELSE
          DO J=1,NCHAN
            SP1R(J)=SP1(J)
            SP2R(J)=SP2(J)
          END DO
        END IF
C------------------------------------------------------------------------------
        DO J=1,NCHAN
          SP1(J)=SP1R(J)
          SP2(J)=SP2R(J)
        END DO
C ajustamos la forma del continuo con un polinomio
30      WRITE(*,*)
        WRITE(*,100)'Polynomial degree (to remove continuum, '//
     +   '-1=NO FIT) '
        NDEG=READILIM('-1',-1,19)
        IF(NDEG.EQ.-1)THEN
          GOTO 40
        END IF
        NDEG=NDEG+1
        IF(NDEG.GT.MAXNTERMS)THEN
          WRITE(*,100)'ERROR: invalid entry. Try again.'
          GOTO 30
        END IF
        CALL POLFIT(X,SP1R,SP1R,NCHAN,NDEG,0,A1,CHISQR)
        DO J=1,NCHAN
          POL=A1(NDEG)
          DO I=NDEG-1,1,-1
            POL=POL*X(J)+A1(I)
          END DO
          POL1(J)=POL
        END DO
        CALL POLFIT(X,SP2R,SP2R,NCHAN,NDEG,0,A2,CHISQR)
        DO J=1,NCHAN
          POL=A2(NDEG)
          DO I=NDEG-1,1,-1
            POL=POL*X(J)+A2(I)
          END DO
          POL2(J)=POL
        END DO
        IF(CSHOW.EQ.'y')THEN
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGPAGE
            CALL PGIDEN_RED
            CALL AUTOPLOT(NCHAN,X,SP1R,1,NCHAN,
     +       'channel','counts','REFERENCE: '//INFILE1,
     +       .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +       0,.FALSE.,'BCNTS','BCNTS',
     +       101,3,
     +       0.,1.,0.5,1.0)
            CALL AUTOPLOT(NCHAN,X,POL1,1,NCHAN,
     +       ' ',' ',' ',
     +       .FALSE.,XMIN,XMAX,YMIN,YMAX,0.00,
     +       0,.FALSE.,' ',' ',
     +       101,4,
     +       0.,1.,0.5,1.0)
            CALL AUTOPLOT(NCHAN,X,SP2R,1,NCHAN,
     +       'channel','counts','PROBLEM: '//INFILE2,
     +       .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +       0,.FALSE.,'BCNTS','BCNTS',
     +       101,2,
     +       0.,1.,0.0,0.5)
            CALL AUTOPLOT(NCHAN,X,POL2,1,NCHAN,
     +       ' ',' ',' ',
     +       .FALSE.,XMIN,XMAX,YMIN,YMAX,0.00,
     +       0,.FALSE.,' ',' ',
     +       101,4,
     +       0.,1.,0.0,0.5)
          END DO
        END IF
        WRITE(*,100)'Repeat polynomial fit '
        CREP(1:1)=READC('n','yn')
        IF(CREP.EQ.'y') GOTO 30
C eliminamos forma del continuo dividiendo por el polinomio
        WRITE(*,100)'Applying polynomial...'
        DO J=1,NCHAN
          SP1(J)=SP1R(J)-POL1(J)
          SP2(J)=SP2R(J)-POL2(J)
        END DO
        WRITE(*,*)
        IF(CSHOW.EQ.'y')THEN
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGPAGE
            CALL PGIDEN_RED
            CALL AUTOPLOT(NCHAN,X,SP1,1,NCHAN,
     +       'channel','counts','REFERENCE: '//INFILE1,
     +       .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +       0,.FALSE.,'BCNTS','BCNTS',
     +       101,3,
     +       0.,1.,0.5,1.0)
            CALL AUTOPLOT(NCHAN,X,SP2,1,NCHAN,
     +       'channel','counts','PROBLEM: '//INFILE2,
     +       .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +       0,.FALSE.,'BCNTS','BCNTS',
     +       101,2,
     +       0.,1.,0.0,0.5)
          END DO
          WRITE(*,100)'Press <CR>...'
          READ(*,*)
        END IF
C------------------------------------------------------------------------------
C aplicamos la campana de coseno otra vez
40      WRITE(*,100)'Applying cosine bell...'
        DO J=1,NCHAN
          SP1(J)=SP1(J)*COSBELL(J)
          SP2(J)=SP2(J)*COSBELL(J)
        END DO
        WRITE(*,*)
        IF(CSHOW.EQ.'y')THEN
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGPAGE
            CALL PGIDEN_RED
            CALL AUTOPLOT(NCHAN,X,SP1,1,NCHAN,
     +       'channel','counts','REFERENCE: '//INFILE1,
     +       .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +       0,.FALSE.,'BCNTS','BCNTS',
     +       101,3,
     +       0.,1.,0.5,1.0)
            CALL AUTOPLOT(NCHAN,X,SP2,1,NCHAN,
     +       'channel','counts','PROBLEM: '//INFILE2,
     +       .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +       0,.FALSE.,'BCNTS','BCNTS',
     +       101,2,
     +       0.,1.,0.0,0.5)
          END DO
          WRITE(*,100)'Press <CR>...'
          READ(*,*)
        END IF
C------------------------------------------------------------------------------
C calculamos matrices de datos a correlacionar
        DO J=1,NCHAN
          DATA1(J)=SP1(J)
          DATA2(J)=SP2(J)
        END DO
        IF(N.GT.NCHAN)THEN
          WRITE(*,100)'Applying zero padding...'
          DO J=1,NZPAD
            DATA1(NCHAN+J)=0.
            DATA2(NCHAN+J)=0.
          END DO
          WRITE(*,*)
        END IF
        IF(CSHOW.EQ.'y')THEN
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGPAGE
            CALL PGIDEN_RED
            CALL AUTOPLOT(N,X,DATA1,1,N,
     +       'channel','counts','REFERENCE: '//INFILE1,
     +       .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +       0,.FALSE.,'BCNTS','BCNTS',
     +       101,3,
     +       0.,1.,0.5,1.0)
            CALL AUTOPLOT(N,X,DATA2,1,N,
     +       'channel','counts','PROBLEM: '//INFILE2,
     +       .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +       0,.FALSE.,'BCNTS','BCNTS',
     +       101,2,
     +       0.,1.,0.0,0.5)
          END DO
          WRITE(*,100)'Press <CR>...'
          READ(*,*)
        END IF
C------------------------------------------------------------------------------
C calculamos correlacion
        CALL FFTCORREL(N,DATA1,DATA2,XCORR,FCORR)
C
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGPAGE
          CALL PGIDEN_RED
          CALL AUTOPLOT(N,XCORR,FCORR,1,N,
     +     'X-offset','correlation value',
     +      'Correlation of the two data sets',
     +     .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +     0,.FALSE.,'BCNTS','BCNTS',
     +     101,5,
     +     0.,1.,0.5,1.0)
          CALL PGMTEXT('T',1.5,0.0,0.0,INFILE1)
          CALL PGMTEXT('T',1.5,1.0,1.0,INFILE2)
        END DO
        CALL FFTCORRZOOM(XCORR,FCORR,N,N-NCHAN,INFILE1,INFILE2,0,
     +   .TRUE.,X0,Y0)
        X0TAB(NTEMP)=X0
        IF(NTEMP.LT.NTEMPTOT) GOTO 10
        IF(NTEMPTOT.GT.1)THEN
          X0MEAN=0.
          DO I=1,NTEMP
            WRITE(*,'(A,I2.2,A,$)')'X0[',I,']= '
            WRITE(*,*)X0TAB(I)
            X0MEAN=X0MEAN+X0TAB(I)
          END DO
          X0MEAN=X0MEAN/REAL(NTEMP)
          IF(NTEMPTOT.EQ.2)THEN
            X0SIGMA=ABS(X0TAB(1)-X0TAB(2))
          ELSE
            X0SIGMA=0.
            DO I=1,NTEMP
              X0SIGMA=X0SIGMA+(X0TAB(I)-X0MEAN)*(X0TAB(I)-X0MEAN)
            END DO
            X0SIGMA=SQRT(X0SIGMA/REAL(NTEMPTOT-1))
          END IF
          WRITE(*,100)'Mean x0 value: '
          WRITE(*,*)X0MEAN
          WRITE(*,100)'Sigma        : '
          WRITE(*,*)X0SIGMA
          X0=X0MEAN
        END IF
C------------------------------------------------------------------------------
C si se solicita, traladamos el espectro PROBLEM para igualarlo al REFERENCE 
        WRITE(*,100)'Shift and save PROBLEM spectrum (y/n) '
        CSHIFT(1:1)=READC('n','yn')
        IF(CSHIFT.EQ.'y')THEN
          WRITE(CDUMMY,*)X0
          WRITE(*,100)'Shift '
          X0=READF(CDUMMY)
          CALL CHREBIN(X0,NCHAN2,SORIG2,SORIGFIN)
          WRITE(*,100)'Output file name'
          OUTFILE=OUTFILEX(30,'@',1,NCHAN2,STWV2,DISP2,1,.FALSE.)
          WRITE(30) (SORIGFIN(J),J=1,NCHAN2)
          CLOSE(30)
          WRITE(*,100)'Are you shifting an error PROBLEM '//
     +     'spectrum (y/n) '
          CSHIFTERR(1:1)=READC('n','yn')
          IF(CSHIFTERR.EQ.'y')THEN
            WRITE(*,100)'Input error file name '
            CALL GUESSEF(INFILE2,ERRFILE)
            ERRFILE=INFILEX(21,ERRFILE,1,NCHAN2,STWV2,DISP2,21,.TRUE.) !..match
            READ(21) (SORIG2(J),J=1,NCHAN2)
            CLOSE(21)
            CALL CHREBIN(X0,NCHAN2,SORIG2,SORIGFIN)
            WRITE(*,100)'Output error file name '
            CALL GUESSEF(OUTFILE,ERRFILE)
            ERRFILE=OUTFILEX(31,ERRFILE,1,NCHAN2,STWV2,DISP2,1,.TRUE.)
            WRITE(31) (SORIGFIN(J),J=1,NCHAN2)
            CLOSE(31)
          END IF
        END IF
C------------------------------------------------------------------------------
        WRITE(*,100)'Write final offset into file (y/n) '
        COUT(1:1)=READC('n','yn')
        IF(COUT.EQ.'y')THEN
          WRITE(*,100)'Output file name'
          FILEXOFF=OUTFILEX(35,'@',0,0,0.,0.,3,.FALSE.)
          WRITE(35,*)X0
          CLOSE(35)
        END IF
C------------------------------------------------------------------------------
        CALL PGEND
        STOP
100     FORMAT(A,$)
101     FORMAT(A)
        END
