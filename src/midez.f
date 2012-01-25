C------------------------------------------------------------------------------
C Version 29-June-1998                                            file: midez.f
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
C Program: midez
C Classification: measurement
C Description: Determines redshifts by using cross correlation.
C
Comment
C
C Calcula la velocidad radial de un espectro a partir de otros espectros
C iniciales (templates) que pueden ensancharse a diferentes velocidades
C
C------------------------------------------------------------------------------
C
        PROGRAM MIDEZ
C
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER TRUELEN
        INTEGER READILIM
        REAL READF
C
        INTEGER NMAX
        PARAMETER (NMAX=8192)      !Si se cambia, cambiar tambien en subrutinas
        INTEGER NPLOTMAX
        PARAMETER (NPLOTMAX=101) !numero maximo de puntos en el intervalo de Rv
C
        INTEGER I,II,J,K
        INTEGER L,LL
        INTEGER NTEMPT                          !numero de templates a utilizar
        INTEGER NCOLOR,NCHAN2
        INTEGER NSIGMAS    !diferentes velocidades de ensanchamiento utilizadas
        INTEGER NC1,NC2,NCHANEFF
        INTEGER K1,K2,K3,K4
        INTEGER NSIZE,NSIDE
        INTEGER NPLOT,NP
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        REAL FACTOR
        REAL TEMPL(NCMAX,NSMAX)
        REAL SP(NCMAX),X(NCMAX)
        REAL S(NMAX),SS(NMAX)
        REAL XCORR(NMAX),FCORR(NMAX)
        REAL VEL(NSMAX)                          !velocidades de ensanchamiento
        REAL SIGMAVEC(NCMAX)
        REAL COSBELL(NCMAX),KFILTER(NMAX)
        REAL RADVEL0,DRADVEL,RADVEL1,RADVEL2
        REAL XMIN,XMAX,YMIN,YMAX,DX,DY
        REAL MEAN,SIGMA
        REAL FL,XOFF_FFT,YOFF_FFT
        REAL XPLOT(NPLOTMAX),YPLOT(NPLOTMAX)
        REAL RADVELFIN(NSMAX)       !velocidad radial deducida de cada template
        DOUBLE PRECISION DMEAN
        CHARACTER*1 CBROAD,COK,CSHOW
        CHARACTER*1 COUT,CXLIM
        CHARACTER*75 CDUMMY
        CHARACTER*75 INFILE(NSMAX),INFILE0,OBJECT0,OUTFILE
        LOGICAL LCOLOR(MAX_ID_RED),LSHOW
C
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
        NSIGMAS=0 !avoid compilation warning
C
        THISPROGRAM='midez'
        CALL WELCOME('29-June-1998')
C abrimos salida grafica
        CALL PIDEGTER(NTERM,IDN,LCOLOR)
C------------------------------------------------------------------------------
        WRITE(*,100)'Spectrum to be measured'
        INFILE0=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        OBJECT0(1:75)=OBJECT
        IF(NSCAN.GT.1)THEN
          CLOSE(20)
          STOP 'FATAL ERROR: NSCAN.GT.1'
        END IF
        READ(20) (SP(J),J=1,NCHAN)
        CLOSE(20)
C------------------------------------------------------------------------------
C Velocidad radial inicial
        WRITE(*,100)'Initial radial velocity  (km/sec)'
        RADVEL0=READF('@')
        WRITE(*,100)'Radial velocity interval (km/sec)'
        DRADVEL=READF('@')
        RADVEL1=RADVEL0-DRADVEL
        RADVEL2=RADVEL0+DRADVEL
C------------------------------------------------------------------------------
        DO J=1,NCHAN
          X(J)=REAL(J)
        END DO
C------------------------------------------------------------------------------
C numero de estrellas template a utilizar
        WRITE(*,100)'No. of template spectra to be used '
        NTEMPT=READILIM('1',1,NSMAX)
C
        DO I=1,NTEMPT
          WRITE(*,'(A,I3.3,A,$)')'Template #',I,' file name '//
     +     '(Rv= 0 km/sec)'
          INFILE(I)=INFILEX(20,'@',NSCAN,NCHAN2,STWV,DISP,1,.FALSE.) !match
          IF(NCHAN2.NE.NCHAN)THEN
            WRITE(*,101)'FATAL ERROR: NCHAN is different in last image!'
            CLOSE(20)
            STOP
          END IF
          IF(NSCAN.NE.1)THEN
            WRITE(*,101)'FATAL ERROR: NSCAN.NE.1'
            CLOSE(20)
            STOP
          END IF
          READ(20) (TEMPL(J,I),J=1,NCHAN)
          CLOSE(20)
        END DO
C------------------------------------------------------------------------------
C normalizamos todos los espectros template
        DO I=1,NTEMPT
          DMEAN=0.D0
          DO J=1,NCHAN
            DMEAN=DMEAN+DBLE(TEMPL(J,I))
          END DO
          DMEAN=DMEAN/DBLE(NCHAN)
          FACTOR=1./REAL(DMEAN)
          DO J=1,NCHAN
            TEMPL(J,I)=TEMPL(J,I)*FACTOR
          END DO
        END DO
C normalizamos espectro problema
        DMEAN=0.
        DO J=1,NCHAN
          DMEAN=DMEAN+DBLE(SP(J))
        END DO
        DMEAN=DMEAN/DBLE(NCHAN)
        FACTOR=1./REAL(DMEAN)
        DO J=1,NCHAN
          SP(J)=SP(J)*FACTOR
        END DO
C------------------------------------------------------------------------------
C representamos graficamente los espectros (sin ensanchar)
10      XMIN=1.
        XMAX=REAL(NCHAN)
        YMIN=TEMPL(1,1)
        YMAX=YMIN
        DO I=1,NTEMPT
          DO J=1,NCHAN
            IF(YMIN.GT.TEMPL(J,I)) YMIN=TEMPL(J,I)
            IF(YMAX.LT.TEMPL(J,I)) YMAX=TEMPL(J,I)
          END DO
        END DO
        CALL RVREBIN(-RADVEL0,NCHAN,SP,SS,STWV,DISP)
        DO J=1,NCHAN
          IF(YMIN.GT.SS(J)) YMIN=SS(J)
          IF(YMAX.LT.SS(J)) YMAX=SS(J)
        END DO
        DX=XMAX-XMIN
        DY=YMAX-YMIN
        XMIN=XMIN-DX/50.
        XMAX=XMAX+DX/50.
        YMIN=YMIN-DY/50.
        YMAX=YMAX+DY/50.
C
12      DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGENV(XMIN,XMAX,YMIN-.5,YMAX+.5,0,0)
          CALL PGIDEN_RED
        END DO
        NCOLOR=1
C
        DO I=1,NTEMPT
          DO J=1,NCHAN
            S(J)=TEMPL(J,I)
          END DO
          NCOLOR=NCOLOR+1
          IF(NCOLOR.GT.14) NCOLOR=1
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(NCOLOR)
            CALL PGBIN(NCHAN,X,S,.TRUE.)
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END DO
        END DO
        CALL RVREBIN(-RADVEL0,NCHAN,SP,SS,STWV,DISP)
        DO J=1,NCHAN
          SS(J)=SS(J)+0.5
        END DO
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGBIN(NCHAN,X,SS,.TRUE.)
        END DO
        CALL RVREBIN(-RADVEL1,NCHAN,SP,SS,STWV,DISP)
        DO J=1,NCHAN
          SS(J)=SS(J)-0.5
        END DO
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGBIN(NCHAN,X,SS,.TRUE.)
        END DO
        CALL RVREBIN(-RADVEL2,NCHAN,SP,SS,STWV,DISP)
        DO J=1,NCHAN
          SS(J)=SS(J)-0.5
        END DO
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGBIN(NCHAN,X,SS,.TRUE.)
          CALL PGLABEL('nchan','normalized no. of counts',CHAR(32))
          CALL PGMTEXT('T',1.5,0.,0.,'File: '//INFILE0)
          CALL PGMTEXT('T',1.5,1.,1.,'Object: '//
     +     OBJECT0(1:TRUELEN(OBJECT0)))
        END DO
C------------------------------------------------------------------------------
        WRITE(*,*)
        WRITE(*,100)'Change X-limits (y/n) '
        CXLIM(1:1)=READC('n','yn')
        IF(CXLIM.EQ.'y')THEN
          WRITE(CDUMMY,*)NCHAN
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          WRITE(*,101)'* Valid range is 1 to '//CDUMMY(1:L)//': '
          WRITE(CDUMMY,*)XMIN
          WRITE(*,100)'Xmin '
          XMIN=READF(CDUMMY)
          WRITE(CDUMMY,*)XMAX
          WRITE(*,100)'Xmax '
          XMAX=READF(CDUMMY)
          GOTO 12
        END IF
C------------------------------------------------------------------------------
        WRITE(*,100)'Is the radial velocity interval ok (y/n) '
        COK(1:1)=READC('y','yn')
        IF(COK.EQ.'n')THEN
          WRITE(*,100)'Radial velocity interval (km/sec)'
          DRADVEL=READF('@')
          RADVEL1=RADVEL0-DRADVEL
          RADVEL2=RADVEL0+DRADVEL
          GOTO 10
        END IF
C------------------------------------------------------------------------------
C velocidades de ensanchamiento de los espectros
        IF(NTEMPT.EQ.1)THEN
          WRITE(*,100)'Are you broadening this spectrum (y/n) '
        ELSE
          WRITE(*,100)'Are you broadening these spectra (y/n) '
        END IF
        CBROAD(1:1)=READC('n','yn')
C
        IF(CBROAD.EQ.'y')THEN
          WRITE(*,100)'No. of velocities '
          NSIGMAS=READILIM('@',1,(NSMAX-NTEMPT)/NTEMPT)
          DO K=1,NSIGMAS
            WRITE(*,'(A,I3,A,$)')'Velocity #',K,' (km/sec)'
            VEL(K)=READF('@')
          END DO
        END IF
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
C determinamos potencia de 2 superior a NCHANEFF
        CALL FFT2POWER(NCHANEFF,NSIZE)
C..............................................................................
C Definimos frecuencias de filtrado
        K1=10
        K2=20
        K3=50
        K4=100
        CALL FFTKFILTER(NSIZE/2,KFILTER,K1,K2,K3,K4)
C..............................................................................
C Definimos la campana de coseno
        FL=0.10
        CALL FFTCOSBELL(NCHANEFF,COSBELL,FL)
C..............................................................................
C Numero de puntos a cada lado del maximo para encontrar el offset
        WRITE(*,*)
        WRITE(*,100)'No. of points at each side of maximum to fit '//
     +   'polynomial '
        NSIDE=READILIM('3',1,NSIZE/2)
C..............................................................................
C Numero de intervalos en velocidad radial
        WRITE(*,100)'No. of intervals in radial velocity (odd) '
        NPLOT=READILIM('7',3,NPLOTMAX)
C..............................................................................
C Si queremos, podemos dibujar
        WRITE(*,100)'Show plots (y/n) '
        CSHOW(1:1)=READC('n','yn')
        LSHOW=(CSHOW.EQ.'y')
C..............................................................................
        WRITE(*,100)'Save results into an output text file (y/n) '
        COUT(1:1)=READC('n','yn')
        IF(COUT.EQ.'y')THEN
          WRITE(*,100)'Output file name'
          OUTFILE=OUTFILEX(30,'@',0,0,0.,0.,3,.FALSE.)
          WRITE(30,150)
          WRITE(30,101)'File..: '//INFILE0(1:TRUELEN(INFILE0))
          WRITE(30,101)'Object: '//OBJECT0(1:TRUELEN(OBJECT0))
          WRITE(30,150)
          WRITE(30,150)
          WRITE(30,110)'First channel: ',NC1
          WRITE(30,110)'Last  channel: ',NC2
          WRITE(30,150)
        END IF
C------------------------------------------------------------------------------
C si se ha solicitado, ensanchamos los espectros template
        IF(CBROAD.EQ.'y')THEN
          WRITE(*,101)'Broadening spectra (wait)...'
          DO K=1,NSIGMAS
            WRITE(CDUMMY,*)VEL(K)
            CALL RMBLANK(CDUMMY,CDUMMY,LL)
            DO I=1,NTEMPT
              II=NTEMPT+(K-1)*NTEMPT+I
              DO J=1,NCHAN
                S(J)=TEMPL(J,I)
              END DO
              DO J=1,NCHAN
                SIGMAVEC(J)=VEL(K)
              END DO
              CALL BROADEN(S,SS,NCHAN,STWV,DISP,SIGMAVEC,.FALSE.)
              DO J=1,NCHAN
                TEMPL(J,II)=SS(J)
              END DO
              L=TRUELEN(INFILE(I))
              INFILE(II)=INFILE(I)(1:L)//' [broad.: '//
     +         CDUMMY(1:LL)//' km/sec]'
              WRITE(*,101)INFILE(II)//' OK!'
            END DO
          END DO
          WRITE(*,*)
        ELSE
          NSIGMAS=0
        END IF
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
        DRADVEL=(RADVEL2-RADVEL1)/REAL(NPLOT-1)
C Correlacion cruzada con las templates sin ensanchar
        DO I=1,NTEMPT
          WRITE(*,*)
          WRITE(*,101)'* Sigma: 0.0'
          WRITE(*,100)'* Template: '
          L=TRUELEN(INFILE(I))
          WRITE(*,101)INFILE(I)(1:L)
          DO NP=1,NPLOT
            RADVEL0=REAL(NP-1)*DRADVEL+RADVEL1
            CALL RVREBIN(-RADVEL0,NCHAN,SP,SS,STWV,DISP)
            CALL FFTPREP(NSIZE,SS,NC1,NC2,.TRUE.,COSBELL,.TRUE.,
     +       KFILTER,LSHOW,INFILE0)
            CALL STOPPLOT(LSHOW)
            DO J=1,NCHAN
              S(J)=TEMPL(J,I)
            END DO
            CALL FFTPREP(NSIZE,S,NC1,NC2,.TRUE.,COSBELL,.TRUE.,
     +       KFILTER,LSHOW,INFILE(I))
            CALL STOPPLOT(LSHOW)
            CALL FFTCORREL(NSIZE,S,SS,XCORR,FCORR)
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
              END DO
            END IF
            WRITE(*,100)'* Radial velocity: '
            WRITE(*,*)RADVEL0
            CALL FFTCORRZOOM(XCORR,FCORR,NSIZE,NSIZE-NCHANEFF,
     +       INFILE0,INFILE(I),NSIDE,LSHOW,XOFF_FFT,YOFF_FFT)
            CALL STOPPLOT(LSHOW)
            XPLOT(NP)=XOFF_FFT
            YPLOT(NP)=RADVEL0
          END DO
          CALL FITZ(NPLOT,XPLOT,YPLOT,RADVEL1,RADVEL2,RADVEL0,
     +     INFILE0,INFILE(I))
          RADVELFIN(I)=RADVEL0
        END DO
C Correlacion cruzada con las templates ensanchadas
        IF(NSIGMAS.GT.0)THEN
          DO K=1,NSIGMAS
            DO I=1,NTEMPT
              WRITE(*,*)
              WRITE(*,100)'* Sigma: '
              WRITE(*,*) VEL(K)
              WRITE(*,100)'* Template: '
              II=NTEMPT+(K-1)*NTEMPT+I
              DO NP=1,NPLOT
                RADVEL0=REAL(NP-1)*DRADVEL+RADVEL1
                CALL RVREBIN(-RADVEL0,NCHAN,SP,SS,STWV,DISP)
                CALL FFTPREP(NSIZE,SS,NC1,NC2,.TRUE.,COSBELL,.TRUE.,
     +           KFILTER,LSHOW,INFILE0)
                CALL STOPPLOT(LSHOW)
                DO J=1,NCHAN
                  S(J)=TEMPL(J,II)
                END DO
                CALL FFTPREP(NSIZE,S,NC1,NC2,.TRUE.,COSBELL,.TRUE.,
     +           KFILTER,LSHOW,INFILE(I))
                CALL STOPPLOT(LSHOW)
                CALL FFTCORREL(NSIZE,S,SS,XCORR,FCORR)
                IF(LSHOW)THEN
                  DO ITERM=NTERM,1,-1
                    CALL PGSLCT(IDN(ITERM))
                    CALL PGPAGE
                    CALL PGIDEN_RED
                    CALL AUTOPLOT(NSIZE,XCORR,FCORR,1,NSIZE,
     +               'X-offset','correlation value',
     +                'Correlation of the two data sets',
     +               .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +               0,.FALSE.,'BCNTS','BCNTS',
     +               101,5,
     +               0.,1.,0.5,1.0)
                  END DO
                END IF
                WRITE(*,100)'* Radial velocity: '
                WRITE(*,*)RADVEL0
                CALL FFTCORRZOOM(XCORR,FCORR,NSIZE,NSIZE-NCHANEFF,
     +           INFILE0,INFILE(I),NSIDE,LSHOW,XOFF_FFT,YOFF_FFT)
                CALL STOPPLOT(LSHOW)
                XPLOT(NP)=XOFF_FFT
                YPLOT(NP)=RADVEL0
              END DO
              CALL FITZ(NPLOT,XPLOT,YPLOT,RADVEL1,RADVEL2,RADVEL0,
     +         INFILE0,INFILE(II))
              RADVELFIN(II)=RADVEL0
            END DO
          END DO
        END IF
C------------------------------------------------------------------------------
        CALL PGEND
C------------------------------------------------------------------------------
C Valores medios para templates sin ensanchar
        WRITE(*,150)
        IF(COUT.EQ.'y') WRITE(30,150)
        MEAN=0.
        DO I=1,NTEMPT
          MEAN=MEAN+RADVELFIN(I)
          WRITE(*,'(A60,$)')INFILE(I)
          WRITE(*,*) RADVELFIN(I)
          IF(COUT.EQ.'y')THEN
            WRITE(30,'(A60,$)')INFILE(I)
            WRITE(30,*) RADVELFIN(I)
          END IF
        END DO
        MEAN=MEAN/REAL(NTEMPT)
        WRITE(*,100)'>>> Mean: '
        WRITE(*,*) MEAN
        IF(COUT.EQ.'y')THEN
          WRITE(30,100)'>>> Mean: '
          WRITE(30,*) MEAN
        END IF
        IF(NTEMPT.GT.1)THEN
          SIGMA=0.
          DO I=1,NTEMPT
            SIGMA=SIGMA+(RADVELFIN(I)-MEAN)*(RADVELFIN(I)-MEAN)
          END DO
          SIGMA=SQRT(SIGMA/REAL(NTEMPT-1))
          WRITE(*,100)'>>> Std.: '
          WRITE(*,*) SIGMA
          IF(COUT.EQ.'y')THEN
            WRITE(30,100)'>>> Std.: '
            WRITE(30,*) SIGMA
          END IF
        END IF
        WRITE(*,150)
        IF(COUT.EQ.'y') WRITE(30,150)
C Valores medios para templates ensanchadas
        IF(NSIGMAS.GT.0)THEN
          DO K=1,NSIGMAS
            MEAN=0.
            DO I=1,NTEMPT
              II=NTEMPT+(K-1)*NTEMPT+I
              MEAN=MEAN+RADVELFIN(II)
              WRITE(*,'(A60,$)')INFILE(II)
              WRITE(*,*) RADVELFIN(II)
              IF(COUT.EQ.'y')THEN
                WRITE(30,'(A60,$)')INFILE(II)
                WRITE(30,*) RADVELFIN(II)
              END IF
            END DO
            MEAN=MEAN/REAL(NTEMPT)
            WRITE(*,100)'>>> Mean: '
            WRITE(*,*) MEAN
            IF(COUT.EQ.'y')THEN
              WRITE(30,100)'>>> Mean: '
              WRITE(30,*) MEAN
            END IF
            IF(NTEMPT.GT.1)THEN
              SIGMA=0.
              DO I=1,NTEMPT
                II=NTEMPT+(K-1)*NTEMPT+I
                SIGMA=SIGMA+(RADVELFIN(II)-MEAN)*(RADVELFIN(II)-MEAN)
              END DO
              SIGMA=SQRT(SIGMA/REAL(NTEMPT-1))
              WRITE(*,100)'>>> Std.: '
              WRITE(*,*) SIGMA
              IF(COUT.EQ.'y')THEN
                WRITE(30,100)'>>> Std.: '
                WRITE(30,*) SIGMA
              END IF
            END IF
            WRITE(*,150)
            IF(COUT.EQ.'y') WRITE(30,150)
          END DO
        END IF
        IF(COUT.EQ.'y') CLOSE(30)
C------------------------------------------------------------------------------
        STOP
100     FORMAT(A,$)
101     FORMAT(A)
110     FORMAT(A,I6)
150     FORMAT(79('-'))
        END
C
C******************************************************************************
C Calcula el valor de velocidad radial al cual corresponde un offset de cero.
C Este valor final se transmite a traves de la variable output RADVEL0
C
        SUBROUTINE FITZ(NPLOT,XPLOT,YPLOT,RADVEL1,RADVEL2,RADVEL0,
     +   INFILE0,INFILE)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INTEGER TRUELEN
C
        INTEGER NPLOT
        REAL XPLOT(NPLOT)
        REAL YPLOT(NPLOT)
        REAL RADVEL1,RADVEL2,RADVEL0
        CHARACTER*(*) INFILE0
        CHARACTER*(*) INFILE
C
        INTEGER I,L
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        REAL XMIN,XMAX,YMIN,YMAX,DX,DY
        REAL B(3),CHISQR                                           !para POLFIT
        REAL X(1000),Y(1000)
        LOGICAL LCOLOR(MAX_ID_RED)
C
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        CALL FINDMM(NPLOT,XPLOT,XMIN,XMAX)
        YMIN=RADVEL1
        YMAX=RADVEL2
        DX=XMAX-XMIN
        DY=YMAX-YMIN
        XMIN=XMIN-DX/50.
        XMAX=XMAX+DX/50.
        YMIN=YMIN-DY/50.
        YMAX=YMAX+DY/50.
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,0)
          CALL PGPOINT(NPLOT,XPLOT,YPLOT,17)
          CALL PGLABEL('Offset from cross correlation (channels)',
     +     'Radial velocity (km/sec)',' ')
          L=TRUELEN(INFILE0)
          CALL PGMTEXT('T',1.5,0.,0.,INFILE0(1:L))
          L=TRUELEN(INFILE)
          CALL PGMTEXT('T',1.5,1.,1.,INFILE(1:L))
        END DO
C ajustamos un polinomio de segundo grado
        CALL POLFIT(XPLOT,YPLOT,YPLOT,NPLOT,3,0,B,CHISQR)
        DO I=1,1000
          X(I)=XMIN+REAL(I-1)/999.*(XMAX-XMIN)
          Y(I)=B(1)+B(2)*X(I)+B(3)*X(I)*X(I)
        END DO
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM)) CALL PGSCI(2)
          CALL PGLINE(1000,X,Y)
          IF(LCOLOR(ITERM)) CALL PGSCI(3)
          CALL PGSLS(2)
          CALL PGMOVE(XMIN,B(1))
          CALL PGDRAW(XMAX,B(1))
          CALL PGMOVE(0.,YMIN)
          CALL PGDRAW(0.,YMAX)
          CALL PGSLS(1)
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
        END DO
C
        RADVEL0=B(1)
        WRITE(*,100)'>>> Radial velocity (km/sec): '
        WRITE(*,*) RADVEL0
C
100     FORMAT(A,$)
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
