C------------------------------------------------------------------------------
C Version 3-July-1998                                          file: tfourier.f
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
C Program: tfourier
C Classification: arithmetic & manipulations
C Description: Filters a spectrum using FFT.
C
Comment
C
C Realizar la transformada de Fourier y permite filtrar frecuencias,
C generando como salida el espectro filtrado
C
        PROGRAM TFOURIER
C
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER READI
        INTEGER READILIM
        REAL READF
C
        INTEGER NMAX
        PARAMETER (NMAX=8192)
C
        INTEGER I,J
        INTEGER NS0
        INTEGER NSHIFT,N
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        INTEGER K1,K2,K3,K4
        REAL SIGNAL0(NCMAX),SIGNALN(NCMAX),SIGNAL(NCMAX)
        REAL MODULO(NMAX),FASE(NMAX)
        REAL COSBELL(NCMAX),FL
        REAL KFILTER(NMAX)
        REAL X(NMAX),YPLOT(NMAX)
        REAL XMIN,XMAX,YMIN,YMAX
        REAL YMIN0,YMAX0,DX,DY
        REAL NEWVALUE
        REAL XR(NMAX),XI(NMAX)
        REAL FACTOR
        DOUBLE PRECISION DMEAN
        CHARACTER*1 CSHIFT,CFILTER,CCL,CSAVE,CNOR
        CHARACTER*75 INFILE,OUTFILE
        LOGICAL LCOLOR(MAX_ID_RED)
C------------------------------------------------------------------------------
        THISPROGRAM='tfourier'
        CALL WELCOME('3-July-1998')
C------------------------------------------------------------------------------
        DO J=1,NMAX
          X(J)=REAL(J)
        END DO
C
        CALL PIDEGTER(NTERM,IDN,LCOLOR)
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGSCH(1.5)
          CALL PGSLW(3)
        END DO
C------------------------------------------------------------------------------
C leemos el espectro de trabajo
        WRITE(*,100)'Input file name'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        IF(NSCAN.GT.1)THEN
          WRITE(*,101)'WARNING: NSCAN.GT.1'
          WRITE(*,100)'Scan no. to be read '
          NS0=READILIM('1',1,NSCAN)
        ELSE
          NS0=1
        END IF
        DO I=1,NS0
          READ(20) (SIGNAL0(J),J=1,NCHAN)
        END DO
        CLOSE(20)
C------------------------------------------------------------------------------
C normalizamos
        WRITE(*,100)'Do you want to normalize this spectrum (y/n) '
        CNOR(1:1)=READC('n','yn')
        IF(CNOR.EQ.'y')THEN
          DMEAN=0.D0
          DO J=1,NCHAN
            DMEAN=DMEAN+SIGNAL0(J)
          END DO
          DMEAN=DMEAN/DBLE(NCHAN)
          FACTOR=1./REAL(DMEAN)
        ELSE
          FACTOR=1.0
        END IF
        DO J=1,NCHAN
          SIGNALN(J)=SIGNAL0(J)*FACTOR
        END DO
C------------------------------------------------------------------------------
C dibujamos el espectro
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(CNOR.EQ.'y')THEN
            CALL AUTOPLOT(NCHAN,X,SIGNALN,1,NCHAN,
     +       'channel','normalized data','normalized original spectrum',
     +       .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +       0,.FALSE.,'BCNTS','BNTS',
     +       101,1,
     +       0.,1.,0.,1.)
          ELSE
            CALL AUTOPLOT(NCHAN,X,SIGNALN,1,NCHAN,
     +       'channel','data','original spectrum',
     +       .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +       0,.FALSE.,'BCNTS','BNTS',
     +       101,1,
     +       0.,1.,0.,1.)
          END IF
          CALL PGBOX(' ',0.0,0,'C',0.0,0)
          CALL PGIDEN_RED
        END DO
C cambiamos los limites del dibujo
        CCL='y'
        DO WHILE(CCL.EQ.'y')
          WRITE(*,100)'Change limits (y/n) '
          CCL(1:1)=READC('n','yn')
          IF(CCL.EQ.'y')THEN
            CALL CHANGEL(XMIN,XMAX,YMIN,YMAX)
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              CALL PGPAGE
              CALL PGIDEN_RED
              IF(CNOR.EQ.'y')THEN
                CALL AUTOPLOT(NCHAN,X,SIGNALN,1,NCHAN,
     +           'channel','normalized data',
     +           'normalized original spectrum',
     +           .FALSE.,XMIN,XMAX,YMIN,YMAX,0.00,
     +           0,.FALSE.,'BCNTS','BNTS',
     +           101,1,
     +           0.,1.,0.,1.)
              ELSE
                CALL AUTOPLOT(NCHAN,X,SIGNALN,1,NCHAN,
     +           'channel','normalized data',
     +           'normalized original spectrum',
     +           .FALSE.,XMIN,XMAX,YMIN,YMAX,0.00,
     +           0,.FALSE.,'BCNTS','BNTS',
     +           101,1,
     +           0.,1.,0.,1.)
              END IF
              CALL PGBOX(' ',0.0,0,'C',0.0,0)
            END DO
          END IF
        END DO
C------------------------------------------------------------------------------
C campana de coseno para bordes
        FL=0.10
        CALL FFTCOSBELL(NCHAN,COSBELL,FL)
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGSCI(5)
          CALL AUTOPLOT(NCHAN,X,COSBELL,1,NCHAN,
     +     ' ',' ',' ',
     +     .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +     0,.FALSE.,' ','CMTSI',
     +     101,5,
     +     0.,1.,0.,1.)
          CALL PGMTEXT('R',3.0,0.5,0.5,'cosine bell')
          CALL PGSCI(1)
        END DO
C------------------------------------------------------------------------------
C multiplicamos datos por campana de coseno
        DO J=1,NCHAN
          SIGNALN(J)=SIGNALN(J)*COSBELL(J)
        END DO
C
        WRITE(*,100)'Press <CR> to continue...'
        READ(*,*)
C
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGPAGE
          CALL PGIDEN_RED
          IF(CNOR.EQ.'y')THEN
            CALL AUTOPLOT(NCHAN,X,SIGNALN,1,NCHAN,
     +       'channel','normalized data',
     +        'original spectrum x cosine bell',
     +       .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +       0,.FALSE.,'BCNTS','BCNTS',
     +       101,1,
     +       0.,1.,0.,1.)
          ELSE
            CALL AUTOPLOT(NCHAN,X,SIGNALN,1,NCHAN,
     +       'channel','data',
     +        'original spectrum x cosine bell',
     +       .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +       0,.FALSE.,'BCNTS','BCNTS',
     +       101,1,
     +       0.,1.,0.,1.)
          END IF
        END DO
C------------------------------------------------------------------------------
C calculamos potencia de 2
        CALL FFT2POWER(NCHAN,N)
C inicializamos matrices para el calculo de la FFT
65      CONTINUE
        DO J=1,NCHAN
          XR(J)=SIGNALN(J)
          XI(J)=0.
        END DO
C
        IF(N.NE.NCHAN)THEN
          DO J=NCHAN+1,N
            XR(J)=0.
            XI(J)=0.
          END DO
        END IF
C calculo de la FFT y del espectro de frecuencias
        CALL CFFT(N,XR,XI,1)
C
        DO J=1,N
          MODULO(J)=SQRT(XR(J)*XR(J)+XI(J)*XI(J))
          FASE(J)=ATAN2(XI(J),XR(J))
        END DO
C
        DO J=1,N
          YPLOT(J)=ALOG10(MODULO(J))
        END DO
C dibujamos el espectro de potencias
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGPAGE
          CALL PGIDEN_RED
          CALL AUTOPLOT(N,X,YPLOT,1,N,
     +     'K axis: discrete transform domain',
     +       'POWER SPECTRUM: Log\\d10\\u(P)',
     +       'frequency (in units of the Nyquist frequency)',
     +     .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +     0,.FALSE.,'BNTS','BNTS',
     +     101,2,
     +     0.,1.,0.,1.)
          CALL PGBOX(' ',0.0,0,'C',0.0,0)
          CALL PGWINDOW(-0.10,2.10,YMIN,YMAX)
          CALL PGBOX('CMTS',0.0,0,' ',0.0,0)
          CALL PGWINDOW(XMIN,XMAX,YMIN,YMAX)
        END DO
C------------------------------------------------------------------------------
C definimos filtrado de frecuencias
        K1=10
        K2=20
        K3=50
        K4=100
        CALL FFTKFILTER(N,KFILTER,K1,K2,K3,K4)
C
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGSCI(3)
          CALL AUTOPLOT(N,X,KFILTER,1,N,
     +     ' ',' ',' ',
     +     .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +     0,.FALSE.,' ','CMTSI',
     +     101,3,
     +     0.,1.,0.,1.)
          CALL PGMTEXT('R',3.0,0.5,0.5,'filter')
          CALL PGSCI(1)
        END DO
C efectuamos el filtrado
        DO J=1,N
          MODULO(J)=MODULO(J)*KFILTER(J)
        END DO
C
        WRITE(*,100)'Press <CR> to continue...'
        READ(*,*)
C------------------------------------------------------------------------------
C calculamos la FFT inversa a partir de las frecuencias filtradas
        DO J=1,N
          XR(J)=MODULO(J)*COS(FASE(J))
          XI(J)=MODULO(J)*SIN(FASE(J))
        END DO
C
        CALL CFFT(N,XR,XI,-1)
C
        DO J=1,NCHAN
          SIGNAL(J)=XR(J)
        END DO
C------------------------------------------------------------------------------
C comparamos espectro inicial con el filtrado
        XMIN=1.
        XMAX=REAL(NCHAN)
        DX=XMAX-XMIN
        XMIN=XMIN-0.05*DX
        XMAX=XMAX+0.05*DX
        CALL FINDMM(NCHAN,SIGNALN,YMIN0,YMAX0)
        CALL FINDMM(NCHAN,SIGNAL,YMIN,YMAX)
        IF(YMIN0.LT.YMIN) YMIN=YMIN0
        IF(YMAX0.GT.YMAX) YMAX=YMAX0
        DY=YMAX-YMIN
        YMIN=YMIN-0.05*DY
        YMAX=YMAX+0.05*DY
C
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGPAGE
          CALL PGIDEN_RED
          CALL AUTOPLOT(NCHAN,X,SIGNALN,1,NCHAN,
     +     'channel','normalized data',
     +      'original vs filtered spectrum',
     +     .FALSE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +     0,.FALSE.,'BCNTS','BCNTS',
     +     101,1,
     +     0.,1.,0.,1.)
          CALL AUTOPLOT(NCHAN,X,SIGNAL,1,NCHAN,
     +     ' ',' ',' ',
     +     .FALSE.,XMIN,XMAX,YMIN,YMAX,0.00,
     +     0,.FALSE.,' ',' ',
     +     101,2,
     +     0.,1.,0.,1.)
        END DO
C------------------------------------------------------------------------------
        WRITE(*,100)'Shift restored data towards the left (y/n) '
        CSHIFT(1:1)=READC('n','yn')
        IF(CSHIFT.EQ.'y')THEN
          WRITE(*,100)'How many pixels'
          NSHIFT=READI('@')
          WRITE(*,100)'New signal value for nonexisting data'
          NEWVALUE=READF('@')
          DO J=1,NCHAN-NSHIFT
            SIGNAL(J)=SIGNAL(J+NSHIFT)
          END DO
          DO J=NCHAN-NSHIFT+1,NCHAN
            SIGNAL(J)=NEWVALUE
          END DO
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(3)
            CALL PGLINE(NCHAN,X,SIGNAL)
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END DO
        END IF
        WRITE(*,100)'Repeat filter (y/n) '
        CFILTER(1:1)=READC('n','yn')
        IF(CFILTER.EQ.'y')GOTO 65
        CALL PGEND
C
        WRITE(*,100)'Save output file (y/n) '
        CSAVE(1:1)=READC('n','yn')
        IF(CSAVE.EQ.'y')THEN
          WRITE(*,100)'Output file name'
          OUTFILE=OUTFILEX(30,'@',1,NCHAN,STWV,DISP,1,.FALSE.)
          DO J=2,NCHAN-1
            SIGNAL(J)=SIGNAL(J)/COSBELL(J)
          END DO
          WRITE(30) (SIGNAL(J),J=1,NCHAN)
          CLOSE(30)
        END IF
C
        STOP
100     FORMAT(A,$)
101     FORMAT(A)
        END
C
C**********************************************************************
C
        SUBROUTINE CHANGEL(XMIN,XMAX,YMIN,YMAX)
        IMPLICIT NONE
        REAL READF
C
        REAL XMIN,XMAX,YMIN,YMAX
        CHARACTER*50 CDUMMY
C----------------------------------------------------------------------
        WRITE(CDUMMY,*)XMIN
        WRITE(*,100)'Xmin '
        XMIN=READF(CDUMMY)
        WRITE(CDUMMY,*)XMAX
        WRITE(*,100)'Xmax '
        XMAX=READF(CDUMMY)
        WRITE(CDUMMY,*)YMIN
        WRITE(*,100)'Ymin '
        YMIN=READF(CDUMMY)
        WRITE(CDUMMY,*)YMAX
        WRITE(*,100)'Ymax '
        YMAX=READF(CDUMMY)
100     FORMAT(A,$)
        END
