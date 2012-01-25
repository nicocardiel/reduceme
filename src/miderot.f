C------------------------------------------------------------------------------
C Version 07-September-2007                                     file: miderot.f
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
C Program: miderot
C Classification: measurement
C Description: Determines the rotation curve of a galaxy by using cross 
C correlation.
C
Comment
C------------------------------------------------------------------------------
C
        PROGRAM MIDEROT
C
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER READILIM
C
        INTEGER NMAX
        PARAMETER (NMAX=8192)      !Si se cambia, cambiar tambien en subrutinas
        INTEGER NSIMULMAX    !numero maximo de simulaciones para hallar errores
        PARAMETER (NSIMULMAX=10000)
        REAL RAIZ2
        PARAMETER (RAIZ2=1.41421356)
        REAL PI2
        PARAMETER (PI2=6.28318531)
        REAL C
        PARAMETER (C=299792.46)
C
        INTEGER I,J,L
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        INTEGER NS0,NS1,NS2,NSUM
        INTEGER NC1,NC2
        INTEGER NCHANEFF,NSIZE
        INTEGER K1,K2,K3,K4,NSIDE
        INTEGER NSEED,NSIMUL,NSIMULT
        REAL A(NCMAX,NSMAX),ERR(NCMAX,NSMAX)
        REAL S0(NMAX),S(NMAX),FNOR0,FNOR
        REAL X(NMAX)
        REAL KFILTER(NMAX),COSBELL(NCMAX)
        REAL XCORR(NMAX),FCORR(NMAX)
        REAL XOFF_FFT,YOFF_FFT
        REAL XOFF_SIM(NSIMULMAX)
        REAL XOFF_FIN1(NSMAX),XOFF_FIN2(NSMAX),XOFF_ERR(NSMAX)
        REAL FL
        REAL RANRED,XRAN1,XRAN2
        REAL XMIN,XMAX,YMIN,YMAX,YMIN0,YMAX0,DX,DY
        REAL XMINC,XMAXC,YMINC,YMAXC
        REAL FMEAN,FSIGMA
        REAL WLMEAN
        CHARACTER*1 CERR,CSHOW
        CHARACTER*50 CDUMMY
        CHARACTER*75 INFILE,ERRFILE
        LOGICAL LCOLOR(MAX_ID_RED)
        LOGICAL IFSCAN(NSMAX),IFCHAN(NCMAX)
        LOGICAL LSHOW
C
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
C------------------------------------------------------------------------------
        OUTFILEX=OUTFILEX
C
        IF(NCMAX.GT.NMAX) STOP 'FATAL ERROR: NCMAX.GT.NMAX'
C
        THISPROGRAM='miderot'
        CALL WELCOME('07-September-2007')
C abrimos salida grafica
        CALL PIDEGTER(NTERM,IDN,LCOLOR)
C------------------------------------------------------------------------------
        NSEED=-1
        DO J=1,NMAX
          X(J)=REAL(J)
        END DO
C------------------------------------------------------------------------------
        WRITE(*,100) 'Work with error images (y/n) '
        CERR(1:1)=READC('y','yn')
C
        WRITE(*,100) 'Image to be measured'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        IF(NSCAN.LE.1)THEN
          CLOSE(20)
          STOP 'FATAL ERROR: NSCAN.LE.1'
        END IF
        DO I=1,NSCAN
          READ(20) (A(J,I),J=1,NCHAN)
        END DO
        CLOSE(20)
C
        IF(CERR.EQ.'y')THEN
          WRITE(*,100) 'Input error file name '
          CALL GUESSEF(INFILE,ERRFILE)
          ERRFILE=INFILEX(21,ERRFILE,NSCAN,NCHAN,STWV,DISP,21,.TRUE.) !match
          DO I=1,NSCAN
            READ(21) (ERR(J,I),J=1,NCHAN)
          END DO
          CLOSE(21)
        END IF
C------------------------------------------------------------------------------
        WRITE(*,100) 'Central scan'
        NS0=READILIM('@',1,NSCAN)
C
        DO I=1,NSCAN
          IFSCAN(I)=.FALSE.
        END DO
10      WRITE(*,100) 'Enter scan region to be measured (0,0=EXIT) '
        CALL READ2I('0,0',NS1,NS2)
        IF((NS1.EQ.0).AND.(NS2.EQ.0)) GOTO 12
        IF((NS1.LT.1).OR.(NS2.GT.NSCAN).OR.(NS1.GT.NS2))THEN
          WRITE(*,101) 'ERROR: Invalid entry. Try again.'
        ELSE
          DO I=NS1,NS2
            IFSCAN(I)=.TRUE.
          END DO
        END IF
        GOTO 10
12      NSUM=0
        DO I=1,NSCAN
          IF(IFSCAN(I)) NSUM=NSUM+1
        END DO
        IF(NSUM.EQ.0)THEN
          WRITE(*,101) 'ERROR: number of scans = 0. Try again.'
          GOTO 10
        END IF
C------------------------------------------------------------------------------
        XMIN=1.
        XMAX=REAL(NCHAN)
        DX=XMAX-XMIN
        XMIN=XMIN-DX/20.
        XMAX=XMAX+DX/20.
C
        DO J=1,NCHAN
          S0(J)=A(J,NS0)
        END DO
        CALL FINDMM(NCHAN,S0,YMIN0,YMAX0)
C
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGSCF(2)
          CALL PGENV(XMIN,XMAX,YMIN0,YMAX0,0,0)
          CALL PGLABEL('channel','No. of counts',INFILE)
          CALL PGIDEN_RED
          CALL PGSCI(3)
          CALL PGBIN(NCHAN,X,S0,.TRUE.)
          CALL PGSCI(1)
        END DO
C------------------------------------------------------------------------------
        WRITE(*,101) ' '
        WRITE(*,101) '* Cross-correlation parameters: '
        WRITE(CDUMMY,'(A2,I10)') '1,',NCHAN
        CALL RMBLANK(CDUMMY,CDUMMY,L)
20      WRITE(*,100) '1st and last channel to be employed '
        CALL READ2I(CDUMMY(1:L),NC1,NC2)
        IF((NC1.LT.1).OR.(NC2.GT.NCHAN).OR.(NC1.GT.NC2))THEN
          WRITE(*,101) 'ERROR: invalid numbers. Try again.'
          GOTO 20
        END IF
        NCHANEFF=NC2-NC1+1
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGSCI(7)
          CALL PGSLS(2)
          CALL PGMOVE(REAL(NC1),YMIN0)
          CALL PGDRAW(REAL(NC1),YMAX0)
          CALL PGMOVE(REAL(NC2),YMIN0)
          CALL PGDRAW(REAL(NC2),YMAX0)
          CALL PGSLS(1)
          CALL PGSCI(1)
        END DO
        WLMEAN=STWV+((REAL(NC1+NC2))/2.-1.)*DISP
        DO J=1,NCHAN
          IFCHAN(J)=.FALSE.
        END DO
        DO J=NC1,NC2
          IFCHAN(J)=.TRUE.
        END DO
        CALL FFT2POWER(NCHANEFF,NSIZE)
        K1=10
        K2=20
        K3=50
        K4=100
        CALL FFTKFILTER(NSIZE/2,KFILTER,K1,K2,K3,K4)
        FL=0.10
        CALL FFTCOSBELL(NCHANEFF,COSBELL,FL)
        WRITE(*,100) 'No. of points at each side of maximum to fit '
        WRITE(*,100) 'polynomial '
        NSIDE=READILIM('3',1,NSIZE/2)
        WRITE(*,100) 'Show plots (y/n) '
        CSHOW(1:1)=READC('y','yn')
        LSHOW=(CSHOW.EQ.'y')
C
        IF(CERR.EQ.'y')THEN
          WRITE(*,100) 'No. of simulations to compute errors '
          NSIMULT=READILIM('100',5,NSIMULMAX)
          NSIMULT=NSIMULT+1        !la simulacion numero 1 no va a simular nada
        ELSE
          NSIMULT=1
        END IF
C------------------------------------------------------------------------------
        DO I=1,NSCAN
          IF(IFSCAN(I))THEN
ccc            IF(I.NE.NS0)THEN
C..............................................................................
              DO NSIMUL=1,NSIMULT
                IF(NSIMUL.EQ.1)THEN                      !esto NO es simulacion
                  FNOR0=0.
                  FNOR=0.
                  DO J=1,NCHAN
                    S0(J)=A(J,NS0)
                    S(J)=A(J,I)
                    IF(IFCHAN(J))THEN
                      FNOR0=FNOR0+S0(J)
                      FNOR=FNOR+S(J)
                    END IF
                  END DO
                  FNOR0=FNOR0/FNOR
                  DO J=1,NCHAN
                    S(J)=S(J)*FNOR0
                  END DO
                ELSE                                     !esto SI es simulacion
                  FNOR0=0.
                  FNOR=0.
                  DO J=1,NCHAN
                    XRAN1=RANRED(NSEED)
                    XRAN2=RANRED(NSEED)
                    S0(J)=A(J,NS0)+RAIZ2*ERR(J,NS0)*
     +               SQRT(-ALOG(1.-XRAN1))*COS(PI2*XRAN2)
                    XRAN1=RANRED(NSEED)
                    XRAN2=RANRED(NSEED)
                    S(J)=A(J,I)+RAIZ2*ERR(J,I)*
     +               SQRT(-ALOG(1.-XRAN1))*COS(PI2*XRAN2)
                    IF(IFCHAN(J))THEN
                      FNOR0=FNOR0+S0(J)
                      FNOR=FNOR+S(J)
                    END IF
                  END DO
                  FNOR0=FNOR0/FNOR
                  DO J=1,NCHAN
                    S(J)=S(J)*FNOR0
                  END DO
                END IF
                IF(LSHOW)THEN
                  CALL FINDMM(NCHAN,S0,YMIN0,YMAX0)
                  CALL FINDMM(NCHAN,S,YMIN,YMAX)
                  IF(YMIN0.LT.YMIN) YMIN=YMIN0
                  IF(YMAX0.GT.YMAX) YMAX=YMAX0
                  DY=YMAX-YMIN
                  YMIN=YMIN-DY/20.
                  YMAX=YMAX+DY/20.
                  DO ITERM=NTERM,1,-1
                    CALL PGSLCT(IDN(ITERM))
                    CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,0)
                    CALL PGLABEL('channel','No. of counts',INFILE)
                    CALL PGIDEN_RED
                    CALL PGSCI(3)
                    CALL PGBIN(NCHAN,X,S0,.TRUE.)
                    CALL PGSCI(2)
                    CALL PGBIN(NCHAN,X,S,.TRUE.)
                    CALL PGSCI(7)
                    CALL PGSLS(2)
                    CALL PGMOVE(REAL(NC1),YMIN)
                    CALL PGDRAW(REAL(NC1),YMAX)
                    CALL PGMOVE(REAL(NC2),YMIN)
                    CALL PGDRAW(REAL(NC2),YMAX)
                    CALL PGSLS(1)
                    CALL PGSCI(1)
                  END DO
                END IF
                CALL STOPPLOT(LSHOW)
                CALL FFTPREP(NSIZE,S0,NC1,NC2,.TRUE.,COSBELL,.TRUE.,
     +           KFILTER,LSHOW,INFILE)
                CALL STOPPLOT(LSHOW)
                CALL FFTPREP(NSIZE,S,NC1,NC2,.TRUE.,COSBELL,.TRUE.,
     +           KFILTER,LSHOW,INFILE)
                CALL STOPPLOT(LSHOW)
                CALL FFTCORREL(NSIZE,S0,S,XCORR,FCORR)
                IF(LSHOW)THEN
                  DO ITERM=NTERM,1,-1
                    CALL PGSLCT(IDN(ITERM))
                    CALL PGPAGE
                    CALL PGIDEN_RED
                    CALL AUTOPLOT(NSIZE,XCORR,FCORR,1,NSIZE,
     +               'X-offset','correlation value',
     +               'Correlation of the two data sets',
     +               .TRUE.,XMINC,XMAXC,YMINC,YMAXC,0.05,
     +               0,.FALSE.,'BCNTS','BCNTS',
     +               101,5,
     +               0.,1.,0.5,1.0)
                  END DO
                END IF
                CALL FFTCORRZOOM(XCORR,FCORR,NSIZE,NSIZE-NCHANEFF,
     +           INFILE,INFILE,NSIDE,LSHOW,XOFF_FFT,YOFF_FFT)
                CALL STOPPLOT(LSHOW)
                XOFF_SIM(NSIMUL)=XOFF_FFT
              END DO
C..............................................................................
              IF(NSIMUL.GT.1)THEN
                FMEAN=0.             !recordar de NSIMUL=1 no es una simulacion
                DO NSIMUL=2,NSIMULT
                  FMEAN=FMEAN+XOFF_SIM(NSIMUL)
                END DO
                FMEAN=FMEAN/REAL(NSIMULT-1)
                FSIGMA=0.
                DO NSIMUL=2,NSIMULT
                  FSIGMA=FSIGMA+(XOFF_SIM(NSIMUL)-FMEAN)**2
                END DO
                FSIGMA=SQRT(FSIGMA/REAL(NSIMULT-2))
                XOFF_FIN1(I)=XOFF_SIM(1)       !valor sin simular (sin errores)
                XOFF_FIN2(I)=FMEAN                       !media en simulaciones
                XOFF_ERR(I)=FSIGMA                       !rms en valor anterior
                WRITE(*,100) 'NSCAN, Xoff, <Xoff> and rms: '
                WRITE(*,*) I,XOFF_SIM(1),FMEAN,FSIGMA
              END IF
C..............................................................................
ccc            END IF
          END IF
        END DO
C------------------------------------------------------------------------------
        IF(NSIMULT.GT.1)THEN
          DO I=1,NSCAN
            IF(IFSCAN(I)) XMAX=REAL(I)
          END DO
          DO I=NSCAN,1,-1
            IF(IFSCAN(I)) XMIN=REAL(I)
          END DO
          DX=XMAX-XMIN
          XMIN=XMIN-DX/20.
          XMAX=XMAX+DX/20.
          DO I=1,NSCAN                                     !pasamos a Angstroms
            XOFF_FIN1(I)=-XOFF_FIN1(I)*DISP               !y cambiamos de signo
            XOFF_FIN2(I)=-XOFF_FIN2(I)*DISP               !y cambiamos de signo
            XOFF_ERR(I)=XOFF_ERR(I)*DISP
          END DO
          YMIN=XOFF_FIN1(I)
          YMAX=YMIN
          DO I=1,NSCAN
            IF(IFSCAN(I))THEN
              IF(XOFF_FIN1(I).LT.YMIN) YMIN=XOFF_FIN1(I)
              IF(XOFF_FIN2(I)-XOFF_ERR(I).LT.YMIN) 
     +         YMIN=XOFF_FIN2(I)-XOFF_ERR(I)
              IF(XOFF_FIN1(I).GT.YMAX) YMAX=XOFF_FIN1(I)
              IF(XOFF_FIN2(I)+XOFF_ERR(I).GT.YMAX) 
     +         YMAX=XOFF_FIN2(I)+XOFF_ERR(I)
            END IF
          END DO
          DY=YMAX-YMIN
          YMIN=YMIN-DY/20.
          YMAX=YMAX+DY/20.
          CALL PGSCH(1.5)
          CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,-2)
          CALL PGBOX('BCTSN',0.0,0,'BTSN',0.0,0)
          CALL PGLABEL('scan','offset (\\A)',INFILE)
          DO I=1,NSCAN
            IF(IFSCAN(I))THEN
              CALL PGSCI(2)
              CALL PGPOINT(1,REAL(I),XOFF_FIN1(I),17)
              CALL PGSCI(3)
              CALL PGPOINT(1,REAL(I),XOFF_FIN2(I),13)
              CALL PGERRY(1,REAL(I),XOFF_FIN2(I)-XOFF_ERR(I),
     +         XOFF_FIN2(I)+XOFF_ERR(I),1.)
              CALL PGSCI(1)
            END IF
          END DO
          YMIN=YMIN/WLMEAN*C                    !ponemos escala derecha en km/s
          YMAX=YMAX/WLMEAN*C
          CALL PGWINDOW(XMIN,XMAX,YMIN,YMAX)
          CALL PGBOX(' ',0.0,0,'CTSM',0.0,0)
        END IF
        WRITE(*,*)
        DO I=1,NSCAN
          IF(IFSCAN(I))THEN
            WRITE(*,*)I,XOFF_FIN2(I)/WLMEAN*C,XOFF_ERR(I)/WLMEAN*C
          END IF
        END DO
C------------------------------------------------------------------------------
        CALL PGEND
C
        STOP
C
100     FORMAT(A,$)
101     FORMAT(A)
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
