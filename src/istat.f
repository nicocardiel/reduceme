C------------------------------------------------------------------------------
C Version 7-December-1996                                         file: istat.f
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
C Program: istat
C Classification: examination & statistics
C Description: Provides some statistics of an image.
C
Comment
C
C Realiza estadistica sobre una imagen (igual que el comando de FIGARO). Si
C se introduce tambien imagen de errores, la estadistica proporciona el valor
C promedio de la relacion senhal/ruido.
C
        PROGRAM ISTAT
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER READI
C
        INTEGER I,J,L,L1,L2
        INTEGER IOPC
        INTEGER NPIX,NPIXT
        INTEGER NC1,NC2,NS1,NS2
        INTEGER NNCMIN,NNCMAX,NNSMIN,NNSMAX
        INTEGER NN
        REAL A(NCMAX,NSMAX),ERR(NCMAX,NSMAX)
        REAL PIXEL(NCMAX*NSMAX)
        REAL FMEDIAN,FMEDIAN1,FMEAN2,FMEANSIGMA
        REAL SMIN,SMAX
        DOUBLE PRECISION MEAN,SIGMA,MEANSN,SIGMASN
        CHARACTER*1 CERR
        CHARACTER*50 CDUMMY
        CHARACTER*75 INFILE,ERRFILE
        LOGICAL IFSCAN(NSMAX),IFCHAN(NCMAX)
C------------------------------------------------------------------------------
        OUTFILEX=OUTFILEX
C evita warnings de compilacion
        NPIXT=0
        MEANSN=0.D0
        SIGMASN=0.D0
C
        THISPROGRAM='istat'
        CALL WELCOME('7-December-1996')
C
        WRITE(*,100)'Work with error images (y/n) '
        CERR(1:1)=READC('n','yn')
C
10      WRITE(*,100)'Input file name'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        DO I=1,NSCAN
          READ(20) (A(J,I),J=1,NCHAN)
        END DO
        CLOSE(20)
        IF(CERR.EQ.'y')THEN
          WRITE(*,100)'Input error file name '
          CALL GUESSEF(INFILE,ERRFILE)
          ERRFILE=INFILEX(22,ERRFILE,NSCAN,NCHAN,STWV,DISP,21,.TRUE.) !match
          DO I=1,NSCAN
            READ(22) (ERR(J,I),J=1,NCHAN)
          END DO
          CLOSE(22)
        END IF
C
20      DO I=1,NSCAN
          IFSCAN(I)=.FALSE.
        END DO
        DO J=1,NCHAN
          IFCHAN(J)=.FALSE.
        END DO
C
        IF(NSCAN.EQ.1)THEN
          IFSCAN(1)=.TRUE.
        ELSE
          WRITE(CDUMMY,'(A,I10)')'1,',NSCAN
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          WRITE(*,101)'Valid region is: '//CDUMMY(1:L)
21        WRITE(*,100)'1st & Last Scan (0,0=EXIT) '
          CALL READ2I('0,0',NS1,NS2)
          IF((NS1.EQ.0).AND.(NS2.EQ.0)) GOTO 22
          IF((NS1.LT.1).OR.(NS2.GT.NSCAN).OR.(NS2.LT.NS1))THEN
            WRITE(*,101)'ERROR: numbers out of range. Try again.'
            GOTO 21
          END IF
          DO I=NS1,NS2
            IFSCAN(I)=.TRUE.
          END DO
          GOTO 21
22        NN=0
          DO I=1,NSCAN
            IF(IFSCAN(I)) NN=NN+1
          END DO
          IF(NN.EQ.0)THEN
            WRITE(*,101)'ERROR: no. of scans to be used = 0!'
            GOTO 21
          END IF
        END IF
C
        IF(NCHAN.EQ.1)THEN
          IFCHAN(1)=.TRUE.
        ELSE
          WRITE(CDUMMY,'(A,I10)')'1,',NCHAN
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          WRITE(*,101)'Valid region is: '//CDUMMY(1:L)
25        WRITE(*,100)'1st & Last Channel (0,0=EXIT) '
          CALL READ2I('0,0',NC1,NC2)
          IF((NC1.EQ.0).AND.(NC2.EQ.0)) GOTO 26
          IF((NC1.LT.1).OR.(NC2.GT.NCHAN).OR.(NC2.LT.NC1))THEN
            WRITE(*,101)'ERROR: numbers out of range. Try again.'
            GOTO 25
          END IF
          DO J=NC1,NC2
            IFCHAN(J)=.TRUE.
          END DO
          GOTO 25
26        NN=0
          DO J=1,NCHAN
            IF(IFCHAN(J)) NN=NN+1
          END DO
          IF(NN.EQ.0)THEN
            WRITE(*,101)'ERROR: no. of channels to be used = 0!'
            GOTO 25
          END IF
        END IF
C
        WRITE(*,100)'Thinking...'
C
        NPIX=0
        DO I=1,NSCAN
          IF(IFSCAN(I))THEN
            DO J=1,NCHAN
              IF(IFCHAN(J))THEN
                NPIX=NPIX+1
                PIXEL(NPIX)=A(J,I)
              END IF
            END DO
          END IF
        END DO
        IF(NPIX.GE.2)THEN
          FMEANSIGMA=FMEAN2(NPIX,PIXEL,3.0)  !media eliminando puntos a 3 sigma
          FMEDIAN=FMEDIAN1(NPIX,PIXEL)                                 !mediana
        ELSE
          FMEANSIGMA=0.
          FMEDIAN=0.
        END IF
C
        DO I=1,NSCAN
          IF(IFSCAN(I))THEN
            DO J=1,NCHAN
              IF(IFCHAN(J))THEN
                NNSMIN=I
                NNSMAX=I
                NNCMIN=J
                NNCMAX=J
                SMIN=A(J,I)
                SMAX=SMIN
                GOTO 30
              END IF
            END DO
          END IF
        END DO
30      CONTINUE
C
        MEAN=0.D0
        DO I=1,NSCAN
          IF(IFSCAN(I))THEN
            DO J=1,NCHAN
              IF(IFCHAN(J))THEN
                MEAN=MEAN+DBLE(A(J,I))
                IF(A(J,I).GT.SMAX)THEN
                  SMAX=A(J,I)
                  NNSMAX=I
                  NNCMAX=J
                END IF
                IF(A(J,I).LT.SMIN)THEN
                  SMIN=A(J,I)
                  NNSMIN=I
                  NNCMIN=J
                END IF
              END IF
            END DO
          END IF
        END DO
        MEAN=MEAN/DBLE(NPIX)
C
        IF(NPIX.LT.2)THEN
          SIGMA=0.D0
          GOTO 40
        END IF
        SIGMA=0.D0
        DO I=1,NSCAN
          IF(IFSCAN(I))THEN
            DO J=1,NCHAN
              IF(IFCHAN(J))THEN
                SIGMA=SIGMA+(DBLE(A(J,I))-MEAN)*(DBLE(A(J,I))-MEAN)
              END IF
            END DO
          END IF
        END DO
        SIGMA=DSQRT(SIGMA/DBLE(NPIX-1))
C
40      CONTINUE
C
        IF(CERR.EQ.'y')THEN
          MEANSN=0.D0
          NPIXT=0
          DO I=1,NSCAN
            IF(IFSCAN(I))THEN
              DO J=1,NCHAN
                IF(IFCHAN(J))THEN
                  IF(ERR(J,I).GT.0.)THEN
                    NPIXT=NPIXT+1
                    MEANSN=MEANSN+DBLE(A(J,I)/ERR(J,I))
                  END IF
                END IF
              END DO
            END IF
          END DO
          MEANSN=MEANSN/DBLE(NPIXT)
          SIGMASN=0.D0
          DO I=1,NSCAN
            IF(IFSCAN(I))THEN
              DO J=1,NCHAN
                IF(IFCHAN(J))THEN
                  IF(ERR(J,I).GT.0.)THEN
                    SIGMASN=SIGMASN+(DBLE(A(J,I)/ERR(J,I))-MEANSN)**2
                  END IF
                END IF
              END DO
            END IF
          END DO
          IF(NPIXT.GT.1)THEN
            SIGMASN=DSQRT(SIGMASN/DBLE(NPIXT-1))
          ELSE
            SIGMASN=0.D0
          END IF
        END IF
C
        WRITE(*,101)' ...OK!'
        WRITE(*,*)
C------------------------------------------------------------------------------
        WRITE(CDUMMY,'(I10,A1,I10)')NPIX,'/',NSCAN*NCHAN
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,101)'* Number of pixels employed/Total number: '
     +   //CDUMMY(1:L)
C
        WRITE(CDUMMY,*) SMAX
        CALL RMBLANK(CDUMMY,CDUMMY,L1)
        WRITE(*,100)'> Maximum.........................: '//
     +   CDUMMY(1:L1)
        WRITE(*,100)'     in pixel (x,y): '
        WRITE(CDUMMY,'(I10,A1,I10)')NNCMAX,',',NNSMAX
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,101)CDUMMY(1:L)
C
        WRITE(CDUMMY,*) SMIN
        CALL RMBLANK(CDUMMY,CDUMMY,L2)
        WRITE(*,100)'> Minimum.........................: '//
     +   CDUMMY(1:L2)
        DO I=1,5-(L2-L1)
          WRITE(*,100)' '
        END DO
        WRITE(*,100)'in pixel (x,y): '
        WRITE(CDUMMY,'(I10,A1,I10)')NNCMIN,',',NNSMIN
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,101)CDUMMY(1:L)
C
        WRITE(CDUMMY,*) REAL(MEAN)
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,101)'> Mean............................: '//CDUMMY(1:L)
C
        WRITE(CDUMMY,*) REAL(SIGMA)
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,101)'> Sigma...........................: '//CDUMMY(1:L)
C
        WRITE(CDUMMY,*) FMEDIAN
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,101)'> Median..........................: '//CDUMMY(1:L)
C
        WRITE(CDUMMY,*) FMEANSIGMA
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,101)'> Mean (removing pixels > 3 sigma): '//CDUMMY(1:L)
C
        WRITE(*,*)
        IF(CERR.EQ.'y')THEN
          WRITE(*,101)'* Signal-to-noise ratio:'
          WRITE(CDUMMY,*) REAL(MEANSN)
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          WRITE(*,101)'> Mean............................: '//
     +     CDUMMY(1:L)
          WRITE(CDUMMY,*) REAL(SIGMASN)
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          WRITE(*,101)'> Sigma...........................: '//
     +     CDUMMY(1:L)
          IF(NPIX.NE.NPIXT) WRITE(*,101)'NOTE: this last statistic is'//
     +     ' uncertain (there are pixels with errors <= 0)'
          WRITE(*,*)
        END IF
C
99      WRITE(*,101)'(1) change file'
        WRITE(*,101)'(2) change limits'
        WRITE(*,101)'(0) STOP'
        WRITE(*,100)'Option '
        IOPC=READI('0')
        IF(IOPC.EQ.0)THEN
          STOP
        ELSE IF(IOPC.EQ.1)THEN
          GOTO 10
        ELSE IF(IOPC.EQ.2)THEN
          GOTO 20
        ELSE
          WRITE(*,101)'ERROR: invalid option. Try again.'
          GOTO 99
        END IF
C
100     FORMAT(A,$)
101     FORMAT(A)
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
