C------------------------------------------------------------------------------
C Version 28-November-1996                                         file:adnsc.f
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
C Program: adnsc
C Classification: arithmetic & manipulations
C Description: Add scans of an image, generating a new image with the same 
C NCHAN and variable NSCAN.
C
Comment
C
        PROGRAM ADNSC
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER READILIM
        REAL READF
C
        INTEGER I,J,K,L
        INTEGER N,N1,N2
        REAL S(NCMAX),A(NCMAX,NSMAX),B(NCMAX,NSMAX)
        REAL ES(NCMAX),ERR(NCMAX,NSMAX),ERRB(NCMAX,NSMAX)
        REAL FACTOR
        CHARACTER*1 CERR,CSINGLE,CNOR,CFHEAD2
        CHARACTER*15 CSCANBINNING
        CHARACTER*75 FILENAME,ERRFILE,HEADFILE,OUTFILE
C------------------------------------------------------------------------------
        THISPROGRAM='adnsc'
        CALL WELCOME('28-November-1996')
        CALL SHOWHLP('explanation')
C
        WRITE(*,100)'Work with error images (y/n) '
        CERR(1:1)=READC('y','yn')
C
        CALL SHOWHLP('single spectrum')
        WRITE(*,100)'Output into a single spectrum (y/n) '
        CSINGLE(1:1)=READC('y','yn')
C
        WRITE(*,100)'Input file name'
        FILENAME=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        DO I=1,NSCAN
          READ(20) (A(J,I),J=1,NCHAN)
        END DO
        CLOSE(20)
        IF(CERR.EQ.'y')THEN
          WRITE(*,100)'Input error file name '
          CALL GUESSEF(FILENAME,ERRFILE)
          ERRFILE=INFILEX(21,ERRFILE,NSCAN,NCHAN,STWV,DISP,21,.TRUE.)!....match
          DO I=1,NSCAN
            READ(21) (ERR(J,I),J=1,NCHAN)
          END DO
          CLOSE(21)
        END IF
        IF(CSINGLE.EQ.'n')THEN
          CALL SHOWHLP('file to index')
          WRITE(*,100)'Create log file (with binning regions) '//
     +     'for index.f (y/n) '
          CFHEAD2(1:1)=READC('n','yn')
          IF(CFHEAD2.EQ.'y')THEN
            WRITE(*,100)'Log file name'
            HEADFILE=OUTFILEX(27,'@',0,0,0.,0.,3,.FALSE.)
          END IF
          GOTO 40
        END IF
C------------------------------------------------------------------------------
C Salida en un unico espectro
        DO J=1,NCHAN
          S(J)=0.
        END DO
        IF(CERR.EQ.'y')THEN
          DO J=1,NCHAN
            ES(J)=0.
          END DO
        END IF
C
        WRITE(*,*)
        CALL SHOWHLP('factors versus normalization')
        WRITE(*,101)'* If you are using different factors for '//
     +   'each spectrum, DO NOT normalize'
        WRITE(*,100)'Normalize final spectrum (y/n) '
        CNOR(1:1)=READC('y','yn')
        IF(CNOR.EQ.'y')THEN
          FACTOR=1.
        END IF
C
        N=0
10      CALL SHOWHLP('enter scan region')
        WRITE(*,100)'Scan region (0,0=EXIT) '
        CALL READ2I('0,0',N1,N2)
        IF((N1.EQ.0).AND.(N2.EQ.0)) GOTO 20
        IF((N1.LT.0).OR.(N2.GT.NSCAN).OR.(N1.GT.N2))THEN
          WRITE(*,101)'ERROR: invalid entry. Try again.'
          GOTO 10
        END IF
        IF(CNOR.EQ.'n')THEN
          WRITE(*,100)'Factor '
          FACTOR=READF('1.0')
        END IF
        DO I=N1,N2
          DO J=1,NCHAN
            S(J)=S(J)+A(J,I)*FACTOR
          END DO
          IF(CERR.EQ.'y')THEN
            DO J=1,NCHAN
              ES(J)=ES(J)+ERR(J,I)*ERR(J,I)*FACTOR*FACTOR
            END DO
          END IF
        END DO
        N=N+(N2-N1+1)
        WRITE(*,110)'>>> Number of scans added......: ',N2-N1+1
        WRITE(*,110)'>>> Total number of scans added: ',N
        GOTO 10
20      IF(N.EQ.0)THEN
          WRITE(*,101)'ERROR: number of scans added = 0. Try again.'
          GOTO 10
        END IF
        IF(CNOR.EQ.'y')THEN
          DO J=1,NCHAN
            S(J)=S(J)/REAL(N)
          END DO
          IF(CERR.EQ.'y')THEN
            DO J=1,NCHAN
              ES(J)=SQRT(ES(J))/REAL(N)
            END DO
          END IF
        ELSE
          IF(CERR.EQ.'y')THEN
            DO J=1,NCHAN
              ES(J)=SQRT(ES(J))
            END DO
          END IF
        END IF
        K=1
        GOTO 70
C------------------------------------------------------------------------------
40      K=0
50      CALL SHOWHLP('scan region with increment')
        WRITE(*,100)'1st Scan (0=Exit) '
        N1=READILIM('@',0,NSCAN)
        IF(N1.EQ.0) GO TO 60
        WRITE(*,100)'No. of Scans '
        N=READILIM('@',1,NSCAN-N1+1)
        WRITE(*,100)'Factor '
        FACTOR=READF('1.0')
        N2=N1+N-1
C
        IF(CFHEAD2.EQ.'y')THEN
          WRITE(CSCANBINNING,'(I7,A1,I7)')N1,',',N2
          CALL RMBLANK(CSCANBINNING,CSCANBINNING,L)
          WRITE(27,'(A15)') CSCANBINNING(1:L)
        END IF
C
        DO J=1,NCHAN
          S(J)=0.0
        END DO
        DO I=N1,N2
          DO J=1,NCHAN
            S(J)=S(J)+A(J,I)*FACTOR
          END DO
        END DO
        IF(CERR.EQ.'y')THEN
          DO J=1,NCHAN
            ES(J)=0.0
          END DO
          DO I=N1,N2
            DO J=1,NCHAN
              ES(J)=ES(J)+ERR(J,I)*ERR(J,I)
            END DO
          END DO
        END IF
C
        K=K+1
        DO J=1,NCHAN
          B(J,K)=S(J)
        END DO
        IF(CERR.EQ.'y')THEN
          DO J=1,NCHAN
            ERRB(J,K)=FACTOR*SQRT(ES(J))
          END DO
        END IF
        IF(K.EQ.NSMAX)THEN
          WRITE(*,101)'WARNING: No. of scans = NSMAX'
          WRITE(*,101)'         Output file must be saved.'
          GOTO 60
        END IF
        GO TO 50
60      IF(CFHEAD2.EQ.'y') CLOSE(27)
C------------------------------------------------------------------------------
C salvamos fichero
70      WRITE(*,100)'Output file name'
        OUTFILE=OUTFILEX(30,'@',K,NCHAN,STWV,DISP,1,.FALSE.)
        IF(CSINGLE.EQ.'y')THEN
          WRITE(30) (S(J),J=1,NCHAN)
        ELSE
          DO I=1,K
            WRITE(30) (B(J,I),J=1,NCHAN)
          END DO
        END IF
        CLOSE(30)
C
        IF(CERR.EQ.'y')THEN
          WRITE(*,100)'Output error file name '
          CALL GUESSEF(OUTFILE,ERRFILE)
          OUTFILE=OUTFILEX(31,ERRFILE,K,NCHAN,STWV,DISP,1,.TRUE.)
          IF(CSINGLE.EQ.'y')THEN
            WRITE(31) (ES(J),J=1,NCHAN)
          ELSE
            DO I=1,K
              WRITE(31) (ERRB(J,I),J=1,NCHAN)
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
