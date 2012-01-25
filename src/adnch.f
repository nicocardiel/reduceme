C------------------------------------------------------------------------------
C Version 28-November-1996                                         file:adnch.f
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
C Program: adnch
C Classification: arithmetic & manipulations
C Description: Add channels of an image, generating a new image with the 
C same NSCAN and variable NCHAN.
C
Comment
C
        PROGRAM ADNCH
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER READILIM
        REAL READF
C
        INTEGER I,J
        INTEGER K
        INTEGER N,N1,N2
        REAL S(NSMAX),A(NCMAX,NSMAX),ERR(NCMAX,NSMAX)
        REAL ES(NSMAX),B(NCMAX,NSMAX),ERRB(NCMAX,NSMAX)
        REAL FACTOR
        CHARACTER*1 CERR
        CHARACTER*75 FILENAME,ERRFILE,OUTFILE
C------------------------------------------------------------------------------
        THISPROGRAM='adnch'
        CALL WELCOME('28-November-1996')
        CALL SHOWHLP('explanation')
C
        WRITE(*,100)'Work with error images (y/n) '
        CERR(1:1)=READC('n','yn')
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
C------------------------------------------------------------------------------
        K=0
C
50      CALL SHOWHLP('selection')
        WRITE(*,100)'1st Channel (0=Exit) '
        N1=READILIM('@',0,NCHAN)
        IF(N1.EQ.0) GO TO 60
        WRITE(*,100)'No. of Channels '
        N=READILIM('@',1,NCHAN-N1+1)
        WRITE(*,100)'Factor'
        FACTOR=READF('@')
        N2=N1+N-1
        DO I=1,NSCAN
          S(I)=0.
        END DO
        DO I=1,NSCAN
          DO J=N1,N2
            S(I)=S(I)+A(J,I)*FACTOR
          END DO
        END DO
        IF(CERR.EQ.'y')THEN
          DO I=1,NSCAN
            ES(I)=0.
          END DO
          DO I=1,NSCAN
            DO J=N1,N2
              ES(I)=ES(I)+ERR(J,I)*ERR(J,I)
            END DO
          END DO
        END IF
        K=K+1
        IF(K.GT.NCMAX)THEN
          WRITE(*,101)'ERROR: maximum number of channels excedeed'
          STOP
        END IF
        DO I=1,NSCAN
          B(K,I)=S(I)
        END DO
        IF(CERR.EQ.'y')THEN
          DO I=1,NSCAN
            ERRB(K,I)=FACTOR*SQRT(ES(I))
          END DO
        END IF
        GO TO 50
  60    CONTINUE
C------------------------------------------------------------------------------
        WRITE(*,100)'Output file name'
        OUTFILE=OUTFILEX(30,'@',NSCAN,K,0.,0.,1,.FALSE.)
        DO I=1,NSCAN
          WRITE(30) (B(J,I),J=1,K)
        END DO
        CLOSE(30)
        IF(CERR.EQ.'y')THEN
          WRITE(*,100)'Output error file name '
          CALL GUESSEF(OUTFILE,ERRFILE)
          OUTFILE=OUTFILEX(31,ERRFILE,NSCAN,K,0.,0.,1,.TRUE.)
          DO I=1,NSCAN
            WRITE(31) (ERRB(J,I),J=1,K)
          END DO
          CLOSE(31)
        END IF
C
        WRITE(*,110)'No. of saved channels: ',K
        STOP
100     FORMAT(A,$)
101     FORMAT(A)
110     FORMAT(A,I6)
        END
