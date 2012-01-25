C------------------------------------------------------------------------------
C Version 28-November-1996                                         file:addnf.f
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
C Program: addnf
C Classification: arithmetic & manipulations
C Description: Add several images, taking into account offsets in the spatial 
C direction.
C
Comment
C
        PROGRAM ADDNF
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER READI
        INTEGER READILIM
C
        INTEGER I,J
        INTEGER N
        INTEGER II,IN
        INTEGER NF
        INTEGER INMAX,INMIN
        REAL S(NCMAX),A(NCMAX,NSMAX),ES(NCMAX),ERR(NCMAX,NSMAX)
        REAL FNF
        CHARACTER*1 CERR,CNOR
        CHARACTER*75 INFILE,ERRFILE,OUTFILE
C------------------------------------------------------------------------------
        THISPROGRAM='addnf'
        CALL WELCOME('28-November-1996')
        CALL SHOWHLP('explanation')
C
        WRITE(*,100)'Work with error images (y/n) '
        CERR(1:1)=READC('n','yn')
C
        DO I=1,NSMAX
          DO J=1,NCMAX
            A(J,I)=0.0
          END DO
        END DO
        IF(CERR.EQ.'y')THEN
          DO I=1,NSMAX
            DO J=1,NCMAX
              ERR(J,I)=0.0
            END DO
          END DO
        END IF
C
        WRITE(*,100)'Total no. of frames to be added'
        NF=READILIM('@',1,9999)
C
        INMAX=0
        INMIN=0
C------------------------------------------------------------------------------
        DO N=1,NF
          WRITE(*,*)
          WRITE(*,'(A,I3,A,$)')'#',N,'   > Input file'
          IF(N.EQ.1)THEN
            INFILE=INFILEX(14,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
          ELSE
            INFILE=INFILEX(14,'@',NSCAN,NCHAN,STWV,DISP,21,.FALSE.) !.....match
          END IF
          IF(CERR.EQ.'y')THEN
            WRITE(*,'(A,I3,A,$)')'#',N,'   > Error file '
            CALL GUESSEF(INFILE,ERRFILE)
            ERRFILE=INFILEX(15,ERRFILE,NSCAN,NCHAN,STWV,DISP,21,.TRUE.)!..match
          END IF
          CALL SHOWHLP('offset')
          WRITE(*,100)'Increment (final scan = old scan - increment) '
          IN=READI('0')
          IF(IN.LE.INMAX) GO TO 50
          INMAX=IN
  50      IF(IN.GE.INMIN) GO TO 60
          INMIN=IN
  60      CONTINUE
          DO 30 I=1,NSCAN
            READ(14) (S(J),J=1,NCHAN) 
            IF(CERR.EQ.'y') READ(15) (ES(J),J=1,NCHAN)
            II=I-IN
            IF(II.LT.1.OR.II.GT.NSCAN) GO TO 30
            DO J=1,NCHAN
              A(J,II)=A(J,II)+S(J)
            END DO
            IF(CERR.EQ.'y')THEN
              DO J=1,NCHAN
                ERR(J,II)=ERR(J,II)+ES(J)*ES(J)
              END DO
            END IF
30          CONTINUE
          CLOSE(14)
          IF(CERR.EQ.'y') CLOSE(15)
        END DO
C------------------------------------------------------------------------------
        IF(INMAX.EQ.0.AND.INMIN.EQ.0) GO TO 70
        IF(INMIN.LT.0) THEN
          DO I=1,IABS(INMIN)
            DO J=1,NCHAN
              A(J,I)=0.0
            END DO
            IF(CERR.EQ.'y')THEN
              DO J=1,NCHAN
                ERR(J,I)=0.0
              END DO
            END IF
          END DO
        END IF
        IF(INMAX.GT.0) THEN
          DO I=NSCAN-INMAX+1,NSCAN
            DO J=1,NCHAN
              A(J,I)=0.0
            END DO
            IF(CERR.EQ.'y')THEN
              DO J=1,NCHAN
                ERR(J,I)=0.0
              END DO
            END IF
          END DO
        END IF
70      CONTINUE
C------------------------------------------------------------------------------
        IF(NF.GT.1)THEN                     !alteramos descriptores de cabecera
          AIRMASS=0.
          TIMEXPOS=0.
          OBJECT='[from addnf]'
          FITSFILE='[from addnf]'
          CALL SHOWHLP('normalize')
          WRITE(*,100)'Normalize the added frame (y/n) '
          CNOR(1:1)=READC('y','yn')
        ELSE
          CNOR='n'
        END IF
        IF(CNOR.EQ.'y')THEN
          FNF=REAL(NF)
          DO I=1,NSCAN
            DO J=1,NCHAN
              A(J,I)=A(J,I)/FNF
            END DO
          END DO
        ELSE
          FNF=1.0
        END IF
C
        WRITE(*,*)
        WRITE(*,100)'Output file name'
        OUTFILE=OUTFILEX(16,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        WRITE(*,100)'Saving file...'
        DO I=1,NSCAN
          WRITE(16) (A(J,I),J=1,NCHAN)
        END DO
        CLOSE(16)
        WRITE(*,101)'  ...OK! File saved and closed.'
C
        IF(CERR.EQ.'y')THEN
          WRITE(*,100)'Output error file name '
          CALL GUESSEF(OUTFILE,ERRFILE)
          OUTFILE=OUTFILEX(17,ERRFILE,NSCAN,NCHAN,STWV,DISP,1,.TRUE.)
          WRITE(*,100)'Saving file...'
          DO I=1,NSCAN
            DO J=1,NCHAN
              ERR(J,I)=SQRT(ERR(J,I))/FNF
            END DO
            WRITE(17) (ERR(J,I),J=1,NCHAN)
          END DO
          CLOSE(17)
          WRITE(*,101)'  ...OK! File saved and closed.'
        END IF
C
        STOP
C
100     FORMAT(A,$)
101     FORMAT(A)
        END
