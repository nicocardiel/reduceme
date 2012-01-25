C------------------------------------------------------------------------------
C Version 7-May-1997                                          file:writeascii.f
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
C Program: writeascii
C Classification: input/output
C Description: Reads a file with REDUCEME format and creates a new file with
C ASCII format.
C
Comment
C
        PROGRAM WRITEASCII
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER READILIM
C
        INTEGER I,J
        INTEGER NS0
        REAL A(NCMAX,NSMAX)
        REAL ERR(NCMAX,NSMAX)
        REAL WV
        CHARACTER*1 CERR
        CHARACTER*75 FILENAME,ERRFILE,OUTFILE
C------------------------------------------------------------------------------
        THISPROGRAM='writeascii'
        CALL WELCOME('7-May-1997')
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
          ERRFILE=INFILEX(21,ERRFILE,NSCAN,NCHAN,STWV,DISP,21,.TRUE.) !...match
          DO I=1,NSCAN
            READ(21) (ERR(J,I),J=1,NCHAN)
          END DO
          CLOSE(21)
        END IF
C------------------------------------------------------------------------------
10      IF(NSCAN.GT.1)THEN
          WRITE(*,100)'Scan to be saved (0=EXIT) '
          NS0=READILIM('0',0,NSCAN)
        ELSE
          NS0=1
        END IF
        IF(NS0.EQ.0) GOTO 90
C salvamos fichero
        WRITE(*,100)'Output file name'
        OUTFILE=OUTFILEX(30,'@',0,0,0.,0.,3,.FALSE.)
        IF(CERR.EQ.'n')THEN
          DO J=1,NCHAN
            WV=STWV+REAL(J-1)*DISP
            WRITE(30,*)WV,A(J,NS0)
          END DO
        ELSE
          DO J=1,NCHAN
            WV=STWV+REAL(J-1)*DISP
            WRITE(30,*)WV,A(J,NS0),ERR(J,NS0)
          END DO
        END IF
        CLOSE(30)
        IF(NSCAN.GT.1) GOTO 10
C------------------------------------------------------------------------------
90      CONTINUE
        STOP
100     FORMAT(A,$) 
        END
