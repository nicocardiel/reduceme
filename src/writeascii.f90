!------------------------------------------------------------------------------
! Version 20-November-2026                                  file:writeascii.f90
!------------------------------------------------------------------------------
! Copyright N. Cardiel & J. Gorgas, Departamento de Astrofisica
! Universidad Complutense de Madrid, 28040-Madrid, Spain
! E-mail: ncl@astrax.fis.ucm.es or fjg@astrax.fis.ucm.es
!------------------------------------------------------------------------------
! This program is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the Free
! Software Foundation; either version 2 of the License, or (at your option) any
! later version. See the file gnu-public-license.txt for details.
!------------------------------------------------------------------------------
!Comment
!
! Program: writeascii
! Classification: input/output
! Description: Reads a file with REDUCEME format and creates a new file with
! ASCII format.
!
!Comment
!
        PROGRAM WRITEASCII
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER READILIM
!
        INTEGER I,J
        INTEGER NS0
        REAL A(NCMAX,NSMAX)
        REAL ERR(NCMAX,NSMAX)
        REAL WV
        CHARACTER*1 CERR
        CHARACTER*255 FILENAME,ERRFILE,OUTFILE
!------------------------------------------------------------------------------
        THISPROGRAM='writeascii'
        CALL WELCOME('20-November-2026')
!
        WRITE(*,100,ADVANCE='NO')'Work with error images (y/n) '
        CERR(1:1)=READC('n','yn')
!
        WRITE(*,100,ADVANCE='NO')'Input file name'
        FILENAME=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        DO I=1,NSCAN
          READ(20) (A(J,I),J=1,NCHAN)
        END DO
        CLOSE(20)
        IF(CERR.EQ.'y')THEN
          WRITE(*,100,ADVANCE='NO')'Input error file name '
          CALL GUESSEF(FILENAME,ERRFILE)
          ERRFILE=INFILEX(21,ERRFILE,NSCAN,NCHAN,STWV,DISP,21,.TRUE.) !...match
          DO I=1,NSCAN
            READ(21) (ERR(J,I),J=1,NCHAN)
          END DO
          CLOSE(21)
        END IF
!------------------------------------------------------------------------------
10      IF(NSCAN.GT.1)THEN
          WRITE(*,100,ADVANCE='NO')'Scan to be saved (0=EXIT) '
          NS0=READILIM('0',0,NSCAN)
        ELSE
          NS0=1
        END IF
        IF(NS0.EQ.0) GOTO 90
! salvamos fichero
        WRITE(*,100,ADVANCE='NO')'Output file name'
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
!------------------------------------------------------------------------------
90      CONTINUE
        STOP
100     FORMAT(A) 
        END
