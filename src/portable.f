C------------------------------------------------------------------------------
C Version 25-March-1997                                         file: portate.f
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
C Program: portable
C Classification: input/output
C Description: Transforms REDUCEME images into a more portable ASCII format.
C
Comment
C
        PROGRAM PORTABLE
        IMPLICIT NONE
C
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
C
        INTEGER I,J
        REAL A(NCMAX,NSMAX),ERR(NCMAX,NSMAX)
        CHARACTER*1 COPC,CERR
        CHARACTER*75 INFILE,ERRFILE,OUTFILE
C
        COMMON/BLKDATA1/NSCAN,NCHAN
        COMMON/BLKDATA2/STWV,DISP
        COMMON/BLKDATA3/A,ERR
C------------------------------------------------------------------------------
        THISPROGRAM='portate'
        CALL WELCOME('25-March-1997')
C
        WRITE(*,100)'Work with error images (y/n) '
        CERR(1:1)=READC('y','yn')
C
        WRITE(*,*)
        WRITE(*,101)'(1) REDUCEME format -> portable format'
        WRITE(*,101)'(2) portable format -> REDUCEME format'
        WRITE(*,100)'Option (1/2) '
        COPC(1:1)=READC('@','12')
        WRITE(*,*)
C------------------------------------------------------------------------------
        IF(COPC.EQ.'1')THEN
          WRITE(*,100)'Input file name'
          INFILE=INFILEX(10,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
          DO I=1,NSCAN
            READ(10) (A(J,I),J=1,NCHAN)
          END DO
          CLOSE(10)
          IF(CERR.EQ.'y')THEN
            WRITE(*,100)'Input error file name '
            CALL GUESSEF(INFILE,ERRFILE)
            ERRFILE=INFILEX(11,ERRFILE,NSCAN,NCHAN,STWV,DISP,21,.TRUE.)!match
            DO I=1,NSCAN
              READ(11) (ERR(J,I),J=1,NCHAN)
            END DO
            CLOSE(11)
          END IF
          WRITE(*,100)'Output file name'
          OUTFILE=OUTFILEX(20,'@',0,0,0.,0.,3,.FALSE.)
          CALL WRITE_ASCII(20,.FALSE.)
          CLOSE(20)
          IF(CERR.EQ.'y')THEN
            WRITE(*,100)'Output error file name '
            CALL GUESSEF(OUTFILE,ERRFILE)
            ERRFILE=OUTFILEX(21,ERRFILE,0,0,0.,0.,3,.FALSE.)
            CALL WRITE_ASCII(21,.TRUE.)
            CLOSE(21)
          END IF
C------------------------------------------------------------------------------
        ELSE
          WRITE(*,100)'Input file name'
          INFILE=INFILEX(10,'@',0,0,0.,0.,3,.FALSE.)
          CALL READ_ASCII(10,.FALSE.)
          CLOSE(10)
          IF(CERR.EQ.'y')THEN
            WRITE(*,100)'Input error file name '
            CALL GUESSEF(INFILE,ERRFILE)
            ERRFILE=INFILEX(11,ERRFILE,0,0,0.,0.,3,.FALSE.)
            CALL READ_ASCII(11,.TRUE.)
            CLOSE(11)
          END IF
          WRITE(*,100)'Output file name'
          OUTFILE=OUTFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
          DO I=1,NSCAN
            WRITE(20) (A(J,I),J=1,NCHAN)
          END DO
          CLOSE(20)
          IF(CERR.EQ.'y')THEN
            WRITE(*,100)'Output error file name '
            CALL GUESSEF(OUTFILE,ERRFILE)
            ERRFILE=OUTFILEX(21,ERRFILE,NSCAN,NCHAN,STWV,DISP,1,.TRUE.)
            DO I=1,NSCAN
              WRITE(21) (ERR(J,I),J=1,NCHAN)
            END DO
            CLOSE(21)
          END IF
C------------------------------------------------------------------------------
        END IF
C------------------------------------------------------------------------------
        STOP
100     FORMAT(A,$)
101     FORMAT(A)
        END
C
C******************************************************************************
C
        SUBROUTINE WRITE_ASCII(UNIT,LERR)
        IMPLICIT NONE
        INTEGER UNIT
        LOGICAL LERR
C
        INCLUDE 'redlib.inc'
        INTEGER TRUELEN
C
        INTEGER I,J
        INTEGER NCHAR
        REAL A(NCMAX,NSMAX),ERR(NCMAX,NSMAX)
C
        COMMON/BLKDATA1/NSCAN,NCHAN
        COMMON/BLKDATA2/STWV,DISP
        COMMON/BLKDATA3/A,ERR
C------------------------------------------------------------------------------
        WRITE(UNIT,*)NSCAN,NCHAN
        WRITE(UNIT,*)STWV,DISP
        WRITE(UNIT,*)AIRMASS 
        WRITE(UNIT,*)TIMEXPOS
        NCHAR=TRUELEN(OBJECT)
        WRITE(UNIT,*) NCHAR
        IF(NCHAR.GT.0)THEN
          WRITE(UNIT,101) OBJECT(1:NCHAR)
        END IF
        IF(LERR)THEN
          WRITE(UNIT,*) 0
        ELSE
          WRITE(UNIT,*) 1
        END IF
        NCHAR=TRUELEN(FITSFILE)
        WRITE(UNIT,*)NCHAR
        IF(NCHAR.GT.0)THEN
          WRITE(UNIT,101) FITSFILE(1:NCHAR)
        END IF
        NCHAR=TRUELEN(COMMENT)
        WRITE(UNIT,*)NCHAR
        IF(NCHAR.GT.0)THEN
          WRITE(UNIT,101) COMMENT(1:NCHAR)
        END IF
        IF(LERR)THEN
          DO I=1,NSCAN
            DO J=1,NCHAN
              WRITE(UNIT,*) ERR(J,I)
            END DO
          END DO
        ELSE
          DO I=1,NSCAN
            DO J=1,NCHAN
              WRITE(UNIT,*) A(J,I)
            END DO
          END DO
        END IF
101     FORMAT(A)
        END
C
C******************************************************************************
C
        SUBROUTINE READ_ASCII(UNIT,LERR)
        IMPLICIT NONE
        INTEGER UNIT
        LOGICAL LERR
C
        INCLUDE 'redlib.inc'
C
        INTEGER I,J
        INTEGER NERROR,NCHAR
        REAL A(NCMAX,NSMAX),ERR(NCMAX,NSMAX)
C
        COMMON/BLKDATA1/NSCAN,NCHAN
        COMMON/BLKDATA2/STWV,DISP
        COMMON/BLKDATA3/A,ERR
C------------------------------------------------------------------------------
        READ(UNIT,*)NSCAN,NCHAN
        READ(UNIT,*)STWV,DISP
        READ(UNIT,*)AIRMASS 
        READ(UNIT,*)TIMEXPOS
        READ(UNIT,*) NCHAR
        IF(NCHAR.GT.0)THEN
          READ(UNIT,101) OBJECT(1:NCHAR)
        END IF
        READ(UNIT,*) NERROR
        IF(NERROR.EQ.0)THEN
          IF(.NOT.LERR)THEN
            WRITE(*,101)'FATAL ERROR1: something is wrong in READ_ASCII'
            STOP
          END IF
        ELSEIF(NERROR.EQ.1)THEN
          IF(LERR)THEN
            WRITE(*,101)'FATAL ERROR2: something is wrong in READ_ASCII'
            STOP
          END IF
        ELSE
          WRITE(*,101)'FATAL ERROR3: something is wrong in READ_ASCII'
          STOP
        END IF
        READ(UNIT,*)NCHAR
        IF(NCHAR.GT.0)THEN
          READ(UNIT,101) FITSFILE(1:NCHAR)
        END IF
        READ(UNIT,*)NCHAR
        IF(NCHAR.GT.0)THEN
          READ(UNIT,101) COMMENT(1:NCHAR)
        END IF
        IF(LERR)THEN
          DO I=1,NSCAN
            DO J=1,NCHAN
              READ(UNIT,*) ERR(J,I)
            END DO
          END DO
        ELSE
          DO I=1,NSCAN
            DO J=1,NCHAN
              READ(UNIT,*) A(J,I)
            END DO
          END DO
        END IF
101     FORMAT(A)
        END
