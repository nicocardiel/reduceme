!------------------------------------------------------------------------------
! Version 20-November-2026                                    file: portate.f90
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
! Program: portable
! Classification: input/output
! Description: Transforms REDUCEME images into a more portable ASCII format.
!
!Comment
!
        PROGRAM PORTABLE
        USE REDLIB_M
        IMPLICIT NONE
!
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
!
        INTEGER I,J
        REAL A(NCMAX,NSMAX),ERR(NCMAX,NSMAX)
        CHARACTER*1 COPC,CERR
        CHARACTER*255 INFILE,ERRFILE,OUTFILE
!
        COMMON/BLKDATA1/NSCAN,NCHAN
        COMMON/BLKDATA2/STWV,DISP
        COMMON/BLKDATA3/A,ERR
!------------------------------------------------------------------------------
        THISPROGRAM='portate'
        CALL WELCOME('20-November-2026')
!
        WRITE(*,100,ADVANCE='NO')'Work with error images (y/n) '
        CERR(1:1)=READC('y','yn')
!
        WRITE(*,*)
        WRITE(*,101)'(1) REDUCEME format -> portable format'
        WRITE(*,101)'(2) portable format -> REDUCEME format'
        WRITE(*,100,ADVANCE='NO')'Option (1/2) '
        COPC(1:1)=READC('@','12')
        WRITE(*,*)
!------------------------------------------------------------------------------
        IF(COPC.EQ.'1')THEN
          WRITE(*,100,ADVANCE='NO')'Input file name'
          INFILE=INFILEX(10,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
          DO I=1,NSCAN
            READ(10) (A(J,I),J=1,NCHAN)
          END DO
          CLOSE(10)
          IF(CERR.EQ.'y')THEN
            WRITE(*,100,ADVANCE='NO')'Input error file name '
            CALL GUESSEF(INFILE,ERRFILE)
            ERRFILE=INFILEX(11,ERRFILE,NSCAN,NCHAN,STWV,DISP,21,.TRUE.)!match
            DO I=1,NSCAN
              READ(11) (ERR(J,I),J=1,NCHAN)
            END DO
            CLOSE(11)
          END IF
          WRITE(*,100,ADVANCE='NO')'Output file name'
          OUTFILE=OUTFILEX(20,'@',0,0,0.,0.,3,.FALSE.)
          CALL WRITE_ASCII(20,.FALSE.)
          CLOSE(20)
          IF(CERR.EQ.'y')THEN
            WRITE(*,100,ADVANCE='NO')'Output error file name '
            CALL GUESSEF(OUTFILE,ERRFILE)
            ERRFILE=OUTFILEX(21,ERRFILE,0,0,0.,0.,3,.FALSE.)
            CALL WRITE_ASCII(21,.TRUE.)
            CLOSE(21)
          END IF
!------------------------------------------------------------------------------
        ELSE
          WRITE(*,100,ADVANCE='NO')'Input file name'
          INFILE=INFILEX(10,'@',0,0,0.,0.,3,.FALSE.)
          CALL READ_ASCII(10,.FALSE.)
          CLOSE(10)
          IF(CERR.EQ.'y')THEN
            WRITE(*,100,ADVANCE='NO')'Input error file name '
            CALL GUESSEF(INFILE,ERRFILE)
            ERRFILE=INFILEX(11,ERRFILE,0,0,0.,0.,3,.FALSE.)
            CALL READ_ASCII(11,.TRUE.)
            CLOSE(11)
          END IF
          WRITE(*,100,ADVANCE='NO')'Output file name'
          OUTFILE=OUTFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
          DO I=1,NSCAN
            WRITE(20) (A(J,I),J=1,NCHAN)
          END DO
          CLOSE(20)
          IF(CERR.EQ.'y')THEN
            WRITE(*,100,ADVANCE='NO')'Output error file name '
            CALL GUESSEF(OUTFILE,ERRFILE)
            ERRFILE=OUTFILEX(21,ERRFILE,NSCAN,NCHAN,STWV,DISP,1,.TRUE.)
            DO I=1,NSCAN
              WRITE(21) (ERR(J,I),J=1,NCHAN)
            END DO
            CLOSE(21)
          END IF
!------------------------------------------------------------------------------
        END IF
!------------------------------------------------------------------------------
        STOP
100     FORMAT(A)
101     FORMAT(A)
        END
!
!******************************************************************************
!
        SUBROUTINE WRITE_ASCII(UNIT,LERR)
        USE REDLIB_M
        IMPLICIT NONE
        INTEGER UNIT
        LOGICAL LERR
!
        INCLUDE 'redlib.inc'
        INTEGER TRUELEN
!
        INTEGER I,J
        INTEGER NCHAR
        REAL A(NCMAX,NSMAX),ERR(NCMAX,NSMAX)
!
        COMMON/BLKDATA1/NSCAN,NCHAN
        COMMON/BLKDATA2/STWV,DISP
        COMMON/BLKDATA3/A,ERR
!------------------------------------------------------------------------------
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
!
!******************************************************************************
!
        SUBROUTINE READ_ASCII(UNIT,LERR)
        USE REDLIB_M
        IMPLICIT NONE
        INTEGER UNIT
        LOGICAL LERR
!
        INCLUDE 'redlib.inc'
!
        INTEGER I,J
        INTEGER NERROR,NCHAR
        REAL A(NCMAX,NSMAX),ERR(NCMAX,NSMAX)
!
        COMMON/BLKDATA1/NSCAN,NCHAN
        COMMON/BLKDATA2/STWV,DISP
        COMMON/BLKDATA3/A,ERR
!------------------------------------------------------------------------------
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
