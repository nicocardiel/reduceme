!------------------------------------------------------------------------------
! Version 20-November-2026                                     file: rotate.f90
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
! Program: rotate
! Classification: arithmetic & manipulations
! Description: Rotates an image.
!
!Comment
!
! Rota e invierte direcciones en una imagen
!
        PROGRAM ROTATE
        USE RED_M
        USE FUTILS_M
        USE REDLIB_M
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
!
        INTEGER I,J
        INTEGER NEW_NSCAN,NEW_NCHAN
        INTEGER IOPC
        REAL A(NCMAX,NSMAX),AERR(NCMAX,NSMAX)
        REAL B(NCMAX,NSMAX),BERR(NCMAX,NSMAX)
        CHARACTER*1 CERR
        CHARACTER*255 INFILE,OUTFILE,ERRFILE
        LOGICAL LROTATE
!------------------------------------------------------------------------------
        NEW_NSCAN=0 !avoid compilation warning
        NEW_NCHAN=0 !avoid compilation warning
        THISPROGRAM='rotate'
        CALL WELCOME('20-November-2026')
!------------------------------------------------------------------------------
        WRITE(*,100,ADVANCE='NO')'Work with error images (y/n) '
        CERR(1:1)=READC('n','yn')
!------------------------------------------------------------------------------
        WRITE(*,100,ADVANCE='NO')'Input file name'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        DO I=1,NSCAN
          READ(20) (A(J,I),J=1,NCHAN)
        END DO
        CLOSE(20)
        IF(CERR.EQ.'y')THEN
          WRITE(*,100,ADVANCE='NO')'Error file name '
          CALL GUESSEF(INFILE,ERRFILE)
          ERRFILE=INFILEX(30,ERRFILE,NSCAN,NCHAN,STWV,DISP,21,.TRUE.) !match
          DO I=1,NSCAN
            READ(30) (AERR(J,I),J=1,NCHAN)
          END DO
          CLOSE(30)
        END IF
!------------------------------------------------------------------------------
! determine whether the frame can be rotated or not
        LROTATE=.TRUE.
        IF(NSCAN.GT.NSMAX)THEN
          WRITE(*,101)'WARNING: this file cannot be rotated'
          WRITE(*,101)'>>> NSCAN.GT.NSMAX.'
          LROTATE=.FALSE.
        END IF
        IF(NCHAN.GT.NCMAX)THEN
          WRITE(*,101)'WARNING: this file cannot be rotated'
          WRITE(*,101)'>>> NCHAN.GT.NCMAX.'
          LROTATE=.FALSE.
        END IF
!
!------------------------------------------------------------------------------
10      WRITE(*,*)
        WRITE(*,101)'(1) Reverse x direction'
        WRITE(*,101)'(2) Reverse y direction'
        IF(LROTATE)THEN
          WRITE(*,101)'(3) Rotate +90 degrees'
          WRITE(*,101)'(4) Rotate -90 degrees'
        END IF
        WRITE(*,101)'(0) Exit & Save'
        WRITE(*,100,ADVANCE='NO')'Option '
        IF(LROTATE)THEN
          IOPC=READILIM('0',0,4)
        ELSE
          IOPC=READILIM('0',0,2)
        END IF
!------------------------------------------------------------------------------
        IF(IOPC.EQ.0)THEN
          GOTO 90
!------------------------------------------------------------------------------
        ELSEIF(IOPC.EQ.1)THEN
          WRITE(*,100,ADVANCE='NO')'Reversing...'
          DO I=1,NSCAN
            DO J=1,NCHAN
              B(J,I)=A(NCHAN+1-J,I)
            END DO
          END DO
          IF(CERR.EQ.'y')THEN
            DO I=1,NSCAN
              DO J=1,NCHAN
                BERR(J,I)=AERR(NCHAN+1-J,I)
              END DO
            END DO
          END IF
          WRITE(*,101)'  ...OK!'
          NEW_NSCAN=NSCAN
          NEW_NCHAN=NCHAN
!------------------------------------------------------------------------------
        ELSEIF(IOPC.EQ.2)THEN
          WRITE(*,100,ADVANCE='NO')'Reversing...'
          DO I=1,NSCAN
            DO J=1,NCHAN
              B(J,I)=A(J,NSCAN-I+1)
            END DO
          END DO
          IF(CERR.EQ.'y')THEN
            DO I=1,NSCAN
              DO J=1,NCHAN
                BERR(J,I)=AERR(J,NSCAN-I+1)
              END DO
            END DO
          END IF
          WRITE(*,101)'  ...OK!'
          NEW_NSCAN=NSCAN
          NEW_NCHAN=NCHAN
!------------------------------------------------------------------------------
        ELSEIF(IOPC.EQ.3)THEN
            WRITE(*,100,ADVANCE='NO')'Rotating...'
          DO I=1,NSCAN
            DO J=1,NCHAN
              B(NSCAN-I+1,J)=A(J,I)
            END DO
          END DO
          IF(CERR.EQ.'y')THEN
            DO I=1,NSCAN
              DO J=1,NCHAN
                BERR(NSCAN-I+1,J)=AERR(J,I)
              END DO
            END DO
            WRITE(*,101)'  ...OK!'
          END IF
          NEW_NSCAN=NCHAN
          NEW_NCHAN=NSCAN
!------------------------------------------------------------------------------
        ELSEIF(IOPC.EQ.4)THEN
            WRITE(*,100,ADVANCE='NO')'Rotating...'
          DO I=1,NSCAN
            DO J=1,NCHAN
              B(I,NCHAN-J+1)=A(J,I)
            END DO
          END DO
          IF(CERR.EQ.'y')THEN
            DO I=1,NSCAN
              DO J=1,NCHAN
                BERR(I,NCHAN-J+1)=AERR(J,I)
              END DO
            END DO
            WRITE(*,101)'  ...OK!'
          END IF
          NEW_NSCAN=NCHAN
          NEW_NCHAN=NSCAN
!------------------------------------------------------------------------------
        END IF
!------------------------------------------------------------------------------
! Update the original frames
        WRITE(*,100,ADVANCE='NO')'Updating...'
        NSCAN=NEW_NSCAN
        NCHAN=NEW_NCHAN
        DO I=1,NSCAN
          DO J=1,NCHAN
            A(J,I)=B(J,I)
          END DO
        END DO
        IF(CERR.EQ.'y')THEN
          DO I=1,NSCAN
            DO J=1,NCHAN
              AERR(J,I)=BERR(J,I)
            END DO
          END DO
        END IF
        WRITE(*,101)'  ...OK!'
        GOTO 10
!------------------------------------------------------------------------------
90      WRITE(*,100,ADVANCE='NO')'Output file name'
        OUTFILE=OUTFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        DO I=1,NSCAN
          WRITE(20) (A(J,I),J=1,NCHAN)
        END DO
        CLOSE(20)
!
        IF(CERR.EQ.'y')THEN
          WRITE(*,100,ADVANCE='NO')'Output error file name '
          CALL GUESSEF(OUTFILE,ERRFILE)
          ERRFILE=OUTFILEX(30,ERRFILE,NSCAN,NCHAN,STWV,DISP,1,.TRUE.)
          DO I=1,NSCAN
            WRITE(30) (AERR(J,I),J=1,NCHAN)
          END DO
          CLOSE(30)
        END IF
!------------------------------------------------------------------------------
        STOP
100     FORMAT(A)
101     FORMAT(A)
        END
