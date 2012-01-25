C------------------------------------------------------------------------------
C Version 8-December-1996                                        file: rotate.f
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
C Program: rotate
C Classification: arithmetic & manipulations
C Description: Rotates an image.
C
Comment
C
C Rota e invierte direcciones en una imagen
C
        PROGRAM ROTATE
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER READILIM
C
        INTEGER I,J
        INTEGER NEW_NSCAN,NEW_NCHAN
        INTEGER IOPC
        REAL A(NCMAX,NSMAX),AERR(NCMAX,NSMAX)
        REAL B(NCMAX,NSMAX),BERR(NCMAX,NSMAX)
        CHARACTER*1 CERR
        CHARACTER*75 INFILE,OUTFILE,ERRFILE
        LOGICAL LROTATE
C------------------------------------------------------------------------------
        NEW_NSCAN=0 !avoid compilation warning
        NEW_NCHAN=0 !avoid compilation warning
        THISPROGRAM='rotate'
        CALL WELCOME('8-December-1996')
C------------------------------------------------------------------------------
        WRITE(*,100)'Work with error images (y/n) '
        CERR(1:1)=READC('n','yn')
C------------------------------------------------------------------------------
        WRITE(*,100)'Input file name'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        DO I=1,NSCAN
          READ(20) (A(J,I),J=1,NCHAN)
        END DO
        CLOSE(20)
        IF(CERR.EQ.'y')THEN
          WRITE(*,100)'Error file name '
          CALL GUESSEF(INFILE,ERRFILE)
          ERRFILE=INFILEX(30,ERRFILE,NSCAN,NCHAN,STWV,DISP,21,.TRUE.) !match
          DO I=1,NSCAN
            READ(30) (AERR(J,I),J=1,NCHAN)
          END DO
          CLOSE(30)
        END IF
C------------------------------------------------------------------------------
C determine whether the frame can be rotated or not
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
C
C------------------------------------------------------------------------------
10      WRITE(*,*)
        WRITE(*,101)'(1) Reverse x direction'
        WRITE(*,101)'(2) Reverse y direction'
        IF(LROTATE)THEN
          WRITE(*,101)'(3) Rotate +90 degrees'
          WRITE(*,101)'(4) Rotate -90 degrees'
        END IF
        WRITE(*,101)'(0) Exit & Save'
        WRITE(*,100)'Option '
        IF(LROTATE)THEN
          IOPC=READILIM('0',0,4)
        ELSE
          IOPC=READILIM('0',0,2)
        END IF
C------------------------------------------------------------------------------
        IF(IOPC.EQ.0)THEN
          GOTO 90
C------------------------------------------------------------------------------
        ELSEIF(IOPC.EQ.1)THEN
          WRITE(*,100)'Reversing...'
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
C------------------------------------------------------------------------------
        ELSEIF(IOPC.EQ.2)THEN
          WRITE(*,100)'Reversing...'
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
C------------------------------------------------------------------------------
        ELSEIF(IOPC.EQ.3)THEN
            WRITE(*,100)'Rotating...'
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
C------------------------------------------------------------------------------
        ELSEIF(IOPC.EQ.4)THEN
            WRITE(*,100)'Rotating...'
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
C------------------------------------------------------------------------------
        END IF
C------------------------------------------------------------------------------
C Update the original frames
        WRITE(*,100)'Updating...'
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
C------------------------------------------------------------------------------
90      WRITE(*,100)'Output file name'
        OUTFILE=OUTFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        DO I=1,NSCAN
          WRITE(20) (A(J,I),J=1,NCHAN)
        END DO
        CLOSE(20)
C
        IF(CERR.EQ.'y')THEN
          WRITE(*,100)'Output error file name '
          CALL GUESSEF(OUTFILE,ERRFILE)
          ERRFILE=OUTFILEX(30,ERRFILE,NSCAN,NCHAN,STWV,DISP,1,.TRUE.)
          DO I=1,NSCAN
            WRITE(30) (AERR(J,I),J=1,NCHAN)
          END DO
          CLOSE(30)
        END IF
C------------------------------------------------------------------------------
        STOP
100     FORMAT(A,$)
101     FORMAT(A)
        END
