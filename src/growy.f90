!------------------------------------------------------------------------------
! Version 20-November-2026                                      file: growy.f90
! @ ncl & fjg
!------------------------------------------------------------------------------
!Comment
!
! Program: growy
! Classification: arithmetic & manipulations
! Description: Expands a single spatial cross section into an image.
!
!Comment
!
! Genera una imagen a partir de un espectro
!
        PROGRAM GROWY
        IMPLICIT NONE
! Include NSMAX,NCMAX,NSCAN,NCHAN
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INTEGER READILIM
!
        INTEGER I,J
        REAL S(NCMAX),A(NCMAX,NSMAX)
        CHARACTER*255 INFILE,OUTFILE
!------------------------------------------------------------------------------
        THISPROGRAM='growy'
        CALL WELCOME('20-November-2026')
!
        WRITE(*,100,ADVANCE='NO')'Input file name (spectrum)'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        IF(NSCAN.GT.1)THEN
          WRITE(*,101)'FATAL ERROR: this file is not a single '//&
     &     'spectrum.'
          CLOSE(20)
          STOP
        END IF
        READ(20) (S(J),J=1,NCHAN)
        CLOSE(20)
!
10      WRITE(*,100,ADVANCE='NO')'Number of scans in output '
        NSCAN=READILIM('@',1,NSMAX)
        IF((NSCAN.LT.1).OR.(NSCAN.GT.NSMAX))THEN
          WRITE(*,101)'ERROR: number out of range. Try again.'
          GOTO 10
        END IF
!
        DO I=1,NSCAN
          DO J=1,NCHAN
            A(J,I)=S(J)
          END DO
        END DO
!
        WRITE(*,100,ADVANCE='NO')'Output file name'
        OUTFILE=OUTFILEX(30,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        DO I=1,NSCAN
          WRITE(30) (A(J,I),J=1,NCHAN)
        END DO
        CLOSE(30)
!
        STOP
100     FORMAT(A)
101     FORMAT(A)
        END
