C------------------------------------------------------------------------------
C Version 6-December-1996                                         file: growy.f
C @ ncl & fjg
C------------------------------------------------------------------------------
Comment
C
C Program: growy
C Classification: arithmetic & manipulations
C Description: Expands a single spatial cross section into an image.
C
Comment
C
C Genera una imagen a partir de un espectro
C
        PROGRAM GROWY
        IMPLICIT NONE
C Include NSMAX,NCMAX,NSCAN,NCHAN
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INTEGER READILIM
C
        INTEGER I,J
        REAL S(NCMAX),A(NCMAX,NSMAX)
        CHARACTER*75 INFILE,OUTFILE
C------------------------------------------------------------------------------
        THISPROGRAM='growy'
        CALL WELCOME('6-December-1996')
C
        WRITE(*,100)'Input file name (spectrum)'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        IF(NSCAN.GT.1)THEN
          WRITE(*,101)'FATAL ERROR: this file is not a single '//
     +     'spectrum.'
          CLOSE(20)
          STOP
        END IF
        READ(20) (S(J),J=1,NCHAN)
        CLOSE(20)
C
10      WRITE(*,100)'Number of scans in output '
        NSCAN=READILIM('@',1,NSMAX)
        IF((NSCAN.LT.1).OR.(NSCAN.GT.NSMAX))THEN
          WRITE(*,101)'ERROR: number out of range. Try again.'
          GOTO 10
        END IF
C
        DO I=1,NSCAN
          DO J=1,NCHAN
            A(J,I)=S(J)
          END DO
        END DO
C
        WRITE(*,100)'Output file name'
        OUTFILE=OUTFILEX(30,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        DO I=1,NSCAN
          WRITE(30) (A(J,I),J=1,NCHAN)
        END DO
        CLOSE(30)
C
        STOP
100     FORMAT(A,$)
101     FORMAT(A)
        END
