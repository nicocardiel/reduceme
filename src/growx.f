C------------------------------------------------------------------------------
C Version 6-December-1996                                         file: growx.f
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
C Program: growx
C Classification: arithmetic & manipulations
C Description: Expands a single spectrum into an image.
C
Comment
C
C Genera una imagen a partir de un corte espacial
C
C Para compilar
C f77pgp growx growx.f Lred.a Lfutils.a
C
        PROGRAM GROWX
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INTEGER READILIM
C
        INTEGER I,J
        REAL S(NSMAX),A(NCMAX,NSMAX)
        CHARACTER*75 INFILE,OUTFILE
C------------------------------------------------------------------------------
        THISPROGRAM='growx'
        CALL WELCOME('6-December-1996')
C
        WRITE(*,100)'Input file name (spatial direction)'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        IF(NCHAN.GT.1)THEN
          WRITE(*,101)'FATAL ERROR: this file is not a single'//
     +     ' spatial direction.'
          CLOSE(20)
          STOP
        END IF
        DO I=1,NSCAN
          READ(20) S(I)
        END DO
        CLOSE(20)
C
ccc10      WRITE(*,100)'Number of channels in output '
        WRITE(*,100)'Number of channels in output '
        NCHAN=READILIM('@',1,NCMAX)
C
        DO I=1,NSCAN
          DO J=1,NCHAN
            A(J,I)=S(I)
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
