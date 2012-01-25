C------------------------------------------------------------------------------
C Version 7-December-1996                                         file: irevx.f
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
C Program: irevx
C Classification: arithmetic & manipulations
C Description: Reverses an image (or spectrum) in the X-direction.
C
Comment
C
        PROGRAM IREVX
        IMPLICIT NONE
C Include NSMAX,NCMAX,NSCAN,NCHAN
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
C generic counters
        INTEGER I,J
C image matrix
        REAL S(NCMAX)
        REAL SS(NCMAX)
C file names
        CHARACTER*80 INFILE,OUTFILE
C----------------------------------------------------------------------
        THISPROGRAM='irevx'
        CALL WELCOME('7-December-1996')
C
        WRITE(*,100)'Input file name'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
C
        WRITE(*,100)'Outpuf file name'
        OUTFILE=OUTFILEX(15,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
C
        DO I=1,NSCAN
          READ(20)(SS(J),J=1,NCHAN)
          DO J=1,NCHAN
            S(J)=SS(NCHAN+1-J)
          END DO
          WRITE(15)(S(J),J=1,NCHAN)
        END DO
        CLOSE(20)
        CLOSE(15)
        STOP
100     FORMAT(A,$)
        END
