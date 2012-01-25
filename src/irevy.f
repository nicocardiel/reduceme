C------------------------------------------------------------------------------
C Version 7-December-1996                                         file: irevy.f
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
C Program: irevy
C Classification: arithmetic & manipulations
C Description: Reverse an image (or spatial direction) in the Y-direction.
C
Comment
C
        PROGRAM IREVY
        IMPLICIT NONE
C Include NSMAX,NCMAX,NSCAN,NCHAN
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
C generic counters
        INTEGER I,J
C image matrix
        REAL A(NCMAX,NSMAX)
        REAL SY(NSMAX)
C file names
        CHARACTER*80 INFILE,OUTFILE
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        THISPROGRAM='irevy'
        CALL WELCOME('7-December-1996')
C fichero de entrada
        WRITE(*,100)'Input file name'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        DO I=1,NSCAN
          READ(20)(A(J,I),J=1,NCHAN)
        END DO
        CLOSE(20)
C damos la vuelta a la imagen
        DO J=1,NCHAN
          DO I=1,NSCAN
            SY(I)=A(J,NSCAN-I+1)
          END DO
          DO I=1,NSCAN
            A(J,I)=SY(I)
          END DO
        END DO
C la salvamos
        WRITE(*,100)'Outpuf file name'
        OUTFILE=OUTFILEX(15,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        DO I=1,NSCAN
          WRITE(15)(A(J,I),J=1,NCHAN)
        END DO
        CLOSE(15)
C------------------------------------------------------------------------------
        STOP
100     FORMAT(A,$)
        END
