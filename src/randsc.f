C------------------------------------------------------------------------------
C Version 07-September-2007                                      file: randsc.f
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
C Program: randsc
C Classification: error handling
C Description: Creates a fake image through bootstraping from an original 
C image and its associated error frame. In each pixel, the original
C signal is randomly modified using the correspoding error.
C
Comment
C
        PROGRAM RANDSC
        IMPLICIT NONE
C Include NSMAX,NCMAX,NSCAN,NCHAN
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
C
        REAL PI2
        PARAMETER (PI2=6.283185307)
C
        INTEGER I,J
        INTEGER NSEED
        REAL A(NCMAX,NSMAX),B(NCMAX,NSMAX)
        REAL ERR(NCMAX,NSMAX)
        REAL R1,R2,SXX
        REAL RANRED
        CHARACTER*75 INFILE,ERRFILE,OUTFILE
C------------------------------------------------------------------------------
        THISPROGRAM='randsc'
        CALL WELCOME('07-September-2007')
C leemos ficheros
        WRITE(*,100)'Input file name'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        DO I=1,NSCAN
          READ(20) (A(J,I),J=1,NCHAN)
        END DO
        CLOSE(20)
        WRITE(*,100)'Input error file name '
        CALL GUESSEF(INFILE,ERRFILE)
        ERRFILE=INFILEX(21,ERRFILE,NSCAN,NCHAN,STWV,DISP,21,.TRUE.)    !match
        DO I=1,NSCAN
          READ(21) (ERR(J,I),J=1,NCHAN)
        END DO
        CLOSE(21)
C------------------------------------------------------------------------------
c calculamos semilla para generar los numeros aleatorios
        NSEED=-1
C------------------------------------------------------------------------------
        WRITE(*,100)'Thinking...'
C
        DO I=1,NSCAN
          DO J=1,NCHAN
            R1=RANRED(NSEED)
            R2=RANRED(NSEED)
            SXX=1.414213562*ERR(J,I)*SQRT(-1.*LOG(1.-R1))*COS(PI2*R2)
            B(J,I)=A(J,I)+SXX
          END DO
        END DO
C
        WRITE(*,101)'OK!'
C------------------------------------------------------------------------------
C salvamos fichero
        WRITE(*,100)'Output file name'
        OUTFILE=OUTFILEX(30,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        DO I=1,NSCAN
          WRITE(30) (B(J,I),J=1,NCHAN)
        END DO
        CLOSE(30)
C
        STOP
100     FORMAT(A,$)
101     FORMAT(A)
        END
