!------------------------------------------------------------------------------
! Version 20-November-2026                                     file: randsc.f90
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
! Program: randsc
! Classification: error handling
! Description: Creates a fake image through bootstraping from an original 
! image and its associated error frame. In each pixel, the original
! signal is randomly modified using the correspoding error.
!
!Comment
!
        PROGRAM RANDSC
        USE RED_M
        USE FSPEC_M
        USE REDLIB_M
        IMPLICIT NONE
! Include NSMAX,NCMAX,NSCAN,NCHAN
        INCLUDE 'redlib.inc'
!
        REAL PI2
        PARAMETER (PI2=6.283185307)
!
        INTEGER I,J
        INTEGER NSEED
        REAL A(NCMAX,NSMAX),B(NCMAX,NSMAX)
        REAL ERR(NCMAX,NSMAX)
        REAL R1,R2,SXX
        CHARACTER*255 INFILE,ERRFILE,OUTFILE
!------------------------------------------------------------------------------
        THISPROGRAM='randsc'
        CALL WELCOME('20-November-2026')
! leemos ficheros
        WRITE(*,100,ADVANCE='NO')'Input file name'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        DO I=1,NSCAN
          READ(20) (A(J,I),J=1,NCHAN)
        END DO
        CLOSE(20)
        WRITE(*,100,ADVANCE='NO')'Input error file name '
        CALL GUESSEF(INFILE,ERRFILE)
        ERRFILE=INFILEX(21,ERRFILE,NSCAN,NCHAN,STWV,DISP,21,.TRUE.)    !match
        DO I=1,NSCAN
          READ(21) (ERR(J,I),J=1,NCHAN)
        END DO
        CLOSE(21)
!------------------------------------------------------------------------------
! calculamos semilla para generar los numeros aleatorios
        NSEED=-1
!------------------------------------------------------------------------------
        WRITE(*,100,ADVANCE='NO')'Thinking...'
!
        DO I=1,NSCAN
          DO J=1,NCHAN
            R1=RANRED(NSEED)
            R2=RANRED(NSEED)
            SXX=1.414213562*ERR(J,I)*SQRT(-1.*LOG(1.-R1))*COS(PI2*R2)
            B(J,I)=A(J,I)+SXX
          END DO
        END DO
!
        WRITE(*,101)'OK!'
!------------------------------------------------------------------------------
! salvamos fichero
        WRITE(*,100,ADVANCE='NO')'Output file name'
        OUTFILE=OUTFILEX(30,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        DO I=1,NSCAN
          WRITE(30) (B(J,I),J=1,NCHAN)
        END DO
        CLOSE(30)
!
        STOP
100     FORMAT(A)
101     FORMAT(A)
        END
