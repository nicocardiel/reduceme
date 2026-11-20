!------------------------------------------------------------------------------
! Version 20-November-2026                                      file: irevy.f90
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
! Program: irevy
! Classification: arithmetic & manipulations
! Description: Reverse an image (or spatial direction) in the Y-direction.
!
!Comment
!
        PROGRAM IREVY
        USE RED_M
        USE FUTILS_M
        USE REDLIB_M
        IMPLICIT NONE
! Include NSMAX,NCMAX,NSCAN,NCHAN
        INCLUDE 'redlib.inc'
! generic counters
        INTEGER I,J
! image matrix
        REAL A(NCMAX,NSMAX)
        REAL SY(NSMAX)
! file names
        CHARACTER*255 INFILE,OUTFILE
!------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
!
        THISPROGRAM='irevy'
        CALL WELCOME('20-November-2026')
! fichero de entrada
        WRITE(*,100,ADVANCE='NO')'Input file name'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        DO I=1,NSCAN
          READ(20)(A(J,I),J=1,NCHAN)
        END DO
        CLOSE(20)
! damos la vuelta a la imagen
        DO J=1,NCHAN
          DO I=1,NSCAN
            SY(I)=A(J,NSCAN-I+1)
          END DO
          DO I=1,NSCAN
            A(J,I)=SY(I)
          END DO
        END DO
! la salvamos
        WRITE(*,100,ADVANCE='NO')'Outpuf file name'
        OUTFILE=OUTFILEX(15,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        DO I=1,NSCAN
          WRITE(15)(A(J,I),J=1,NCHAN)
        END DO
        CLOSE(15)
!------------------------------------------------------------------------------
        STOP
100     FORMAT(A)
        END
