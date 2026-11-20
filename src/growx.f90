!------------------------------------------------------------------------------
! Version 20-November-2026                                      file: growx.f90
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
! Program: growx
! Classification: arithmetic & manipulations
! Description: Expands a single spectrum into an image.
!
!Comment
!
! Genera una imagen a partir de un corte espacial
!
! Para compilar
! f77pgp growx growx.f Lred.a Lfutils.a
!
        PROGRAM GROWX
        USE RED_M
        USE FUTILS_M
        USE REDLIB_M
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
!
        INTEGER I,J
        REAL S(NSMAX)
        REAL, ALLOCATABLE :: A(:,:)
        CHARACTER*255 INFILE,OUTFILE
!------------------------------------------------------------------------------
        THISPROGRAM='growx'
        CALL WELCOME('20-November-2026')
!
        WRITE(*,100,ADVANCE='NO')'Input file name (spatial direction)'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        IF(NCHAN.GT.1)THEN
          WRITE(*,101)'FATAL ERROR: this file is not a single'//&
     &     ' spatial direction.'
          CLOSE(20)
          STOP
        END IF
        DO I=1,NSCAN
          READ(20) S(I)
        END DO
        CLOSE(20)
!
!cc10      WRITE(*,100)'Number of channels in output '
        WRITE(*,100,ADVANCE='NO')'Number of channels in output '
        NCHAN=READILIM('@',1,NCMAX)
        ! allocate with the actual image dimensions (static arrays
        ! were limited by NCMAX,NSMAX); initialised as before
        IF(ALLOCATED(A)) DEALLOCATE(A)
        ALLOCATE(A(NCHAN,NSCAN))
        A=0.0
!
        DO I=1,NSCAN
          DO J=1,NCHAN
            A(J,I)=S(I)
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
