!------------------------------------------------------------------------------
! Version 24-November-2026                                      file: irevx.f90
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
! Program: irevx
! Classification: arithmetic & manipulations
! Description: Reverses an image (or spectrum) in the X-direction.
!
!Comment
!
        PROGRAM IREVX
        USE RED_M
        USE REDLIB_M
        IMPLICIT NONE
! Include NSMAX,NCMAX,NSCAN,NCHAN
        INCLUDE 'redlib.inc'
! generic counters
        INTEGER I,J
! image matrix
        REAL S(NCMAX)
        REAL SS(NCMAX)
! file names
        CHARACTER*255 INFILE,OUTFILE
!----------------------------------------------------------------------
        THISPROGRAM='irevx'
        CALL WELCOME('24-November-2026')
!
        WRITE(*,100,ADVANCE='NO')'Input file name'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
!
        WRITE(*,100,ADVANCE='NO')'Outpuf file name'
        OUTFILE=OUTFILEX(15,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
!
        DO I=1,NSCAN
          READ(20)(SS(J),J=1,NCHAN)
          DO J=1,NCHAN
            S(J)=SS(NCHAN+1-J)
          END DO
          WRITE(15)(S(J),J=1,NCHAN)
        END DO
        CLOSE(20)
        CALL CLOSEFILEX(15)
        STOP
100     FORMAT(A)
        END
