!------------------------------------------------------------------------------
! Version 24-November-2026                                    file: rvshift.f90
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
! Program: rvshift
! Classification: wavelengths
! Description: Applies a radial velocity shift to an image.
!
!Comment
!
! Traslada un espectro una cierta velocidad, manteniendolo en escala lineal
! con la misma STWV y DISP que el espectro de entrada.
!
        PROGRAM RVSHIFT
        USE RED_M
        USE FSPEC_M
        USE FUTILS_M
        USE REDLIB_M
!
        IMPLICIT NONE
!
        INCLUDE 'redlib.inc'
!
        INTEGER I,J
        REAL S(NCMAX),SS(NCMAX)
        REAL RADVEL
        CHARACTER*1 CERR
        CHARACTER*255 INFILE,ERRINFILE,OUTFILE,ERROUTFILE
!------------------------------------------------------------------------------
        THISPROGRAM='rvshift'
        CALL WELCOME('24-November-2026')
!------------------------------------------------------------------------------
        WRITE(*,100,ADVANCE='NO')'Work with error images (y/n) '
        CERR(1:1)=READC('y','yn')
!
        WRITE(*,100,ADVANCE='NO')'Input file name'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        IF((STWV.EQ.0.0).AND.(DISP.EQ.0.0))THEN
          WRITE(*,101)'FATAL ERROR: STWV=0.0 and DISP=0.0!'
          CLOSE(20)
          STOP
        END IF
!
        IF(CERR.EQ.'y')THEN
          WRITE(*,100,ADVANCE='NO')'Input error file name '
          CALL GUESSEF(INFILE,ERRINFILE)
          ERRINFILE=&
     &     INFILEX(21,ERRINFILE,NSCAN,NCHAN,STWV,DISP,21,.TRUE.) !........match
        END IF
!
        WRITE(*,100,ADVANCE='NO')'Output file name'
        OUTFILE=OUTFILEX(30,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        IF(CERR.EQ.'y')THEN
          WRITE(*,100,ADVANCE='NO')'Output error file name '
          CALL GUESSEF(OUTFILE,ERROUTFILE)
          ERROUTFILE=&
     &     OUTFILEX(31,ERROUTFILE,NSCAN,NCHAN,STWV,DISP,1,.TRUE.) !.......match
        END IF
        
!
        WRITE(*,100,ADVANCE='NO')'Radial velocity (km/sec)'
        RADVEL=READF('@')
!
        WRITE(*,100,ADVANCE='NO')'Wait...'
        DO I=1,NSCAN
          READ(20) (S(J),J=1,NCHAN)
          CALL RVREBIN(RADVEL,NCHAN,S,SS,STWV,DISP)
          WRITE(30) (SS(J),J=1,NCHAN)
          IF(CERR.EQ.'y')THEN
            READ(21) (S(J),J=1,NCHAN)
            CALL RVREBIN(RADVEL,NCHAN,S,SS,STWV,DISP)
            WRITE(31) (SS(J),J=1,NCHAN)
          END IF
        END DO
        WRITE(*,101)'   ...OK!'
!------------------------------------------------------------------------------
        CLOSE(20)
        CALL CLOSEFILEX(30)
        IF(CERR.EQ.'y')THEN
          CLOSE(21)
          CALL CLOSEFILEX(31)
        END IF
!------------------------------------------------------------------------------
        STOP
100     FORMAT(A)
101     FORMAT(A)
!
        END
