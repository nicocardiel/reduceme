!------------------------------------------------------------------------------
! Version 20-November-2026                                   file: rvshiftrot.f
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
! Program: rvshiftrot
! Classification: wavelengths
! Description: Applies radial velocity shifts to the different scans of
! an image. This program can be used to correct 2D spectroscopic images
! from rotation curves. It is similar to rvshift, but using different
! radial velocities for each scan (these velocities are read from an
! ASCII file).
!
!Comment
!
        PROGRAM RVSHIFTROT
!
        IMPLICIT NONE
!
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
!
        INTEGER I,J
        REAL S(NCMAX),SS(NCMAX)
        REAL RADVEL
        CHARACTER*1 CERR,CSIGN
        CHARACTER*255 INFILE,ERRINFILE,OUTFILE,ERROUTFILE
        CHARACTER*255 RADVELFILE
!------------------------------------------------------------------------------
        THISPROGRAM='rvshiftrot'
        CALL WELCOME('20-November-2026')
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
        WRITE(*,100,ADVANCE='NO')&
     &     'File name with radial velocity for each scan '
        RADVELFILE=INFILEX(22,'@',0,0,.0,.0,3,.FALSE.)
        WRITE(*,100,ADVANCE='NO')'Change sign of velocities (y/n) '
        CSIGN(1:1)=READC('n','yn')
!
        WRITE(*,100,ADVANCE='NO')'Wait...'
        DO I=1,NSCAN
          READ(20) (S(J),J=1,NCHAN)
          READ(22,*) RADVEL
          IF(CSIGN.EQ.'y') RADVEL=-RADVEL
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
        CLOSE(22)
        CLOSE(30)
        IF(CERR.EQ.'y')THEN
          CLOSE(21)
          CLOSE(31)
        END IF
!------------------------------------------------------------------------------
        STOP
100     FORMAT(A)
101     FORMAT(A)
!
        END
