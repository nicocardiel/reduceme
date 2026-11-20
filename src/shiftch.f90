!------------------------------------------------------------------------------
! Version 20-November-2026                                      file: shiftch.f
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
! Program: shiftch
! Classification: wavelengths
! Description: Applies a constant shift (channels) to an image.
!
!Comment
!
! Traslada todos los espectros de una imagen un factor constante en la
! direccion X (l.d.o.)
!
        PROGRAM SHIFTCH
!
        IMPLICIT NONE
!
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        REAL READF
!
        INTEGER I,J
        REAL S(NCMAX),SS(NCMAX)
        REAL CHANSHIFT
        CHARACTER*1 CERR
        CHARACTER*255 INFILE,ERRINFILE,OUTFILE,ERROUTFILE
!------------------------------------------------------------------------------
        THISPROGRAM='shiftch'
        CALL WELCOME('20-November-2026')
!------------------------------------------------------------------------------
        WRITE(*,100,ADVANCE='NO')'Work with error images (y/n) '
        CERR(1:1)=READC('y','yn')
!
        WRITE(*,100,ADVANCE='NO')'Input file name'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
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
        WRITE(*,100,ADVANCE='NO')'Enter channel shift'
        CHANSHIFT=READF('@')
!
        WRITE(*,100,ADVANCE='NO')'Wait...'
        DO I=1,NSCAN
          READ(20) (S(J),J=1,NCHAN)
          CALL CHREBIN(CHANSHIFT,NCHAN,S,SS)
          WRITE(30) (SS(J),J=1,NCHAN)
          IF(CERR.EQ.'y')THEN
            READ(21) (S(J),J=1,NCHAN)
            CALL CHREBIN(CHANSHIFT,NCHAN,S,SS)
            WRITE(31) (SS(J),J=1,NCHAN)
          END IF
        END DO
        WRITE(*,101)'   ...OK!'
!------------------------------------------------------------------------------
        CLOSE(20)
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
