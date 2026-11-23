!------------------------------------------------------------------------------
! Version 24-November-2026                                   file: resample.f90
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
! Program: resample
! Classification: arithmetic & manipulations
! Description: Transforms an image with an initial STWV and DISP into another 
! image with a different STWV and DISP.
!
!Comment
!
! Este programa transforma una imagen con una determinada STWV y DISP a otra
! imagen con una nueva STWV y DISP.
!
        PROGRAM RESAMPLE
        USE RED_M
        USE FSPEC_M
        USE FUTILS_M
        USE REDLIB_M
        IMPLICIT NONE
!
        INCLUDE 'redlib.inc'
!
        INTEGER I,J
        INTEGER NCHAN2
        REAL, ALLOCATABLE :: A(:,:),ERR(:,:)
        REAL, ALLOCATABLE :: AUX(:,:)
        REAL STWV2,DISP2
        REAL S(NCMAX),SS(NCMAX),X(NCMAX),XX(NCMAX)
        REAL FACTOR
        CHARACTER*1 CERR,CFLUX
        CHARACTER*50 CDUMMY
        CHARACTER*255 INFILE,ERRFILE,OUTFILE
!------------------------------------------------------------------------------
        THISPROGRAM='resample'
        CALL WELCOME('24-November-2026')
! imagen de entrada (y errores)
        WRITE(*,100,ADVANCE='NO')'Work with error images (y/n) '
        CERR(1:1)=READC('n','yn')
        WRITE(*,100,ADVANCE='NO')'Preserve flux (y/n) '
        CFLUX(1:1)=READC('n','yn')
!
        WRITE(*,100,ADVANCE='NO')'Input file name'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        ! allocate with the actual image dimensions (static arrays
        ! were limited by NCMAX,NSMAX); initialised as before
        IF(ALLOCATED(A)) DEALLOCATE(A)
        IF(ALLOCATED(ERR)) DEALLOCATE(ERR)
        ALLOCATE(A(NCHAN,NSCAN),ERR(NCHAN,NSCAN))
        A=0.0
        ERR=0.0
        DO I=1,NSCAN
          READ(20) (A(J,I),J=1,NCHAN)
        END DO
        CLOSE(20)
        IF(CERR.EQ.'y')THEN
          CALL GUESSEF(INFILE,ERRFILE)
          WRITE(*,100,ADVANCE='NO')'Error file name '
          ERRFILE=INFILEX(21,ERRFILE,NSCAN,NCHAN,STWV,DISP,21,.TRUE.) !match
          DO I=1,NSCAN
            READ(21) (ERR(J,I),J=1,NCHAN)
          END DO
          CLOSE(21)
        END IF
!------------------------------------------------------------------------------
! nuevos parametros de la imagen de salida
        WRITE(CDUMMY,*) STWV
        WRITE(*,100,ADVANCE='NO')'New STWV '
        STWV2=READF(CDUMMY)
!
        WRITE(CDUMMY,*) DISP
        WRITE(*,100,ADVANCE='NO')'New DISP '
        DISP2=READF(CDUMMY)
!
        WRITE(CDUMMY,*) NCHAN
        WRITE(*,100,ADVANCE='NO')'New NCHAN '
        NCHAN2=READILIM(CDUMMY,1,NCMAX)
! the resampled spectra (NCHAN2 channels) overwrite A and ERR: enlarge
! them if needed, keeping their contents (and zeros elsewhere)
        IF(NCHAN2.GT.NCHAN)THEN
          ALLOCATE(AUX(NCHAN2,NSCAN))
          AUX=0.0
          AUX(1:NCHAN,:)=A
          CALL MOVE_ALLOC(AUX,A)
          ALLOCATE(AUX(NCHAN2,NSCAN))
          AUX=0.0
          AUX(1:NCHAN,:)=ERR
          CALL MOVE_ALLOC(AUX,ERR)
        END IF
!------------------------------------------------------------------------------
        IF(CFLUX.EQ.'y')THEN
          FACTOR=DISP2/DISP
        ELSE
          FACTOR=1.0
        END IF
!------------------------------------------------------------------------------
! creamos los nuevos espectros
        DO J=1,NCHAN
          X(J)=STWV+REAL(J-1)*DISP
        END DO
!
        DO I=1,NSCAN
          DO J=1,NCHAN
            S(J)=A(J,I)
          END DO
          CALL REBINING(X,S,NCHAN,XX,SS,NCHAN2,STWV2,DISP2)
          DO J=1,NCHAN2
            A(J,I)=SS(J)*FACTOR
          END DO
        END DO
!
        IF(CERR.EQ.'y')THEN
          DO I=1,NSCAN
            DO J=1,NCHAN
              S(J)=ERR(J,I)
            END DO
            CALL REBINING(X,S,NCHAN,XX,SS,NCHAN2,STWV2,DISP2)
            DO J=1,NCHAN2
              ERR(J,I)=SS(J)*FACTOR
            END DO
          END DO
        END IF
!------------------------------------------------------------------------------
! salvamos la nueva imagen
        WRITE(*,100,ADVANCE='NO')'Output file name'
        OUTFILE=OUTFILEX(30,'@',NSCAN,NCHAN2,STWV2,DISP2,1,.FALSE.)
        DO I=1,NSCAN
          WRITE(30) (A(J,I),J=1,NCHAN2)
        END DO
        CALL CLOSEFILEX(30)
        IF(CERR.EQ.'y')THEN
          CALL GUESSEF(OUTFILE,ERRFILE)
          WRITE(*,100,ADVANCE='NO')'Output error file name '
          OUTFILE=OUTFILEX(31,ERRFILE,NSCAN,NCHAN2,STWV2,DISP2,1,.TRUE.)
          DO I=1,NSCAN
            WRITE(31) (ERR(J,I),J=1,NCHAN2)
          END DO
          CALL CLOSEFILEX(31)
        END IF
!------------------------------------------------------------------------------
        STOP
100     FORMAT(A)
        END
