!------------------------------------------------------------------------------
! Version 24-November-2026                                     file: prfcal.f90
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
! Program: prfcal
! Classification: flux calibration
! Description: Creates a mean flux calibration curve from individual curves, 
! and generates an image containing the mean curve (as the first spectrum) and 
! the individual ones (in successive spectra).
!
!Comment
!
! Spectra (one in each column) shared by the main program and PLOT
! (they were in COMMON /BLKPLOT2/)
        MODULE PRFCAL_DATA_M
        IMPLICIT NONE
        REAL, ALLOCATABLE :: S(:,:)
        END MODULE PRFCAL_DATA_M
!
! Calcula la media de varias curvas respuesta para la calibracion en
! flujo, y genera una "imagen" que contiene, como primer espectro la
! curva promediada, y luego todas las curvas individuales. El programa
! que mide los indices lee este fichero y calcula el error debido a calibracion
! en flujo mediante el calculo de los indices con todas las curvas.
! Antes de calcular el promedio, el programa renormaliza todas las curvas
! a uno para obtener un promedio que tenga sentido.
!
        PROGRAM PRFCAL
        USE RED_M
        USE FUTILS_M
        USE REDLIB_M
        USE PRFCAL_DATA_M
        IMPLICIT NONE
! Include NSMAX,NCMAX,NSCAN,NCHAN
        INCLUDE 'redlib.inc'
!
        INTEGER I,J,L
        INTEGER NF
        INTEGER NC1,NC2
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        REAL SS(NCMAX)
        REAL X(NCMAX)
        REAL FACTOR(NSMAX)
        CHARACTER*50 CDUMMY
        CHARACTER*255 INFILE,OUTFILE
        LOGICAL LCOLOR(MAX_ID_RED)
!
        COMMON/BLKPLOT1/NF,NCHAN
        COMMON/BLKPLOT2/SS,X
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
!------------------------------------------------------------------------------
        THISPROGRAM='prfcal'
        CALL WELCOME('24-November-2026')
!------------------------------------------------------------------------------
        CALL PIDEGTER(NTERM,IDN,LCOLOR)
!------------------------------------------------------------------------------
        WRITE(*,100,ADVANCE='NO')'No. of spectra '
        NF=READILIM('@',1,NSMAX)
!
        DO I=1,NF
          WRITE(*,'(A,I2,A)',ADVANCE='NO')'Spectrum #',I,' file name'
          IF(I.EQ.1)THEN
            INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
! allocate with the actual dimensions (static array limited by
! NCMAX,NSMAX); initialised as before
            ALLOCATE(S(NCHAN,NF))
            S=0.0
          ELSE
            INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,21,.FALSE.)  !match
          END IF
          READ(20) (S(J,I),J=1,NCHAN)
          CLOSE(20)
        END DO
!
        DO J=1,NCHAN
          X(J)=REAL(J)
        END DO
!------------------------------------------------------------------------------
! plot different spectra
        CALL PLOT
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGLABEL(CHAR(32),CHAR(32),'Before renormalization')
        END DO
!------------------------------------------------------------------------------
10      WRITE(*,100,ADVANCE='NO')&
     & 'Channel region to renormalized spectra '
        WRITE(CDUMMY,'(A,I10)')'1,',NCHAN
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        CALL READ2I(CDUMMY(1:L),NC1,NC2)
        IF((NC1.LT.1).OR.(NC2.GT.NCHAN).OR.(NC1.GT.NC2))THEN
          WRITE(*,101)'ERROR: invalid entry. Try again.'
          GOTO 10
        END IF
!
        DO I=1,NF
          FACTOR(I)=0.
          DO J=NC1,NC2
            FACTOR(I)=FACTOR(I)+S(J,I)
          END DO
          FACTOR(I)=FACTOR(I)/REAL(NC2-NC1+1)
          WRITE(*,'(A,I3,A)',ADVANCE='NO')'Image #',I,'   --> Factor: '
          WRITE(*,*) FACTOR(I)
        END DO
!
        DO I=1,NF
          DO J=1,NCHAN
            S(J,I)=S(J,I)/FACTOR(I)
          END DO
        END DO
!
        CALL PLOT
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGLABEL(CHAR(32),CHAR(32),'After renormalization')
        END DO
!------------------------------------------------------------------------------
! mean spectrum
        DO J=1,NCHAN
          SS(J)=0.
        END DO
        DO I=1,NF
          DO J=1,NCHAN
            SS(J)=SS(J)+S(J,I)
          END DO
        END DO
        DO J=1,NCHAN
          SS(J)=SS(J)/REAL(NF)
        END DO
!
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGSLW(4)
          CALL PGBIN(NCHAN,X,SS,.TRUE.)
          CALL PGSLW(1)
        END DO
!
        OBJECT='[from prfcal]'
        WRITE(*,100,ADVANCE='NO')&
     &     'Output file name for mean+individual curves'
        OUTFILE=OUTFILEX(30,'@',NF+1,NCHAN,STWV,DISP,1,.FALSE.)
        WRITE(30) (SS(J),J=1,NCHAN)
        DO I=1,NF
          WRITE(30) (S(J,I),J=1,NCHAN)
        END DO
        CALL CLOSEFILEX(30)
!------------------------------------------------------------------------------
        WRITE(*,100,ADVANCE='NO')'Output file name for mean curve'
        OUTFILE=OUTFILEX(30,'@',1,NCHAN,STWV,DISP,1,.FALSE.)
        WRITE(30) (SS(J),J=1,NCHAN)
        CALL CLOSEFILEX(30)
!------------------------------------------------------------------------------
! create file with differences relative to the mean spectrum
        DO I=1,NF
          DO J=1,NCHAN
            S(J,I)=S(J,I)/SS(J)
          END DO
        END DO
!
        CALL PLOT
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGSLW(4)
          CALL PGMOVE(1.0,1.0)
          CALL PGDRAW(REAL(NCHAN),1.0)
          CALL PGSLW(1)
          CALL PGLABEL(CHAR(32),CHAR(32),&
     &     'Differences relative to mean spectrum')
        END DO
!
        WRITE(*,100,ADVANCE='NO') 'Output file name for differences'
        OUTFILE=OUTFILEX(30,'@',NF,NCHAN,STWV,DISP,1,.FALSE.)
        DO I=1,NF
          WRITE(30) (S(J,I),J=1,NCHAN)
        END DO
        CALL CLOSEFILEX(30)
!------------------------------------------------------------------------------
        CALL PGEND
!
        STOP
100     FORMAT(A)
101     FORMAT(A)
        END
!
!******************************************************************************
!
        SUBROUTINE PLOT
        USE RED_M
        USE FUTILS_M
        USE REDLIB_M
        USE PRFCAL_DATA_M
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
!
        INTEGER I,J
        INTEGER NF
        INTEGER NCOLOR
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        REAL SS(NCMAX),X(NCMAX)
        REAL XMAX,XMIN,YMAX,YMIN,DX,DY
        CHARACTER*1 CYLIMIT
        CHARACTER*50 CDUMMY
        LOGICAL LCOLOR(MAX_ID_RED)
!
        COMMON/BLKPLOT1/NF,NCHAN
        COMMON/BLKPLOT2/SS,X
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
!------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
!
        XMIN=1.
        XMAX=REAL(NCHAN)
        DX=XMAX-XMIN
        XMIN=XMIN-DX/30.
        XMAX=XMAX+DX/30.
        YMIN=S(1,1)
        YMAX=YMIN
        DO I=1,NF
          DO J=1,NCHAN
            IF(S(J,I).LT.YMIN) YMIN=S(J,I)
            IF(S(J,I).GT.YMAX) YMAX=S(J,I)
          END DO
        END DO
        DY=YMAX-YMIN
        YMIN=YMIN-DY/20.
        YMAX=YMAX+DY/20.
!
10      DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,0)
          CALL PGIDEN_RED
          CALL PGLABEL('channel','response',CHAR(32))
        END DO
        DO I=1,NF
          NCOLOR=MOD(I,12)+2
          DO J=1,NCHAN
            SS(J)=S(J,I)
          END DO
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(NCOLOR)
            CALL PGBIN(NCHAN,X,SS,.TRUE.)
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END DO
        END DO
!
        WRITE(*,100,ADVANCE='NO')'Change Y-limits (y/n) '
        CYLIMIT(1:1)=READC('n','yn')
        IF(CYLIMIT.EQ.'y')THEN
          WRITE(CDUMMY,*)YMIN
          WRITE(*,100,ADVANCE='NO')'Ymin '
          YMIN=READF(CDUMMY)
          WRITE(CDUMMY,*)YMAX
          WRITE(*,100,ADVANCE='NO')'Ymax '
          YMAX=READF(CDUMMY)
          GOTO 10
        END IF
!
100     FORMAT(A)
        END
