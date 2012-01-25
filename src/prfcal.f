C------------------------------------------------------------------------------
C Version 23-March-2004                                          file: prfcal.f
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
C Program: prfcal
C Classification: flux calibration
C Description: Creates a mean flux calibration curve from individual curves, 
C and generates an image containing the mean curve (as the first spectrum) and 
C the individual ones (in successive spectra).
C
Comment
C
C Calcula la media de varias curvas respuesta para la calibracion en
C flujo, y genera una "imagen" que contiene, como primer espectro la
C curva promediada, y luego todas las curvas individuales. El programa
C que mide los indices lee este fichero y calcula el error debido a calibracion
C en flujo mediante el calculo de los indices con todas las curvas.
C Antes de calcular el promedio, el programa renormaliza todas las curvas
C a uno para obtener un promedio que tenga sentido.
C
        PROGRAM PRFCAL
        IMPLICIT NONE
C Include NSMAX,NCMAX,NSCAN,NCHAN
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INTEGER READILIM
C
        INTEGER I,J,L
        INTEGER NF
        INTEGER NC1,NC2
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        REAL S(NCMAX,NSMAX),SS(NCMAX)
        REAL X(NCMAX)
        REAL FACTOR(NSMAX)
        CHARACTER*50 CDUMMY
        CHARACTER*75 INFILE,OUTFILE
        LOGICAL LCOLOR(MAX_ID_RED)
C
        COMMON/BLKPLOT1/NF,NCHAN
        COMMON/BLKPLOT2/S,SS,X
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
C------------------------------------------------------------------------------
        THISPROGRAM='prfcal'
        CALL WELCOME('23-March-2004')
C------------------------------------------------------------------------------
        CALL PIDEGTER(NTERM,IDN,LCOLOR)
C------------------------------------------------------------------------------
        WRITE(*,100)'No. of spectra '
        NF=READILIM('@',1,NSMAX)
C
        DO I=1,NF
          WRITE(*,'(A,I2,A,$)')'Spectrum #',I,' file name'
          IF(I.EQ.1)THEN
            INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
          ELSE
            INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,21,.FALSE.)  !match
          END IF
          READ(20) (S(J,I),J=1,NCHAN)
          CLOSE(20)
        END DO
C
        DO J=1,NCHAN
          X(J)=REAL(J)
        END DO
C------------------------------------------------------------------------------
C plot different spectra
        CALL PLOT
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGLABEL(CHAR(32),CHAR(32),'Before renormalization')
        END DO
C------------------------------------------------------------------------------
10      WRITE(*,100)'Channel region to renormalized spectra '
        WRITE(CDUMMY,'(A,I10)')'1,',NCHAN
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        CALL READ2I(CDUMMY(1:L),NC1,NC2)
        IF((NC1.LT.1).OR.(NC2.GT.NCHAN).OR.(NC1.GT.NC2))THEN
          WRITE(*,101)'ERROR: invalid entry. Try again.'
          GOTO 10
        END IF
C
        DO I=1,NF
          FACTOR(I)=0.
          DO J=NC1,NC2
            FACTOR(I)=FACTOR(I)+S(J,I)
          END DO
          FACTOR(I)=FACTOR(I)/REAL(NC2-NC1+1)
          WRITE(*,'(A,I3,A,$)')'Image #',I,'   --> Factor: '
          WRITE(*,*) FACTOR(I)
        END DO
C
        DO I=1,NF
          DO J=1,NCHAN
            S(J,I)=S(J,I)/FACTOR(I)
          END DO
        END DO
C
        CALL PLOT
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGLABEL(CHAR(32),CHAR(32),'After renormalization')
        END DO
C------------------------------------------------------------------------------
C mean spectrum
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
C
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGSLW(4)
          CALL PGBIN(NCHAN,X,SS,.TRUE.)
          CALL PGSLW(1)
        END DO
C
        OBJECT='[from prfcal]'
        WRITE(*,100)'Output file name for mean+individual curves'
        OUTFILE=OUTFILEX(30,'@',NF+1,NCHAN,STWV,DISP,1,.FALSE.)
        WRITE(30) (SS(J),J=1,NCHAN)
        DO I=1,NF
          WRITE(30) (S(J,I),J=1,NCHAN)
        END DO
        CLOSE(30)
C------------------------------------------------------------------------------
        WRITE(*,100)'Output file name for mean curve'
        OUTFILE=OUTFILEX(30,'@',1,NCHAN,STWV,DISP,1,.FALSE.)
        WRITE(30) (SS(J),J=1,NCHAN)
        CLOSE(30)
C------------------------------------------------------------------------------
C create file with differences relative to the mean spectrum
        DO I=1,NF
          DO J=1,NCHAN
            S(J,I)=S(J,I)/SS(J)
          END DO
        END DO
C
        CALL PLOT
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGSLW(4)
          CALL PGMOVE(1.0,1.0)
          CALL PGDRAW(REAL(NCHAN),1.0)
          CALL PGSLW(1)
          CALL PGLABEL(CHAR(32),CHAR(32),
     +     'Differences relative to mean spectrum')
        END DO
C
        WRITE(*,100) 'Output file name for differences'
        OUTFILE=OUTFILEX(30,'@',NF,NCHAN,STWV,DISP,1,.FALSE.)
        DO I=1,NF
          WRITE(30) (S(J,I),J=1,NCHAN)
        END DO
        CLOSE(30)
C------------------------------------------------------------------------------
        CALL PGEND
C
        STOP
100     FORMAT(A,$)
101     FORMAT(A)
        END
C
C******************************************************************************
C
        SUBROUTINE PLOT
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'futils.inc'
        REAL READF
C
        INTEGER I,J
        INTEGER NF
        INTEGER NCOLOR
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        REAL S(NCMAX,NSMAX),SS(NCMAX),X(NCMAX)
        REAL XMAX,XMIN,YMAX,YMIN,DX,DY
        CHARACTER*1 CYLIMIT
        CHARACTER*50 CDUMMY
        LOGICAL LCOLOR(MAX_ID_RED)
C
        COMMON/BLKPLOT1/NF,NCHAN
        COMMON/BLKPLOT2/S,SS,X
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
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
C
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
C
        WRITE(*,100)'Change Y-limits (y/n) '
        CYLIMIT(1:1)=READC('n','yn')
        IF(CYLIMIT.EQ.'y')THEN
          WRITE(CDUMMY,*)YMIN
          WRITE(*,100)'Ymin '
          YMIN=READF(CDUMMY)
          WRITE(CDUMMY,*)YMAX
          WRITE(*,100)'Ymax '
          YMAX=READF(CDUMMY)
          GOTO 10
        END IF
C
100     FORMAT(A,$)
        END
