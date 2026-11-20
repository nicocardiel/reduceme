!------------------------------------------------------------------------------
! Version 20-November-2026                                       file: stplot.f
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
! Program: stplot
! Classification: graphic display
! Description: Plots the spectrophotometric spectra.
!
!Comment
!
! Programa para representar graficamente las tablas que contienen las curvas
! de las estrellas estandard calibradas (es decir, dos columnas, precedidas
! por un entero que es el numero de datos).
!
        PROGRAM STPLOT
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        REAL READF
!
        INTEGER NMAX
        PARAMETER(NMAX=4000)                           !Numero maximo de puntos
        INTEGER I,NPT
        REAL X(NMAX),Y(NMAX)
        REAL XMIN,XMAX,YMIN,YMAX,DX,DY
        REAL YOFFSET
        CHARACTER*1 COVER
        CHARACTER*255 INFILE
!------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
        OUTFILEX=OUTFILEX
!
        THISPROGRAM='stplot'
        CALL WELCOME('20-November-2026')
!------------------------------------------------------------------------------
        COVER='n'
10      WRITE(*,100,ADVANCE='NO')'Standard file name'
        INFILE=INFILEX(20,'@',0,0,0.,0.,3,.FALSE.)
        READ(20,*) NPT
        DO I=1,NPT
          READ(20,*) X(I),Y(I)
        END DO
        CLOSE(20)
        IF(COVER.EQ.'n')THEN
          CALL FINDMM(NPT,X,XMIN,XMAX)
          CALL FINDMM(NPT,Y,YMAX,YMIN)
          DX=XMAX-XMIN
          DY=YMAX-YMIN
          XMIN=XMIN-DX/50.
          XMAX=XMAX+DX/50.
          YMIN=YMIN-DY/50.
          YMAX=YMAX+DY/50.
          CALL PGBEGIN(0,'?',1,1)
          CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,0)
          CALL PGLABEL('wavelength','Flux',INFILE)
        END IF
        IF(COVER.EQ.'y')THEN
          WRITE(*,100,ADVANCE='NO')'Y-offset '
          YOFFSET=READF('0.0')
          DO I=1,NPT
            Y(I)=Y(I)+YOFFSET
          END DO
        END IF
        CALL PGBIN(NPT,X,Y,.TRUE.)
        WRITE(*,100,ADVANCE='NO')'Overplot (y/n) '
        COVER(1:1)=READC('n','yn')
        IF(COVER.EQ.'y') GOTO 10
        CALL PGEND
        STOP
100     FORMAT(A)
        END
