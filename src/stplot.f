C------------------------------------------------------------------------------
C Version 8-December-1996                                        file: stplot.f
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
C Program: stplot
C Classification: graphic display
C Description: Plots the spectrophotometric spectra.
C
Comment
C
C Programa para representar graficamente las tablas que contienen las curvas
C de las estrellas estandard calibradas (es decir, dos columnas, precedidas
C por un entero que es el numero de datos).
C
        PROGRAM STPLOT
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        REAL READF
C
        INTEGER NMAX
        PARAMETER(NMAX=4000)                           !Numero maximo de puntos
        INTEGER I,NPT
        REAL X(NMAX),Y(NMAX)
        REAL XMIN,XMAX,YMIN,YMAX,DX,DY
        REAL YOFFSET
        CHARACTER*1 COVER
        CHARACTER*75 INFILE
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
        OUTFILEX=OUTFILEX
C
        THISPROGRAM='stplot'
        CALL WELCOME('8-December-1996')
C------------------------------------------------------------------------------
        COVER='n'
10      WRITE(*,100)'Standard file name'
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
          WRITE(*,100)'Y-offset '
          YOFFSET=READF('0.0')
          DO I=1,NPT
            Y(I)=Y(I)+YOFFSET
          END DO
        END IF
        CALL PGBIN(NPT,X,Y,.TRUE.)
        WRITE(*,100)'Overplot (y/n) '
        COVER(1:1)=READC('n','yn')
        IF(COVER.EQ.'y') GOTO 10
        CALL PGEND
        STOP
100     FORMAT(A,$)
        END
