C------------------------------------------------------------------------------
C Version 7-December-1996                                   file: outfitshead.f
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
C Program: outfitshead
C Classification: miscellany
C Description: Outputs the whole header of a FITS file into a text file.
C
Comment
C
C Este programa lee un fichero en formato FITS normal y escribe toda la
C cabecera en un fichero de texto 
C
        PROGRAM OUTFITSHEAD
        IMPLICIT NONE
C
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER TRUELEN
C
        INTEGER I,I0
        INTEGER IHEADER
        CHARACTER*75 INFILE,OUTFILE
        CHARACTER*2880 HDU
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
        OUTFILEX=OUTFILEX
C
        THISPROGRAM='outfitshead'
        CALL WELCOME('7-December-1996')
C
        WRITE(*,100)'FITS file name'
        INFILE=INFILEX(20,'@',0,0,0.,0.,4,.FALSE.)
C
        WRITE(*,100)'Output file name (header information)'
        OUTFILE(1:75)=READC('@','@')
        OPEN(30,FILE=OUTFILE,STATUS='UNKNOWN',FORM='FORMATTED')
5       READ(30,'(A2880)',END=6) HDU
        GOTO 5
6       CONTINUE
        WRITE(30,150)
        WRITE(30,101)'This file is: '//INFILE(1:TRUELEN(INFILE))
        WRITE(30,151)
C Leemos el primer HDU (Header Data Unit)
        IHEADER=0
20      IHEADER=IHEADER+1
        WRITE(*,110)'Next header is #',IHEADER
        READ(20,'(A2880)',REC=IHEADER) HDU
        DO I=1,36
          I0=(I-1)*80
          WRITE(30,'(A80)') HDU(I0+1:I0+80)
ccc          WRITE(*,'(A80)') HDU(I0+1:I0+80)
          IF(HDU(I0+1:I0+3).EQ.'END') GOTO 22                     !fin de HDU's
        END DO
        GOTO 20
C
22      WRITE(*,110)'> No. of HDUs read: ',IHEADER
        CLOSE(20)
        CLOSE(30)
C------------------------------------------------------------------------------
        STOP
C
100     FORMAT(A,$)
101     FORMAT(A)
110     FORMAT(A,I6)
150     FORMAT(79('*'))
151     FORMAT(79('-'))
C
        END
