!------------------------------------------------------------------------------
! Version 20-November-2026                                  file: outfitshead.f
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
! Program: outfitshead
! Classification: miscellany
! Description: Outputs the whole header of a FITS file into a text file.
!
!Comment
!
! Este programa lee un fichero en formato FITS normal y escribe toda la
! cabecera en un fichero de texto 
!
        PROGRAM OUTFITSHEAD
        IMPLICIT NONE
!
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        CHARACTER*255 READC_FILE
        INTEGER TRUELEN
!
        INTEGER I,I0
        INTEGER IHEADER
        CHARACTER*255 OUTFILE
        CHARACTER*255 INFILE
        CHARACTER*2880 HDU
!------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
        OUTFILEX=OUTFILEX
!
        THISPROGRAM='outfitshead'
        CALL WELCOME('20-November-2026')
!
        WRITE(*,100,ADVANCE='NO')'FITS file name'
        INFILE=INFILEX(20,'@',0,0,0.,0.,4,.FALSE.)
!
        WRITE(*,100,ADVANCE='NO')'Output file name (header information)'
        OUTFILE=READC_FILE('@','@')
        OPEN(30,FILE=OUTFILE,STATUS='UNKNOWN',FORM='FORMATTED')
5       READ(30,'(A2880)',END=6) HDU
        GOTO 5
6       CONTINUE
        WRITE(30,150)
        WRITE(30,101)'This file is: '//INFILE(1:TRUELEN(INFILE))
        WRITE(30,151)
! Leemos el primer HDU (Header Data Unit)
        IHEADER=0
20      IHEADER=IHEADER+1
        WRITE(*,110)'Next header is #',IHEADER
        READ(20,'(A2880)',REC=IHEADER) HDU
        DO I=1,36
          I0=(I-1)*80
          WRITE(30,'(A80)') HDU(I0+1:I0+80)
!cc          WRITE(*,'(A80)') HDU(I0+1:I0+80)
          IF(HDU(I0+1:I0+3).EQ.'END') GOTO 22                     !fin de HDU's
        END DO
        GOTO 20
!
22      WRITE(*,110)'> No. of HDUs read: ',IHEADER
        CLOSE(20)
        CLOSE(30)
!------------------------------------------------------------------------------
        STOP
!
100     FORMAT(A)
101     FORMAT(A)
110     FORMAT(A,I6)
150     FORMAT(79('*'))
151     FORMAT(79('-'))
!
        END
