C------------------------------------------------------------------------------
C Version 6-December-1996                                       file: fitstex.f
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
C Program: fitstex
C Classification: miscellany
C Description: Reads the output from fitshead and generates an output LaTeX 
C file with a table.
C
Comment
C
C Este programa lee un fichero en formato ASCII generado con el programa
C fitshead, y produce un fichero en LaTeX con una tabla.
C
        PROGRAM FITSTEX
        IMPLICIT NONE
C
        INCLUDE 'futils.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'redlib.inc'
        INTEGER TRUELEN
        INTEGER READILIM
C
        INTEGER I,L,LL
        INTEGER NSKIP,NAXIS
        CHARACTER*1 CFORMAT
ccc        CHARACTER*10 FITSFILE
ccc        CHARACTER*20 OBJECT
        CHARACTER*75 INFILE,OUTFILE
        CHARACTER*255 CLINEA,CCOMMENT
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        THISPROGRAM='fitstex'
        CALL WELCOME('6-December-1996')
C
        WRITE(*,100)'Input file (output from fitshead)'
        INFILE=INFILEX(20,'@',0,0,0.,0.,3,.FALSE.)
        WRITE(*,100)'Format of previous file is [l]ong or [s]hort (l/s)'
        CFORMAT(1:1)=READC('@','ls')
        WRITE(*,100)'NAXIS '
        NAXIS=READILIM('2',1,3)
        WRITE(*,100)'Number of lines to be skipped '
        NSKIP=READILIM('1',0,9999)
        WRITE(*,100)'Ouput file (LaTeX)'
        OUTFILE=OUTFILEX(30,'@',0,0,0.,0.,3,.FALSE.)
        WRITE(*,100)'Comment'
        CCOMMENT=READC('@','@')
C------------------------------------------------------------------------------
C Escribimos cabecera del fichero LaTeX
        WRITE(30,101)'\\documentstyle[12pt]{article}'
        WRITE(30,101)'\\topmargin -0.5in'
        WRITE(30,101)'\\oddsidemargin -0.2in'
        WRITE(30,101)'\\textwidth   7.27in'
        WRITE(30,101)'\\textheight 11.15in'
        WRITE(30,101)'\\parindent 0mm'
        WRITE(30,101)'\\pagestyle{headings}'
        WRITE(30,101)
        WRITE(30,101)'\\begin{document}'
        WRITE(30,101)'{\\Large'
C------------------------------------------------------------------------------
C Saltamos lineas a ignorar
        IF(NSKIP.GE.1)THEN
          DO I=1,NSKIP
            READ(20,101) CLINEA
          END DO
        END IF
C------------------------------------------------------------------------------
        I=0
C Leemos nombre del fichero FITS y nombre del objeto
20      IF(CFORMAT.EQ.'l')THEN
          IF(NAXIS.EQ.1)THEN
            READ(20,'(A10,10X,1(6X),A20)',END=70)FITSFILE,OBJECT
          ELSEIF(NAXIS.EQ.2)THEN
            READ(20,'(A10,10X,2(6X),A20)',END=70)FITSFILE,OBJECT
          ELSEIF(NAXIS.EQ.3)THEN
            READ(20,'(A10,10X,3(6X),A20)',END=70)FITSFILE,OBJECT
          ELSE
            WRITE(*,100)'FATAL ERROR: invalid NAXIS.'
            CLOSE(20)
            CLOSE(30)
            STOP
          END IF
        ELSE
          IF(NAXIS.EQ.1)THEN
            READ(20,'(A10,1(6X),A20)',END=70)FITSFILE,OBJECT
          ELSEIF(NAXIS.EQ.2)THEN
            READ(20,'(A10,2(6X),A20)',END=70)FITSFILE,OBJECT
          ELSEIF(NAXIS.EQ.3)THEN
            READ(20,'(A10,3(6X),A20)',END=70)FITSFILE,OBJECT
          ELSE
            WRITE(*,100)'FATAL ERROR: invalid NAXIS.'
            CLOSE(20)
            CLOSE(30)
            STOP
          END IF
        END IF
        I=I+1
        IF(I.EQ.31)THEN
          WRITE(30,101)'\\end{tabular}'
          WRITE(30,101)
          WRITE(30,101)'\\newpage'
          I=1
        END IF
        IF(I.EQ.1)THEN
          WRITE(30,101)
          WRITE(30,101)'\\begin{center}'
          WRITE(30,101)'{\\normalsize'
          WRITE(30,101)CCOMMENT
          WRITE(30,101)'} %end of size'
          WRITE(30,101)'\\end{center}'
          WRITE(30,101)'\\begin{tabular}{|l|l|c|c|c|c|c|} \\hline'
          WRITE(30,101)
     +     '\\makebox[25mm][c]{\\large FITS file} '//
     +     '& \\makebox[30mm][c]{\\large Object} '//
     +     '& \\makebox[20mm][c]{}'//
     +     '& \\makebox[20mm][c]{}'//
     +     '& \\makebox[20mm][c]{}'//
     +     '& \\makebox[20mm][c]{}'//
     +     '& \\makebox[20mm][c]{} \\\\ \\hline'
        END IF
C------------------------------------------------------------------------------
        WRITE(30,100) '{\\normalsize '//FITSFILE//'}'//' & '//
     +   '{\\normalsize '
        LL=TRUELEN(OBJECT)
        IF(LL.GT.0)THEN
          DO L=1,LL
            IF(OBJECT(L:L).EQ.'$')THEN
              WRITE(30,100)'\\$'
            ELSEIF(OBJECT(L:L).EQ.'&')THEN
              WRITE(30,100)'\\&'
            ELSEIF(OBJECT(L:L).EQ.'%')THEN
              WRITE(30,100)'\\%'
            ELSEIF(OBJECT(L:L).EQ.'#')THEN
              WRITE(30,100)'\\#'
            ELSEIF(OBJECT(L:L).EQ.'_')THEN
              WRITE(30,100)'\\_'
            ELSEIF(OBJECT(L:L).EQ.'{')THEN
              WRITE(30,100)'\\{'
            ELSEIF(OBJECT(L:L).EQ.'}')THEN
              WRITE(30,100)'\\}'
            ELSE
              WRITE(30,100) OBJECT(L:L)
            END IF
          END DO
        END IF
        WRITE(30,101) '}'//' & & & & & \\\\ \\hline'
        GOTO 20
C------------------------------------------------------------------------------
70      CONTINUE
C Final del fichero LaTeX
        WRITE(30,101)'\\end{tabular}'
        WRITE(30,101)
        WRITE(30,101)'} %end of font size'
        WRITE(30,101)'\\end{document}'
        CLOSE(20)
        CLOSE(30)
C------------------------------------------------------------------------------
        STOP
100     FORMAT(A,$)
101     FORMAT(A)
        END
