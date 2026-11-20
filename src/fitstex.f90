!------------------------------------------------------------------------------
! Version 20-November-2026                                      file: fitstex.f
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
! Program: fitstex
! Classification: miscellany
! Description: Reads the output from fitshead and generates an output LaTeX 
! file with a table.
!
!Comment
!
! Este programa lee un fichero en formato ASCII generado con el programa
! fitshead, y produce un fichero en LaTeX con una tabla.
!
        PROGRAM FITSTEX
        IMPLICIT NONE
!
        INCLUDE 'futils.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'redlib.inc'
        INTEGER TRUELEN
        INTEGER READILIM
!
        INTEGER I,L,LL
        INTEGER NSKIP,NAXIS
        CHARACTER*1 CFORMAT
!cc        CHARACTER*10 FITSFILE
!cc        CHARACTER*20 OBJECT
        CHARACTER*255 INFILE,OUTFILE
        CHARACTER*255 CLINEA,CCOMMENT
!------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
!
        THISPROGRAM='fitstex'
        CALL WELCOME('20-November-2026')
!
        WRITE(*,100,ADVANCE='NO')'Input file (output from fitshead)'
        INFILE=INFILEX(20,'@',0,0,0.,0.,3,.FALSE.)
        WRITE(*,100,ADVANCE='NO')&
     &     'Format of previous file is [l]ong or [s]hort (l/s)'
        CFORMAT(1:1)=READC('@','ls')
        WRITE(*,100,ADVANCE='NO')'NAXIS '
        NAXIS=READILIM('2',1,3)
        WRITE(*,100,ADVANCE='NO')'Number of lines to be skipped '
        NSKIP=READILIM('1',0,9999)
        WRITE(*,100,ADVANCE='NO')'Ouput file (LaTeX)'
        OUTFILE=OUTFILEX(30,'@',0,0,0.,0.,3,.FALSE.)
        WRITE(*,100,ADVANCE='NO')'Comment'
        CCOMMENT=READC('@','@')
!------------------------------------------------------------------------------
! Escribimos cabecera del fichero LaTeX
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
!------------------------------------------------------------------------------
! Saltamos lineas a ignorar
        IF(NSKIP.GE.1)THEN
          DO I=1,NSKIP
            READ(20,101) CLINEA
          END DO
        END IF
!------------------------------------------------------------------------------
        I=0
! Leemos nombre del fichero FITS y nombre del objeto
20      IF(CFORMAT.EQ.'l')THEN
          IF(NAXIS.EQ.1)THEN
            READ(20,'(A10,10X,1(6X),A20)',END=70)FITSFILE,OBJECT
          ELSEIF(NAXIS.EQ.2)THEN
            READ(20,'(A10,10X,2(6X),A20)',END=70)FITSFILE,OBJECT
          ELSEIF(NAXIS.EQ.3)THEN
            READ(20,'(A10,10X,3(6X),A20)',END=70)FITSFILE,OBJECT
          ELSE
            WRITE(*,100,ADVANCE='NO')'FATAL ERROR: invalid NAXIS.'
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
            WRITE(*,100,ADVANCE='NO')'FATAL ERROR: invalid NAXIS.'
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
          WRITE(30,101)&
     &     '\\makebox[25mm][c]{\\large FITS file} '//&
     &     '& \\makebox[30mm][c]{\\large Object} '//&
     &     '& \\makebox[20mm][c]{}'//&
     &     '& \\makebox[20mm][c]{}'//&
     &     '& \\makebox[20mm][c]{}'//&
     &     '& \\makebox[20mm][c]{}'//&
     &     '& \\makebox[20mm][c]{} \\\\ \\hline'
        END IF
!------------------------------------------------------------------------------
        WRITE(30,100,ADVANCE='NO')&
     &     '{\\normalsize '//FITSFILE//'}'//' & '//&
     &   '{\\normalsize '
        LL=TRUELEN(OBJECT)
        IF(LL.GT.0)THEN
          DO L=1,LL
            IF(OBJECT(L:L).EQ.'$')THEN
              WRITE(30,100,ADVANCE='NO')'\\$'
            ELSEIF(OBJECT(L:L).EQ.'&')THEN
              WRITE(30,100,ADVANCE='NO')'\\&'
            ELSEIF(OBJECT(L:L).EQ.'%')THEN
              WRITE(30,100,ADVANCE='NO')'\\%'
            ELSEIF(OBJECT(L:L).EQ.'#')THEN
              WRITE(30,100,ADVANCE='NO')'\\#'
            ELSEIF(OBJECT(L:L).EQ.'_')THEN
              WRITE(30,100,ADVANCE='NO')'\\_'
            ELSEIF(OBJECT(L:L).EQ.'{')THEN
              WRITE(30,100,ADVANCE='NO')'\\{'
            ELSEIF(OBJECT(L:L).EQ.'}')THEN
              WRITE(30,100,ADVANCE='NO')'\\}'
            ELSE
              WRITE(30,100,ADVANCE='NO') OBJECT(L:L)
            END IF
          END DO
        END IF
        WRITE(30,101) '}'//' & & & & & \\\\ \\hline'
        GOTO 20
!------------------------------------------------------------------------------
70      CONTINUE
! Final del fichero LaTeX
        WRITE(30,101)'\\end{tabular}'
        WRITE(30,101)
        WRITE(30,101)'} %end of font size'
        WRITE(30,101)'\\end{document}'
        CLOSE(20)
        CLOSE(30)
!------------------------------------------------------------------------------
        STOP
100     FORMAT(A)
101     FORMAT(A)
        END
