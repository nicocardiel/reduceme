C------------------------------------------------------------------------------
C Version 13-November-1998                                      file:creahelp.f
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
C Program: creahelp
C Classification: miscellany
C Description: Creates the help file 'helpred.txt', and the auxiliary files
C 'allfiles.tex', and libraries.tex.
C
Comment
C
        PROGRAM CREAHELP
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER TRUELEN
        INTEGER READILIM
C
        INTEGER NCLASES
        PARAMETER(NCLASES=13)
        INTEGER NMAXFILES
        PARAMETER (NMAXFILES=200)
C
        INTEGER I,J,K,LL,LLL,LL1!,LL2
        INTEGER L1,L2
        INTEGER NFILES
        INTEGER NLINES
        INTEGER IOPC
        CHARACTER*1 CCONT,CSEP,CLAST
        CHARACTER*50 FILENAME(NMAXFILES),FILENAME0
        CHARACTER*50 CLASES(NCLASES),CLASE0
        CHARACTER*50 CP,CC
        CHARACTER*75 OUTFILE,LIBNAME
        CHARACTER*255 LINEA
        CHARACTER*800 DESCRIPCION
        LOGICAL LOGFILE,LP,LC,LD,LOK
        LOGICAL IFFILE(NMAXFILES)
C------------------------------------------------------------------------------
        DATA (CLASES(I),I=1,NCLASES)/
     +   'input/output',
     +   'examination & statistics',
     +   'arithmetic & manipulations',
     +   'distortion',
     +   'graphic display',
     +   'error handling',
     +   'cosmic rays',
     +   'wavelengths',
     +   'flux calibration',
     +   'extinction correction',
     +   'sky subtraction',
     +   'measurement',
     +   'miscellany'/
C
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
        INFILEX=INFILEX
C------------------------------------------------------------------------------
        WRITE(*,101)'(1) helpred.txt, helpred.html and helpreds.html'
        WRITE(*,101)'(2) library --routine header-- (LaTeX)'
        WRITE(*,101)'(3) library --routine header-- (HTML)'
        WRITE(*,100)'Option (1/2/3)'
        IOPC=READILIM('@',1,3)
C------------------------------------------------------------------------------
        INQUIRE(FILE='Allsrc.lis',EXIST=LOGFILE)
        IF(LOGFILE)THEN
          OPEN(20,FILE='Allsrc.lis',STATUS='OLD',FORM='FORMATTED')
          J=0
5         READ(20,101,END=7) FILENAME(J+1)
          IFFILE(J+1)=.FALSE.
          J=J+1
          IF(J+1.GT.NMAXFILES)THEN
            WRITE(*,101)'FATAL ERROR: redim NMAXFILES.'
            CLOSE(20)
            STOP
          END IF
          GOTO 5
7         CLOSE(20)
          NFILES=J
          WRITE(*,100)'No. of source files read: '
          WRITE(*,*) NFILES
        ELSE
          WRITE(*,101)'ERROR: file Allsrc.lis has not been found.'
          STOP
        END IF
C
        IF(IOPC.EQ.2) GOTO 50
        IF(IOPC.EQ.3) GOTO 60
C
        INQUIRE(FILE='helpred.txt',EXIST=LOGFILE)
        IF(LOGFILE)THEN
          WRITE(*,101)'ERROR: file helpred.txt already exist.'
          WRITE(*,101)'>> file helpred.txt will not be recreated.'
          STOP
        END IF
        OPEN(30,FILE='helpred.txt',STATUS='NEW',FORM='FORMATTED')
C
        INQUIRE(FILE='helpred.html',EXIST=LOGFILE)
        IF(LOGFILE)THEN
          WRITE(*,101)'ERROR: file helpred.html already exist.'
          WRITE(*,101)'>> file helpred.html will not be recreated.'
          STOP
        END IF
        OPEN(33,FILE='helpred.html',STATUS='NEW',FORM='FORMATTED')
C
        WRITE(33,101)'<HTML>'
        WRITE(33,101)'<HEAD>'
        WRITE(33,101)'   <TITLE>REDUCEME: Classified list</TITLE>'
        WRITE(33,100)'   <meta HTTP-EQUIV="content-type" '
        WRITE(33,101)'CONTENT="text/html; charset=ISO-8859-1">'
        WRITE(33,101)'</HEAD>'
        WRITE(33,101)'<BODY BGCOLOR="#CCCCCC" TEXT="#333333"'
        WRITE(33,101)' LINK="#000066" ALINK="#000066" VLINK="#000044">'
        WRITE(33,*)
        WRITE(33,101)'<A NAME="indice_temas"></A>'
        WRITE(33,101)'<TABLE WIDTH="100%" BORDER=1>'
        WRITE(33,101)'<TR>'
        WRITE(33,101)'  <TD ALIGN="RIGHT">'
        WRITE(33,101)'    <A HREF="reduceme.html">'
        WRITE(33,100)'    <IMG src="gif/home.gif" BORDER="0" '
        WRITE(33,101)'ALT="logo home page"><IMG '
        WRITE(33,100)'    src="gif/logo.gif" '
        WRITE(33,101)'BORDER="0" ALT="logo home page"></A>'
        WRITE(33,101)'  </TD>'
        WRITE(33,101)'</TR>'
        WRITE(33,101)'</TABLE>'
        WRITE(33,*)
        WRITE(33,101)'<HR>'
        WRITE(33,101)'<CENTER>'
        WRITE(33,101)'<A HREF="#INPUT">Input/Output</A> |'
        WRITE(33,100)'<A HREF="#EXAMI">Examination &amp '
        WRITE(33,101)' Statistics</A> |'
        WRITE(33,100)'<A HREF="#ARITH">Arithmetic &amp '
        WRITE(33,101)'Manipulations</A> |'
        WRITE(33,101)'<A HREF="#DISTO">Distortion</A> |'
        WRITE(33,101)'<A HREF="#GRAPH">Graphic Display</A> |'
        WRITE(33,101)'<A HREF="#ERROR">Error Handling</A> |'
        WRITE(33,101)'<A HREF="#COSMI">Cosmic Rays</A> |'
        WRITE(33,101)'<A HREF="#WAVEL">Wavelengths</A> |'
        WRITE(33,101)'<A HREF="#FLUXC">Flux Calibration</A> |'
        WRITE(33,101)'<A HREF="#EXTIN">Extinction Correction</A> |'
        WRITE(33,101)'<A HREF="#SKYSU">Sky Subtraction</A> |'
        WRITE(33,101)'<A HREF="#MEASU">Measurement</A> |'
        WRITE(33,101)'<A HREF="#MISCE">Miscellany</A>'
        WRITE(33,101)'</CENTER>'
        WRITE(33,*)
C------------------------------------------------------------------------------
        CALL PRINTOUT(30,'CLASSIFIED LIST OF COMMANDS',1)
        WRITE(*,100)'Working...'
        DO I=1,NCLASES
ccc       WRITE(*,100)'Next class is: '
ccc       WRITE(*,101)CLASES(I)
          CLASE0=CLASES(I)
          CALL CHUPPER(CLASE0)
          CALL PRINTOUT(30,CLASE0,2)
          CALL PRINTOUT(33,CLASE0,3)
          WRITE(33,101)'<UL>'
          DO J=1,NFILES
            FILENAME0=FILENAME(J)
            INQUIRE(FILE=FILENAME0,EXIST=LOGFILE)
            IF(LOGFILE)THEN
              OPEN(40,FILE=FILENAME0,STATUS='OLD',FORM='FORMATTED')
              LP=.FALSE.
              LC=.FALSE.
              LD=.FALSE.
              LOK=.FALSE.
              DO K=1,LEN(DESCRIPCION)
                DESCRIPCION(K:K)=' '
              END DO
12            READ(40,101,END=18)LINEA
              IF(.NOT.LOK)THEN
                LOK=((LINEA(1:7).EQ.'Comment').OR.
     +           (LINEA(1:9).EQ.'/*Comment'))       !buscamos primer comentario
                IF(.NOT.LOK) GOTO 12
              END IF
              IF(LP)THEN
                IF(LC)THEN
                  IF(LD)THEN
                    IF((LINEA(1:7).EQ.'Comment').OR.
     +               (LINEA(1:9).EQ.'/*Comment'))THEN        !ultimo comentario
                      WRITE(30,*)
                      WRITE(33,*)
                      CALL PRINTDESCRIPTION(0,DESCRIPCION,79)
                      IFFILE(J)=.TRUE.
                      GOTO 18
                    END IF
                    L1=TRUELEN(LINEA)
                    IF(LINEA(1:2).EQ.'/*') L1=TRUELEN(LINEA(1:L1-2))
                    IF(L1.LT.3)THEN
                    ELSE
                      L2=TRUELEN(DESCRIPCION)
                      DESCRIPCION(L2+1:L2+1)=' '
                      DESCRIPCION(L2+2:L2+2+L1-1-2)=LINEA(1+2:L1)
                    END IF
                  ELSE
                    IF(CC.EQ.CLASES(I))THEN
                      IF((LINEA(1:15).EQ.'C Description: ').OR.
     +                 (LINEA(1:15).EQ.'/*Description: '))THEN
                        LD=.TRUE.
                        L1=TRUELEN(LINEA)
                        IF(LINEA(1:2).EQ.'/*') L1=TRUELEN(LINEA(1:L1-2))
                        L2=TRUELEN(DESCRIPCION)
                        DESCRIPCION(L2+1:L2+1)=' '
                        DESCRIPCION(L2+1+1:L2+1+L1-16+1)=LINEA(16:L1)
                      END IF
                    ELSE
                      GOTO 18
                    END IF
                  END IF
                ELSE
                  IF((LINEA(1:18).EQ.'C Classification: ').OR.
     +             (LINEA(1:18).EQ.'/*Classification: '))THEN
                    LC=.TRUE.
                    L1=TRUELEN(LINEA)
                    IF(LINEA(1:2).EQ.'/*') L1=TRUELEN(LINEA(1:L1-2))
                    CC=LINEA(19:L1)
                  END IF
                END IF
              ELSE
                IF((LINEA(1:11).EQ.'C Program: ').OR.
     +           (LINEA(1:11).EQ.'/*Program: '))THEN
                  LP=.TRUE.
                  L1=TRUELEN(LINEA)
                  IF(LINEA(1:2).EQ.'/*') L1=TRUELEN(LINEA(1:L1-2))
                  CP=LINEA(12:L1)
                  L1=TRUELEN(CP)
                  DESCRIPCION(1:L1)=CP(1:L1)
                  DESCRIPCION(L1+1:L1+1)=':'
                END IF
              END IF
              GOTO 12
18            CLOSE(40)
C..............................................................................
            ELSE
              WRITE(*,*)
              WRITE(*,100)'ERROR: file '
              WRITE(*,100)FILENAME0(1:TRUELEN(FILENAME0))
              WRITE(*,101)' does not exist.'
              WRITE(*,100)'Do you want to continue (y/n) '
              CCONT(1:1)=READC('y','yn')
              IF(CCONT.EQ.'n')THEN
                CLOSE(20)
                CLOSE(30)
                CLOSE(33)
                STOP
              END IF
            END IF
          END DO
          WRITE(33,101)'</UL>'
        END DO
C------------------------------------------------------------------------------
        WRITE(30,*)
        WRITE(33,*)
        WRITE(33,101)'<P>'
        WRITE(33,101)'<TABLE WIDTH="100%" BORDER=1>'
        WRITE(33,101)'<TR>'
        WRITE(33,101)'  <TD ALIGN="RIGHT">'
        WRITE(33,101)'    <A HREF="reduceme.html">'
        WRITE(33,100)'    <IMG src="gif/home.gif" BORDER="0" '
        WRITE(33,101)'ALT="logo home page"><IMG '
        WRITE(33,100)'    src="gif/logo.gif" '
        WRITE(33,101)'BORDER="0" ALT="logo home page"></A>'
        WRITE(33,101)'  </TD>'
        WRITE(33,101)'</TR>'
        WRITE(33,101)'</TABLE>'
        WRITE(33,*)
        WRITE(33,101)'</BODY>'
        WRITE(33,101)'</HTML>'
        WRITE(30,100)'                                        '
        WRITE(30,101)'                      (end of document)'
        WRITE(30,100)'----------------------------------------'
        WRITE(30,101)'---------------------------------------'
        CLOSE(30)
        CLOSE(33)
        WRITE(*,101)'   ...OK!'
C------------------------------------------------------------------------------
        DO J=1,NFILES
          IF(.NOT.IFFILE(J))THEN
            WRITE(*,100)'WARNING: something was wrong with file: '
            WRITE(*,101)FILENAME(J)(1:TRUELEN(FILENAME(J)))
          END IF
        END DO
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
        INQUIRE(FILE='allfiles.tex',EXIST=LOGFILE)
        IF(LOGFILE)THEN
          WRITE(*,101)'ERROR: file allfiles.tex already exist.'
          WRITE(*,101)'>> file allfiles.tex will not be recreated.'
          STOP
        END IF
        OPEN(32,FILE='allfiles.tex',STATUS='NEW',FORM='FORMATTED')
C
        INQUIRE(FILE='helpreds.html',EXIST=LOGFILE)
        IF(LOGFILE)THEN
          WRITE(*,101)'ERROR: file helpreds.html already exist.'
          WRITE(*,101)'>> file helpreds.html will not be recreated.'
          STOP
        END IF
        OPEN(36,FILE='helpreds.html',STATUS='NEW',FORM='FORMATTED')
C
        WRITE(36,101)'<HTML>'
        WRITE(36,101)'<HEAD>'
        WRITE(36,101)'   <TITLE>REDUCEME: Sorted list</TITLE>'
        WRITE(36,101)'</HEAD>'
        WRITE(36,101)'<BODY BGCOLOR="#CCCCCC" TEXT="#333333"'
        WRITE(36,101)' LINK="#000066" ALINK="#000066" VLINK="#000044">'
        WRITE(36,*)
        WRITE(36,101)'<A NAME="alfabeto"></A>'
        WRITE(36,101)'<TABLE WIDTH="100%" BORDER=1>'
        WRITE(36,101)'<TR>'
        WRITE(36,101)'  <TD ALIGN="RIGHT">'
        WRITE(36,101)'    <A HREF="reduceme.html">'
        WRITE(36,100)'    <IMG src="gif/home.gif" BORDER="0" '
        WRITE(36,101)'ALT="logo home page"><IMG '
        WRITE(36,100)'    src="gif/logo.gif" '
        WRITE(36,101)'BORDER="0" ALT="logo home page"></A>'
        WRITE(36,101)'  </TD>'
        WRITE(36,101)'</TR>'
        WRITE(36,101)'</TABLE>'
        WRITE(36,*)
        WRITE(36,101)'<HR>'
        WRITE(36,101)'<H3 ALIGN=CENTER>'
        WRITE(36,101)'<A HREF="#a">a</A> | '
        WRITE(36,101)'<A HREF="#b">b</A> | '
        WRITE(36,101)'<A HREF="#c">c</A> | '
        WRITE(36,101)'<A HREF="#d">d</A> | '
        WRITE(36,101)'<A HREF="#e">e</A> | '
        WRITE(36,101)'<A HREF="#f">f</A> | '
        WRITE(36,101)'<A HREF="#g">g</A> | '
        WRITE(36,101)'<A HREF="#h">h</A> | '
        WRITE(36,101)'<A HREF="#i">i</A> | '
        WRITE(36,101)'<A HREF="#j">j</A> | '
        WRITE(36,101)'<A HREF="#k">k</A> | '
        WRITE(36,101)'<A HREF="#l">l</A> | '
        WRITE(36,101)'<A HREF="#m">m</A> | '
        WRITE(36,101)'<A HREF="#n">n</A> | '
        WRITE(36,101)'<A HREF="#o">o</A> | '
        WRITE(36,101)'<A HREF="#p">p</A> | '
        WRITE(36,101)'<A HREF="#q">q</A> | '
        WRITE(36,101)'<A HREF="#r">r</A> | '
        WRITE(36,101)'<A HREF="#s">s</A> | '
        WRITE(36,101)'<A HREF="#t">t</A> | '
        WRITE(36,101)'<A HREF="#u">u</A> | '
        WRITE(36,101)'<A HREF="#v">v</A> | '
        WRITE(36,101)'<A HREF="#w">w</A> | '
        WRITE(36,101)'<A HREF="#x">x</A> | '
        WRITE(36,101)'<A HREF="#y">y</A> | '
        WRITE(36,101)'<A HREF="#z">z</A>'
        WRITE(36,101)'</H3><HR>'
        WRITE(36,*)
        WRITE(36,101)'<UL>'
C------------------------------------------------------------------------------
        CLAST=' '
        WRITE(*,100)'Working...'
        DO J=1,NFILES
          FILENAME0=FILENAME(J)
          OPEN(40,FILE=FILENAME0,STATUS='OLD',FORM='FORMATTED')
          LP=.FALSE.
          LC=.FALSE.
          LD=.FALSE.
          LOK=.FALSE.
          DO K=1,LEN(DESCRIPCION)
            DESCRIPCION(K:K)=' '
          END DO
22        READ(40,101,END=28)LINEA
          IF(.NOT.LOK)THEN
            LOK=((LINEA(1:7).EQ.'Comment').OR.
     +       (LINEA(1:9).EQ.'/*Comment'))           !buscamos primer comentario
            IF(.NOT.LOK) GOTO 22
          END IF
          IF(LP)THEN
            IF(LC)THEN
              IF(LD)THEN
                IF((LINEA(1:7).EQ.'Comment').OR.
     +           (LINEA(1:9).EQ.'/*Comment'))THEN            !ultimo comentario
                  CALL PRINTDESCRIPTION(1,DESCRIPCION,79)
                  GOTO 28
                END IF
                L1=TRUELEN(LINEA)
                IF(LINEA(1:2).EQ.'/*') L1=TRUELEN(LINEA(1:L1-2))
                IF(L1.LT.3)THEN
                ELSE
                  L2=TRUELEN(DESCRIPCION)
                  DESCRIPCION(L2+1:L2+1)=' '
                  DESCRIPCION(L2+2:L2+2+L1-1-2)=LINEA(1+2:L1)
                END IF
              ELSE
                IF((LINEA(1:15).EQ.'C Description: ').OR.
     +           (LINEA(1:15).EQ.'/*Description: '))THEN
                  LD=.TRUE.
                  L1=TRUELEN(LINEA)
                  IF(LINEA(1:2).EQ.'/*') L1=TRUELEN(LINEA(1:L1-2))
                  L2=TRUELEN(DESCRIPCION)
                  DESCRIPCION(1:L1-16+1)=LINEA(16:L1)
                END IF
              END IF
            ELSE
              IF((LINEA(1:18).EQ.'C Classification: ').OR.
     +         (LINEA(1:18).EQ.'/*Classification: '))THEN
                LC=.TRUE.
                L1=TRUELEN(LINEA)
                IF(LINEA(1:2).EQ.'/*') L1=TRUELEN(LINEA(1:L1-2))
                CC=LINEA(19:L1)
              END IF
            END IF
          ELSE
            IF((LINEA(1:11).EQ.'C Program: ').OR.
     +       (LINEA(1:11).EQ.'/*Program: '))THEN
              LP=.TRUE.
              L1=TRUELEN(LINEA)
              IF(LINEA(1:2).EQ.'/*') L1=TRUELEN(LINEA(1:L1-2))
              CP=LINEA(12:L1)
              L1=TRUELEN(CP)
              WRITE(32,*)
              WRITE(32,100)'\\subsection{'
              WRITE(32,100)CP(1:L1)
              WRITE(32,101)'}'
              WRITE(36,*)
              IF(CP(1:1).EQ.CLAST)THEN
                WRITE(36,100)'<LI>'
              ELSE
                CLAST=CP(1:1)
                WRITE(36,101)'<P>'
                WRITE(36,*)
                WRITE(36,100)'<A NAME="'
                WRITE(36,100)CLAST
                WRITE(36,101)'"></A><HR>'
                WRITE(36,101)'<BR>'
                WRITE(36,101)'<TABLE>'
                WRITE(36,101)'<TR>'
                WRITE(36,100)'  <TD ALIGN="CENTER">'
                WRITE(36,100)'<H2><FONT COLOR="#008000">'
                WRITE(36,100)CLAST
                WRITE(36,101)'</FONT></H2></TD>'
                WRITE(36,100)'  <TD ALIGN="RIGHT"><A HREF='
                WRITE(36,101)'"#alfabeto"><IMG src="gif/up.gif" '
                WRITE(36,101)'   ALT="up.gif"></A></TD>'
                WRITE(36,101)'</TR>'
                WRITE(36,101)'</TABLE>'
                WRITE(36,101)'<BR>'
                WRITE(36,*)
                WRITE(36,100)'<LI>'
              END IF
              LLL=INDEX(CP(1:L1),'.')
              WRITE(36,100)'<A HREF="programs/'
              IF(LLL.NE.0)THEN
                WRITE(36,100)CP(1:LLL-1)
              ELSE
                WRITE(36,100)CP(1:L1)
              END IF
              WRITE(36,100)'.html">'
              WRITE(36,100)CP(1:L1)
              WRITE(36,101)'</A>: '
            END IF
          END IF
          GOTO 22
28        CLOSE(40)
        END DO
C------------------------------------------------------------------------------
        WRITE(36,101)'</UL>'
        WRITE(36,*)
        WRITE(36,101)'<P>'
        WRITE(36,101)'<TABLE WIDTH="100%" BORDER=1>'
        WRITE(36,101)'<TR>'
        WRITE(36,101)'  <TD ALIGN="RIGHT">'
        WRITE(36,101)'    <A HREF="reduceme.html">'
        WRITE(36,100)'    <IMG src="gif/home.gif" BORDER="0" '
        WRITE(36,101)'ALT="logo home page"><IMG '
        WRITE(36,100)'    src="gif/logo.gif" '
        WRITE(36,101)'BORDER="0" ALT="logo home page"></A>'
        WRITE(36,101)'  </TD>'
        WRITE(36,101)'</TR>'
        WRITE(36,101)'</TABLE>'
        WRITE(36,*)
        WRITE(36,101)'</BODY>'
        WRITE(36,101)'</HTML>'
        CLOSE(32)
        CLOSE(36)
        WRITE(*,101)'   ...OK!'
        STOP
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
50      CONTINUE
        WRITE(*,100)'Ouput LaTeX file'
        OUTFILE=OUTFILEX(30,'@',0,0,0.,0.,3,.FALSE.)
C------------------------------------------------------------------------------
        WRITE(*,100)'Working...'
        DO J=1,NFILES
          FILENAME0=FILENAME(J)
          OPEN(40,FILE=FILENAME0,STATUS='OLD',FORM='FORMATTED')
          LOK=.FALSE.
          NLINES=0
          DO K=1,LEN(DESCRIPCION)
            DESCRIPCION(K:K)=' '
          END DO
52        READ(40,101,END=58)LINEA
          CSEP='+'
          L1=TRUELEN(LINEA)
          IF(L1.GT.0)THEN
            DO WHILE(INDEX(LINEA(1:L1),CSEP).NE.0)       !buscamos un separador
              IF(CSEP.EQ.'+')THEN
                CSEP='('
              ELSEIF(CSEP.EQ.'(')THEN
                CSEP=')'
              ELSEIF(CSEP.EQ.')')THEN
                CSEP='/'
              ELSEIF(CSEP.EQ.'/')THEN
                CSEP='%'
              ELSEIF(CSEP.EQ.'%')THEN
                CSEP='$'
              ELSEIF(CSEP.EQ.'$')THEN
                CSEP='#'
              ELSE
                WRITE(*,101)'FATAL ERROR: invalid separator.'
                CLOSE(30)
                CLOSE(40)
                STOP
              END IF
            END DO
          END IF
          IF(.NOT.LOK)THEN
            LOK=((LINEA(1:7).EQ.'Comment').OR.
     +       (LINEA(1:9).EQ.'/*Comment'))           !buscamos primer comentario
            IF(.NOT.LOK) GOTO 52
            READ(40,*)                                      !saltamos una linea
            READ(40,101,END=58)LINEA
            WRITE(30,*)
            WRITE(30,101)'\\begin{minipage}{\\textwidth}'
            WRITE(30,100)'{\\large $\\rhd$ \\verb'
            WRITE(30,100)CSEP
            LL=TRUELEN(LINEA)
            IF(LINEA(LL-1:LL).EQ.'*/') LL=TRUELEN(LINEA(1:LL-2))
            WRITE(30,100)LINEA(3:LL)                    !nombre de la subrutina
            WRITE(30,100)CSEP
            WRITE(30,101)'}'
            WRITE(30,*)
            WRITE(30,101)'\\vspace{2mm}'
            WRITE(30,*)
            WRITE(30,101)'{\\normalsize'
            GOTO 52
          END IF
          IF((LINEA(1:7).EQ.'Comment').OR.
     +       (LINEA(1:9).EQ.'/*Comment'))THEN                !ultimo comentario
            WRITE(30,101)'}%end of normalsize'
            WRITE(30,101)'\\end{minipage}'
            WRITE(30,*)
            GOTO 58
          ELSE
            NLINES=NLINES+1
            IF(MOD(NLINES,2).EQ.0)THEN
              WRITE(30,101)'\\makebox[10mm][l]{'//
     +        '\\raisebox{0pt}[0pt][0pt]{\\vrule height13pt depth2pt}'//
     +        '\\hfill '//
     +        '$\\bigcirc$'//
     +        '\\hfill '//
     +        '\\raisebox{0pt}[0pt][0pt]{\\vrule height7pt depth0pt}'//
     +        '}'
            ELSE
              WRITE(30,101)'\\makebox[10mm][l]{'//
     +        '\\raisebox{0pt}[0pt][0pt]{\\vrule height13pt depth2pt}'//
     +        '\\hfill '//
     +        '\\raisebox{0pt}[0pt][0pt]{\\vrule height7pt depth0pt}'//
     +        '}'
            END IF
            IF(L1.GT.0)THEN
              WRITE(30,100)'\\verb'
              WRITE(30,100)CSEP
              WRITE(30,100)LINEA(1:L1)
              WRITE(30,101)CSEP
            END IF
            WRITE(30,101)'\\hfill'
            IF(MOD(NLINES,2).EQ.0)THEN
              WRITE(30,101)'\\makebox[10mm][r]{'//
     +        '\\raisebox{0pt}[0pt][0pt]{\\vrule height7pt depth0pt}'//
     +        '\\hfill '//
     +        '$\\bigcirc$'//
     +        '\\hfill '//
     +        '\\raisebox{0pt}[0pt][0pt]{\\vrule height13pt depth2pt}'//
     +        '}'
            ELSE
              WRITE(30,101)'\\makebox[10mm][r]{'//
     +        '\\raisebox{0pt}[0pt][0pt]{\\vrule height7pt depth0pt}'//
     +        '\\hfill '//
     +        '\\raisebox{0pt}[0pt][0pt]{\\vrule height13pt depth2pt}'//
     +        '}'
            END IF
            WRITE(30,101)'\\\\'
          END IF
          GOTO 52
58        CLOSE(40)
        END DO
C------------------------------------------------------------------------------
        CLOSE(30)
        WRITE(*,101)'   ...OK!'
        STOP
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
60      CONTINUE
        WRITE(*,100)'Library description'
        LIBNAME(1:75)=READC('@','@')
        WRITE(*,100)'Ouput HTML file'
        OUTFILE=OUTFILEX(30,'@',0,0,0.,0.,3,.FALSE.)
        WRITE(30,101)'<HTML>'
        WRITE(30,101)'<HEAD>'
        IF(LIBNAME(1:TRUELEN(LIBNAME)).EQ.'libbutton.a')THEN
          WRITE(30,100)'   <TITLE>BUTTON: '
        ELSE
          WRITE(30,100)'   <TITLE>REDUCEME: '
        END IF
        WRITE(30,101)'sorted list of subroutines</TITLE>'
        WRITE(30,101)'   <STYLE type="text/css">'
        WRITE(30,101)'     PRE.subname {color: rgb(105,0,0);'
        WRITE(30,101)'                   font-size: medium}'
        WRITE(30,101)'     PRE.subtext {color: rgb(0,0,0)}'
        WRITE(30,101)'   </STYLE>'
        WRITE(30,101)'</HEAD>'
        WRITE(30,*)
        WRITE(30,101)'<BODY BGCOLOR="#CCCCCC" TEXT="#333333"'
        WRITE(30,101)' LINK="#000066" ALINK="#000066" VLINK="#000044">'
        WRITE(30,*)
        WRITE(30,101)'<A NAME="indice_rutinas"></A>'
        WRITE(30,101)'<TABLE WIDTH="100%" BORDER=1>'
        WRITE(30,101)'<TR>'
        WRITE(30,101)'  <TD ALIGN="RIGHT">'
        WRITE(30,101)'    <A HREF="reduceme.html">'
        WRITE(30,100)'    <IMG src="gif/home.gif" BORDER="0" '
        WRITE(30,101)'ALT="logo home page"><IMG '
        WRITE(30,100)'    src="gif/logo.gif" '
        WRITE(30,101)'BORDER="0" ALT="logo home page"></A>'
        WRITE(30,101)'  </TD>'
        WRITE(30,101)'</TR>'
        WRITE(30,101)'</TABLE>'
        WRITE(30,*)
        WRITE(30,101)'<P>'
        WRITE(30,101)'<FONT COLOR="#666666">'
        WRITE(30,100)'<CENTER><H1>'
        WRITE(30,100)LIBNAME(1:TRUELEN(LIBNAME))
        WRITE(30,101)'</H1></CENTER></FONT>'
        WRITE(30,101)'<P>'
        WRITE(30,101)'<HR>'
        WRITE(30,101)'<CENTER>'
        DO J=1,NFILES
          FILENAME0=FILENAME(J)
          L1=TRUELEN(FILENAME0)
          IF(FILENAME0(L1-1:L1).EQ.'.f') L1=L1-2
          IF(FILENAME0(L1-1:L1).EQ.'.c') L1=L1-2
          IF(J.EQ.1)THEN
            WRITE(30,100)'<A HREF="#'
            WRITE(30,100)FILENAME0(1:L1)
            WRITE(30,100)'">'
            WRITE(30,100)FILENAME0(1:L1)
            WRITE(30,101)'</A>'
          ELSE
            WRITE(30,100)' | '
            WRITE(30,100)'<A HREF="#'
            WRITE(30,100)FILENAME0(1:L1)
            WRITE(30,100)'">'
            WRITE(30,100)FILENAME0(1:L1)
            WRITE(30,101)'</A>'
          END IF
        END DO
        WRITE(30,101)'</CENTER>'
        WRITE(30,101)'<HR>'
        WRITE(30,101)'<P>'
        WRITE(30,*)
C------------------------------------------------------------------------------
        WRITE(*,100)'Working...'
        DO J=1,NFILES
          FILENAME0=FILENAME(J)
          OPEN(40,FILE=FILENAME0,STATUS='OLD',FORM='FORMATTED')
          LOK=.FALSE.
          NLINES=0
          DO K=1,LEN(DESCRIPCION)
            DESCRIPCION(K:K)=' '
          END DO
62        READ(40,101,END=68)LINEA
          L1=TRUELEN(LINEA)
          IF(.NOT.LOK)THEN
            LOK=((LINEA(1:7).EQ.'Comment').OR.
     +       (LINEA(1:9).EQ.'/*Comment'))           !buscamos primer comentario
            IF(.NOT.LOK) GOTO 62
            READ(40,*)                                      !saltamos una linea
            READ(40,101,END=68)LINEA
            WRITE(30,*)
            WRITE(30,100)'<A NAME="'
            LL1=TRUELEN(FILENAME0)
            IF(FILENAME0(LL1-1:LL1).EQ.'.f') LL1=LL1-2
            IF(FILENAME0(LL1-1:LL1).EQ.'.c') LL1=LL1-2
ccc         LL1=INDEX(LINEA(3:),' ')
ccc         LL2=INDEX(LINEA(3:),'(')
ccc         IF(LL2.EQ.0)LL2=INDEX(LINEA(LL1+3:),' ')
ccc         WRITE(30,100)LINEA(LL1+3:LL2+1)
            WRITE(30,100)FILENAME0(1:LL1)
            WRITE(30,100)'">'
            IF(J.NE.1) WRITE(30,100)'<HR>'
            WRITE(30,101)'</A><P>'
            WRITE(30,100)'<A HREF="#indice_rutinas">'
            WRITE(30,100)'<IMG src="gif/up.gif" ALT="up.gif" '
            WRITE(30,101)'ALIGN="RIGHT"></A>'
            WRITE(30,101)'<PRE class="subname">'
            LL=TRUELEN(LINEA)
            IF(LINEA(LL-1:LL).EQ.'*/') LL=TRUELEN(LINEA(1:LL-2))
            WRITE(30,101)LINEA(3:LL)                    !nombre de la subrutina
            READ(40,101,END=68)LINEA
            LL=TRUELEN(LINEA)
            IF(LINEA(LL-1:LL).EQ.'*/') LL=TRUELEN(LINEA(1:LL-2))
            IF(LL.GT.2)THEN    !hay mas de una linea con el nombre de la rutina
              DO WHILE(LL.GT.2)
                WRITE(30,101)LINEA(3:LL)   !continuacion nombre de la subrutina
                READ(40,101,END=68)LINEA
                LL=TRUELEN(LINEA)
                IF(LINEA(LL-1:LL).EQ.'*/') LL=TRUELEN(LINEA(1:LL-2))
              END DO
            END IF
            WRITE(30,101)'</PRE>'
            WRITE(30,*)
            WRITE(30,101)'<PRE class="subtext">'
            GOTO 62
          END IF
          IF((LINEA(1:7).EQ.'Comment').OR.
     +       (LINEA(1:9).EQ.'/*Comment'))THEN                !ultimo comentario
            WRITE(30,101)'</PRE>'
            WRITE(30,*)
            GOTO 68
          ELSE
ccc         NLINES=NLINES+1
ccc         IF(MOD(NLINES,2).EQ.0)THEN
ccc         ELSE
ccc         END IF
            IF(L1.GT.0)THEN
              IF(LINEA(L1-1:L1).EQ.'*/') L1=TRUELEN(LINEA(1:L1-2))
              WRITE(30,101)LINEA(3:L1)
            END IF
ccc         IF(MOD(NLINES,2).EQ.0)THEN
ccc         ELSE
ccc         END IF
          END IF
          GOTO 62
68        CLOSE(40)
          IF(J.NE.NFILES)THEN
            IF(NLINES.GT.0)THEN
              WRITE(30,101)'<HR>'
              WRITE(30,101)'<P>'
            END IF
          END IF
        END DO
C------------------------------------------------------------------------------
        WRITE(30,101)'<TABLE WIDTH="100%" BORDER=1>'
        WRITE(30,101)'<TR>'
        WRITE(30,101)'  <TD ALIGN="RIGHT">'
        WRITE(30,101)'    <A HREF="reduceme.html">'
        WRITE(30,100)'    <IMG src="gif/home.gif" BORDER="0" '
        WRITE(30,101)'ALT="logo home page"><IMG '
        WRITE(30,100)'    src="gif/logo.gif" '
        WRITE(30,101)'BORDER="0" ALT="logo home page"></A>'
        WRITE(30,101)'  </TD>'
        WRITE(30,101)'</TR>'
        WRITE(30,101)'</TABLE>'
        CLOSE(30)
        WRITE(*,101)'   ...OK!'
        STOP
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
        END
C
C******************************************************************************
C muestra una cadena centrada en la unidad UNIT
C MODE:
C       1 con linea de '=' arriba y abajo
C       2 con linea de '=' abajo y linea completa de '-' arriba
C       3 con anchor para fichero html
C
        SUBROUTINE PRINTOUT(UNIT,CADENA,MODE)
        IMPLICIT NONE
        INTEGER TRUELEN
C
        INTEGER UNIT
        CHARACTER*(*) CADENA
        INTEGER MODE
C
        INTEGER I,L,LL,K
        CHARACTER*80 BLANCO,DOBLE
        CHARACTER*80 CDUMMY
C------------------------------------------------------------------------------
        DO I=1,80
          BLANCO(I:I)=' '
          DOBLE(I:I)='='
        END DO
C
        L=TRUELEN(CADENA)
        LL=(80-L)/2
C------------------------------------------------------------------------------
        IF(MODE.EQ.1)THEN
          WRITE(UNIT,100)BLANCO(1:LL)
          WRITE(UNIT,101)DOBLE(1:L)
          WRITE(UNIT,100)BLANCO(1:LL)
          WRITE(UNIT,101)CADENA(1:L)
          WRITE(UNIT,100)BLANCO(1:LL)
          WRITE(UNIT,101)DOBLE(1:L)
C------------------------------------------------------------------------------
        ELSEIF(MODE.EQ.2)THEN
          WRITE(UNIT,*)
          WRITE(UNIT,150)
          WRITE(UNIT,100)BLANCO(1:LL)
          WRITE(UNIT,101)CADENA(1:L)
          WRITE(UNIT,100)BLANCO(1:LL)
          WRITE(UNIT,101)DOBLE(1:L)
C------------------------------------------------------------------------------
        ELSEIF(MODE.EQ.3)THEN
          CDUMMY=CADENA
          CALL RMBLANK(CDUMMY,CDUMMY,K)
          IF(K.GT.5)THEN
            WRITE(UNIT,101)'<!---------------------------------------'//
     +       '------------------------------------->'
            WRITE(UNIT,100)'<A NAME="'
            WRITE(UNIT,100)CDUMMY(1:5)
            WRITE(UNIT,100)'"></A>'
            WRITE(UNIT,*)
            WRITE(UNIT,101)'<HR>'
            WRITE(UNIT,*)
            WRITE(UNIT,101)'<TABLE WIDTH="100%">'
            WRITE(UNIT,101)'<TR>'
            WRITE(UNIT,100)'  <TD ALIGN="LEFT">'
            WRITE(UNIT,100)'<H2><FONT COLOR="#008000">'
            WRITE(UNIT,100)CADENA(1:L)
            WRITE(UNIT,101)'</FONT></H2></TD>'
            WRITE(UNIT,100)'  <TD ALIGN="RIGHT">'
            WRITE(UNIT,100)'<A HREF="#indice_temas">'
            WRITE(UNIT,100)'<IMG src="gif/up.gif" ALT="up.gif">'
            WRITE(UNIT,101)'</A></TD>'
            WRITE(UNIT,101)'</TR>'
            WRITE(UNIT,101)'</TABLE>'
            WRITE(UNIT,*)
          ELSE
            WRITE(*,101)'FATAL ERROR: in subroutine PRINTOUT:'
            WRITE(*,100)'K='
            WRITE(*,*) K
            WRITE(*,101)'(K out of range)'
            STOP
          END IF
C------------------------------------------------------------------------------
        ELSE
          WRITE(*,101)'FATAL ERROR: in subroutine PRINTOUT:'
          WRITE(*,100)'MODE='
          WRITE(*,*) MODE
          WRITE(*,101)'(MODE out of range)'
          STOP
        END IF
C
100     FORMAT(A,$)
101     FORMAT(A)
150     FORMAT(79('-'))
        END
C
C******************************************************************************
C Muestra la cadena en la unidad UNIT con lineas de LWIDTH caracteres
C Nota: si UNIT=0 escribe en las unidades 30 y 33
C Nota: si UNIT=1 escribe en las unidades 32 y 36
        SUBROUTINE PRINTDESCRIPTION(UNIT,CADENA,LWIDTH)
        IMPLICIT NONE
        INTEGER TRUELEN
C
        INTEGER UNIT
        CHARACTER*(*) CADENA
        INTEGER LWIDTH
C
        INTEGER I,L,LL,LLL
        LOGICAL LLOOP
C------------------------------------------------------------------------------
        LL=0 !avoid warning compilation
C
        L=TRUELEN(CADENA)
        IF(L.LE.LWIDTH)THEN
          IF(UNIT.EQ.0)THEN
            WRITE(30,101) CADENA(1:L)
            WRITE(33,101)
            LL=INDEX(CADENA(1:L),':')
            LLL=INDEX(CADENA(1:LL),'.')
            WRITE(33,100)'<LI><A HREF="programs/'
            IF(LLL.NE.0)THEN
              WRITE(33,100)CADENA(1:LLL-1)
            ELSE
              WRITE(33,100)CADENA(1:LL-1)
            END IF
            WRITE(33,100)'.html">'
            WRITE(33,100)CADENA(1:LL)
            WRITE(33,101)'</A>'
            WRITE(33,100)CADENA(LL+1:L)
            WRITE(33,101)'</LI>'
          ELSE
            WRITE(32,101) CADENA(1:L)
            WRITE(36,101) CADENA(1:L)
          END IF
          RETURN
        END IF
C
        IF(UNIT.EQ.0)THEN
          WRITE(33,101)
          LL=INDEX(CADENA(1:L),':')
          LLL=INDEX(CADENA(1:LL),'.')
          WRITE(33,100)'<LI><A HREF="programs/'
          IF(LLL.NE.0)THEN
            WRITE(33,100)CADENA(1:LLL-1)
          ELSE
            WRITE(33,100)CADENA(1:LL-1)
          END IF
          WRITE(33,100)'.html">'
        END IF
        DO WHILE(L.GT.LWIDTH)
          LLOOP=.TRUE.
          I=LWIDTH+1
          DO WHILE(LLOOP)
            I=I-1
            IF(I.EQ.0)THEN
              LLOOP=.FALSE.
            ELSE
              IF(CADENA(I:I).EQ.' ') LLOOP=.FALSE.
            END IF
          END DO
          IF(I.EQ.0)THEN
            WRITE(*,101)'FATAL ERROR: in subroutine PRINTDESCRIPTION.'
            STOP
          END IF
          IF(UNIT.EQ.0)THEN
            IF(LL.NE.0)THEN
              WRITE(33,100) CADENA(1:LL)
              WRITE(33,101)'</A>'
              WRITE(33,101)CADENA(LL+1:I-1)
              LL=0
            ELSE
              WRITE(33,101) CADENA(1:I-1)
            END IF
            WRITE(30,101) CADENA(1:I-1)
          ELSE
            WRITE(32,101) CADENA(1:I-1)
            WRITE(36,101) CADENA(1:I-1)
          END IF
          CADENA(1:L-I)=CADENA(I+1:L)
          L=L-I
          IF(L.LE.LWIDTH)THEN
            IF(UNIT.EQ.0)THEN
              WRITE(30,101) CADENA(1:L)
              WRITE(33,101) CADENA(1:L)
            ELSE
              WRITE(32,101) CADENA(1:L)
              WRITE(36,101) CADENA(1:L)
            END IF
          END IF
        END DO
        IF(UNIT.EQ.0)THEN
          WRITE(33,101)'</LI>'
        END IF
100     FORMAT(A,$)
101     FORMAT(A)
        END
