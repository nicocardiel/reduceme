C------------------------------------------------------------------------------
C Version 28-July-1998                                          File: showhlp.f
C------------------------------------------------------------------------------
C Copyright N. Cardiel & J. Gorgas, Departamento de Astrofisica
C Universidad Complutense de Madrid, 28040-Madrid, Spain
C E-mail: ncl@astrax.fis.ucm.es or fjg@astrax.fis.ucm.es
C------------------------------------------------------------------------------
C This routine is free software; you can redistribute it and/or modify it
C under the terms of the GNU General Public License as published by the Free
C Software Foundation; either version 2 of the License, or (at your option) any
C later version. See the file gnu-public-license.txt for details.
C------------------------------------------------------------------------------
Comment
C
C SUBROUTINE SHOWHLP(CADENA)
C
C Input: CADENA
C Output: CADENA
C
C Show additional help in programs at running time. This routine is employed
C for maintenance purposes (i.e. the creation of the help WEB page). When
C running any REDUCEME program, this routine checks whether any of the
C following two files exist:
C .running_HLP     -> generates HELP info in the terminal
C .running_HLPHTML -> generates HELP info to create a WEB page
C
C CHARACTER*(*) CADENA -> character string to identify the piece of information
C                         to be extracted from the help file, which must be
C                         located in the directory $reduceme_dir/help/programs
C
Comment
C------------------------------------------------------------------------------
        SUBROUTINE SHOWHLP(CADENA)
        IMPLICIT NONE
        CHARACTER*(*) CADENA
C
        INCLUDE 'redlib.inc'
        INTEGER TRUELEN
C
        INTEGER IO,I
        INTEGER LRED,LPROGRAM,LCADENA,L1,L2
        CHARACTER*255 LOCALFILE,REDUCEMEDIR,LINEA
        LOGICAL LOK1,LOK2
        LOGICAL LOGFILE,LSEARCHIO,LEXIST,LOPENED
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        INQUIRE(FILE='.running_HLP',EXIST=LOK1)
        INQUIRE(FILE='.running_HLPHTML',EXIST=LOK2)
        IF((.NOT.LOK1).AND.(.NOT.LOK2)) RETURN
C
        CALL GETENV('reduceme_dir',REDUCEMEDIR)
        LRED=TRUELEN(REDUCEMEDIR)
        IF(LRED.EQ.0) RETURN
C
        LPROGRAM=TRUELEN(THISPROGRAM)
        IF(LPROGRAM.EQ.0) RETURN
C
        L1=1
        L2=LRED
        LOCALFILE(L1:L2)=REDUCEMEDIR(1:LRED)
        L1=L2+1
        L2=L2+15
        LOCALFILE(L1:L2)='/help/programs/'
        L1=L2+1
        L2=L2+LPROGRAM
        LOCALFILE(L1:L2)=THISPROGRAM(1:LPROGRAM)
        L1=L2+1
        L2=L2+4
        LOCALFILE(L1:L2)='.hlp'
C
        INQUIRE(FILE=LOCALFILE(1:L2),EXIST=LOGFILE)
        IF(.NOT.LOGFILE) RETURN
C
C buscamos un numero de unidad de fichero no utilizada
        LSEARCHIO=.TRUE.
        IO=49
        DO WHILE(LSEARCHIO)
          IO=IO+1
          IF(IO.EQ.91) RETURN
          INQUIRE(UNIT=IO,EXIST=LEXIST,OPENED=LOPENED)
          LSEARCHIO=(.NOT.LEXIST).OR.(LOPENED)
        END DO
C
        LCADENA=TRUELEN(CADENA)
C
        OPEN(IO,FILE=LOCALFILE(1:L2),STATUS='OLD',FORM='FORMATTED')
10      READ(IO,'(A)',END=20) LINEA
        IF(LINEA(1:5).EQ.'#KEY:')THEN
          L2=TRUELEN(LINEA)
          IF(LINEA(7:L2).EQ.CADENA(1:LCADENA))THEN
            IF(LOK1)THEN
              WRITE(*,111)
            ELSE
              IF(CADENA(1:LCADENA).EQ.'explanation')THEN
                WRITE(*,'(A)')'<BLOCKQUOTE><I><B><FONT COLOR="#006666">'
              ELSE
                WRITE(*,'(A)')'<BLOCKQUOTE><I><FONT COLOR="#000066">'
              END IF
            END IF
12          READ(IO,'(A)',END=20) LINEA
            IF(LINEA(1:5).EQ.'#KEY:') GOTO 20
            L1=TRUELEN(LINEA)
            IF(LOK1)THEN
              WRITE(*,'(8X,A,A,$)')'| ',LINEA(1:L1)
              IF(L1.LT.60)THEN
                DO I=L1+1,60
                  WRITE(*,'(A,$)')' '
                END DO
              END IF
              WRITE(*,'(A)')'|'
            ELSE
              WRITE(*,'(A)')LINEA(1:L1)
            END IF
            GOTO 12
          END IF
        END IF
        GOTO 10
20      CLOSE(IO)
        IF(LOK1)THEN
          WRITE(*,111)
        ELSE
          IF(CADENA(1:LCADENA).EQ.'explanation')THEN
            WRITE(*,'(A)')'</FONT></B></I></BLOCKQUOTE>'
          ELSE
            WRITE(*,'(A)')'</FONT></I></BLOCKQUOTE>'
          END IF
        END IF
C
111     FORMAT(8X,'|',61('-'),'|')
        END
