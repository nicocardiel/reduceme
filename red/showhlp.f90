!------------------------------------------------------------------------------
! Version 20-November-2026                                    File: showhlp.f90
!------------------------------------------------------------------------------
! Copyright N. Cardiel & J. Gorgas, Departamento de Astrofisica
! Universidad Complutense de Madrid, 28040-Madrid, Spain
! E-mail: ncl@astrax.fis.ucm.es or fjg@astrax.fis.ucm.es
!------------------------------------------------------------------------------
! This routine is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the Free
! Software Foundation; either version 2 of the License, or (at your option) any
! later version. See the file gnu-public-license.txt for details.
!------------------------------------------------------------------------------
!Comment
!
! SUBROUTINE SHOWHLP(CADENA)
!
! Input: CADENA
! Output: CADENA
!
! Show additional help in programs at running time. This routine is employed
! for maintenance purposes (i.e. the creation of the help WEB page). When
! running any REDUCEME program, this routine checks whether any of the
! following two files exist:
! .running_HLP     -> generates HELP info in the terminal
! .running_HLPHTML -> generates HELP info to create a WEB page
!
! CHARACTER*(*) CADENA -> character string to identify the piece of information
!                         to be extracted from the help file, which must be
!                         located in the directory $reduceme_dir/help/programs
!
!Comment
!------------------------------------------------------------------------------
        SUBROUTINE SHOWHLP(CADENA)
        IMPLICIT NONE
        CHARACTER*(*) CADENA
!
        INCLUDE 'redlib.inc'
        INTEGER TRUELEN
!
        INTEGER IO,I
        INTEGER LRED,LPROGRAM,LCADENA,L1,L2
        CHARACTER*255 LOCALFILE,REDUCEMEDIR,LINEA
        LOGICAL LOK1,LOK2
        LOGICAL LOGFILE,LSEARCHIO,LEXIST,LOPENED
!------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
!
        INQUIRE(FILE='.running_HLP',EXIST=LOK1)
        INQUIRE(FILE='.running_HLPHTML',EXIST=LOK2)
        IF((.NOT.LOK1).AND.(.NOT.LOK2)) RETURN
!
        CALL GET_ENVIRONMENT_VARIABLE('reduceme_dir',REDUCEMEDIR)
        LRED=TRUELEN(REDUCEMEDIR)
        IF(LRED.EQ.0) RETURN
!
        LPROGRAM=TRUELEN(THISPROGRAM)
        IF(LPROGRAM.EQ.0) RETURN
!
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
!
        INQUIRE(FILE=LOCALFILE(1:L2),EXIST=LOGFILE)
        IF(.NOT.LOGFILE) RETURN
!
! buscamos un numero de unidad de fichero no utilizada
        LSEARCHIO=.TRUE.
        IO=49
        DO WHILE(LSEARCHIO)
          IO=IO+1
          IF(IO.EQ.91) RETURN
          INQUIRE(UNIT=IO,EXIST=LEXIST,OPENED=LOPENED)
          LSEARCHIO=(.NOT.LEXIST).OR.(LOPENED)
        END DO
!
        LCADENA=TRUELEN(CADENA)
!
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
              WRITE(*,'(8X,A,A)',ADVANCE='NO')'| ',LINEA(1:L1)
              IF(L1.LT.60)THEN
                DO I=L1+1,60
                  WRITE(*,'(A)',ADVANCE='NO')' '
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
!
111     FORMAT(8X,'|',61('-'),'|')
        END
