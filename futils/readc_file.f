C------------------------------------------------------------------------------
C Version 20-November-2026                                   File: readc_file.f
C------------------------------------------------------------------------------
C Copyright N. Cardiel & J. Gorgas, Departamento de Astrofisica
C Universidad Complutense de Madrid, 28040-Madrid, Spain
C------------------------------------------------------------------------------
C This routine is free software; you can redistribute it and/or modify it
C under the terms of the GNU General Public License as published by the Free
C Software Foundation; either version 2 of the License, or (at your option) any
C later version. See the file gnu-public-license.txt for details.
C------------------------------------------------------------------------------
Comment
C
C CHARACTER*255 FUNCTION READC_FILE(CDEF,CVAL)
C
C Input: CDEF,CVAL
C Output: READC_FILE (function)
C
C Return a file name entered by the user through the keyboard. The interface
C and behaviour are the same as in READC, with the following differences:
C
C - If the package has been compiled with readline support (see configure
C   option --with-readline) and both stdin and stdout are terminals, the line
C   is read with readline: TAB completes the file name up to the longest
C   unambiguous prefix, a second TAB lists the possible completions, and
C   previous entries can be recalled with the arrow keys. Otherwise (e.g.
C   input redirected from a file) the line is read with a Fortran READ, as
C   in READC.
C - If the entered string contains "*" or "?", the matching files are listed
C   (through a call to "ls") and the user is prompted again.
C
C CHARACTER*(*) CDEF -> character with default value for READC_FILE
C               ('@' if there is no default)
C CHARACTER*(*) CVAL -> character string with valid characters
C               ('@' if all characters are valid)
C
Comment
C------------------------------------------------------------------------------
        CHARACTER*255 FUNCTION READC_FILE(CDEF,CVAL)
        USE ISO_C_BINDING, ONLY: C_CHAR,C_INT,C_NULL_CHAR
        USE ISO_FORTRAN_ENV, ONLY: OUTPUT_UNIT
        IMPLICIT NONE
        CHARACTER*(*) CDEF
        CHARACTER*(*) CVAL
C
C C routine in freadline.c
        INTERFACE
          INTEGER(C_INT) FUNCTION FREADLINE(PROMPT,BUFFER,BUFLEN)
     +     BIND(C,NAME='freadline')
          IMPORT C_CHAR,C_INT
          CHARACTER(KIND=C_CHAR) PROMPT(*)
          CHARACTER(KIND=C_CHAR) BUFFER(*)
          INTEGER(C_INT),VALUE :: BUFLEN
          END FUNCTION FREADLINE
        END INTERFACE
C
        INTEGER I,L1,L2,LP
        INTEGER TRUEBEG,TRUELEN
        INTEGER NERR
        INTEGER ISTATUS
        INTEGER ISYSTEM
        INTEGER SYSTEMFUNCTION
        CHARACTER*255 CADENA
        CHARACTER*512 CPROMPT
        LOGICAL LRUN,LMANUAL,LHTML
C------------------------------------------------------------------------------
        NERR=0
C build the prompt (same text as in READC); LP is its length, including
C the trailing blank
10      LP=0
        IF(CDEF.NE.'@')THEN
          L1=TRUEBEG(CDEF)
          IF(L1.NE.0)THEN
            L2=TRUELEN(CDEF)
            CPROMPT='['//CDEF(L1:L2)//'] ? '
            LP=L2-L1+6
          END IF
        ELSE
          CPROMPT='? '
          LP=2
          L1=0 !evita un WARNING de compilacion
        END IF
C flush pending Fortran output (e.g. a question written with FORMAT(A,$))
C before readline writes the prompt through the C standard I/O library
        FLUSH(OUTPUT_UNIT)
        ISTATUS=FREADLINE(CPROMPT(1:LP)//C_NULL_CHAR,CADENA,
     +   LEN(CADENA))
        IF(ISTATUS.EQ.0)THEN             !readline not used: as in READC
          IF(LP.GT.0) WRITE(*,100,ADVANCE='NO')CPROMPT(1:LP)
          READ(*,'(A)',ERR=20) CADENA
        ELSEIF(ISTATUS.LT.0)THEN                  !EOF (Ctrl-D) in a tty
          WRITE(*,*)
          GOTO 20
        END IF
        IF(TRUELEN(CADENA).GE.10) CALL ENDPROGRAM(CADENA)
        IF(CVAL.EQ.'@')THEN
          IF(TRUELEN(CADENA).EQ.0)THEN
            IF(CDEF.EQ.'@')THEN
              GOTO 10
            END IF
            CADENA=CDEF(L1:L2)
          END IF
        ELSE
          IF(TRUELEN(CADENA).EQ.0)THEN
            IF(CDEF.EQ.'@')THEN
              GOTO 10
            END IF
            CADENA=CDEF(L1:L2)
          ELSE
            DO I=1,TRUELEN(CADENA)
              IF(INDEX(CVAL,CADENA(I:I)).EQ.0)THEN
                WRITE(*,101)'ERROR: invalid character(s). Try again.'
                NERR=NERR+1
                IF(NERR.GT.10) STOP 'FATAL ERROR: too many errors.'
                GOTO 10
              END IF
            END DO
          END IF
        END IF
        READC_FILE=CADENA
        CALL LRUNX(LRUN,LMANUAL,LHTML)
        IF(LRUN)THEN
          CALL RMBLANK(CADENA,CADENA,L2)
          WRITE(*,101) ACHAR(27)//'[1;46m'//CADENA(1:L2)//
     +     ACHAR(27)//'[1;m'
        END IF
        IF(LMANUAL)THEN
          CALL RMBLANK(CADENA,CADENA,L2)
          WRITE(*,101) '\\ttshade{'//CADENA(1:L2)//'}'
        END IF
        IF(LHTML)THEN
          CALL RMBLANK(CADENA,CADENA,L2)
          WRITE(*,101) '<FONT COLOR="#FF0000">'//CADENA(1:L2)//'</FONT>'
        END IF
C wildcards: list the matching files and prompt again
        IF((INDEX(READC_FILE,'*').NE.0).OR.
     +   (INDEX(READC_FILE,'?').NE.0))THEN
          ISYSTEM=SYSTEMFUNCTION('ls '//
     +     READC_FILE(1:TRUELEN(READC_FILE)))
          GOTO 10
        END IF
        RETURN
20      WRITE(*,101)'ERROR: invalid entry. Try again.'
        NERR=NERR+1
        IF(NERR.GT.10) STOP 'FATAL ERROR: too many errors.'
        GOTO 10
100     FORMAT(A)
101     FORMAT(A)
        END
