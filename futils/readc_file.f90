!------------------------------------------------------------------------------
! Version 20-November-2026                                   File: readc_file.f
!------------------------------------------------------------------------------
! Copyright N. Cardiel & J. Gorgas, Departamento de Astrofisica
! Universidad Complutense de Madrid, 28040-Madrid, Spain
!------------------------------------------------------------------------------
! This routine is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the Free
! Software Foundation; either version 2 of the License, or (at your option) any
! later version. See the file gnu-public-license.txt for details.
!------------------------------------------------------------------------------
!Comment
!
! CHARACTER*255 FUNCTION READC_FILE(CDEF,CVAL)
!
! Input: CDEF,CVAL
! Output: READC_FILE (function)
!
! Return a file name entered by the user through the keyboard. The interface
! and behaviour are the same as in READC, with the following differences:
!
! - If the package has been compiled with readline support (see configure
!   option --with-readline) and both stdin and stdout are terminals, the line
!   is read with readline: TAB completes the file name up to the longest
!   unambiguous prefix, a second TAB lists the possible completions, and
!   previous entries can be recalled with the arrow keys. Otherwise (e.g.
!   input redirected from a file) the line is read with a Fortran READ, as
!   in READC.
! - If the entered string contains "*" or "?", the matching files are listed
!   (through a call to "ls") and the user is prompted again.
!
! CHARACTER*(*) CDEF -> character with default value for READC_FILE
!               ('@' if there is no default)
! CHARACTER*(*) CVAL -> character string with valid characters
!               ('@' if all characters are valid)
!
!Comment
!------------------------------------------------------------------------------
        CHARACTER*255 FUNCTION READC_FILE(CDEF,CVAL)
        USE ISO_C_BINDING, ONLY: C_CHAR,C_INT,C_NULL_CHAR
        USE ISO_FORTRAN_ENV, ONLY: OUTPUT_UNIT
        IMPLICIT NONE
        CHARACTER*(*) CDEF
        CHARACTER*(*) CVAL
!
! C routine in freadline.c
        INTERFACE
          INTEGER(C_INT) FUNCTION FREADLINE(PROMPT,BUFFER,BUFLEN)&
     &     BIND(C,NAME='freadline')
          IMPORT C_CHAR,C_INT
          CHARACTER(KIND=C_CHAR) PROMPT(*)
          CHARACTER(KIND=C_CHAR) BUFFER(*)
          INTEGER(C_INT),VALUE :: BUFLEN
          END FUNCTION FREADLINE
        END INTERFACE
!
        INTEGER I,L1,L2,LP
        INTEGER TRUEBEG,TRUELEN
        INTEGER NERR
        INTEGER ISTATUS
        INTEGER ISYSTEM
        INTEGER SYSTEMFUNCTION
        CHARACTER*255 CADENA
        CHARACTER*512 CPROMPT
        LOGICAL LRUN,LMANUAL,LHTML
!------------------------------------------------------------------------------
        NERR=0
! build the prompt (same text as in READC); LP is its length, including
! the trailing blank
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
! flush pending Fortran output (e.g. a question written with FORMAT(A,$))
! before readline writes the prompt through the C standard I/O library
        FLUSH(OUTPUT_UNIT)
        ISTATUS=FREADLINE(CPROMPT(1:LP)//C_NULL_CHAR,CADENA,&
     &   LEN(CADENA))
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
          WRITE(*,101) ACHAR(27)//'[1;46m'//CADENA(1:L2)//&
     &     ACHAR(27)//'[1;m'
        END IF
        IF(LMANUAL)THEN
          CALL RMBLANK(CADENA,CADENA,L2)
          WRITE(*,101) '\\ttshade{'//CADENA(1:L2)//'}'
        END IF
        IF(LHTML)THEN
          CALL RMBLANK(CADENA,CADENA,L2)
          WRITE(*,101) '<FONT COLOR="#FF0000">'//CADENA(1:L2)//'</FONT>'
        END IF
! wildcards: list the matching files and prompt again
        IF((INDEX(READC_FILE,'*').NE.0).OR.&
     &   (INDEX(READC_FILE,'?').NE.0))THEN
          ISYSTEM=SYSTEMFUNCTION('ls '//&
     &     READC_FILE(1:TRUELEN(READC_FILE)))
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
