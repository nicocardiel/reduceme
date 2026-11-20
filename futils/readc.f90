!------------------------------------------------------------------------------
! Version 20-November-2026                                      File: readc.f90
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
! CHARACTER*(*) FUNCTION READC(CDEF,CVAL)
!
! Input: CDEF,CVAL
! Output: READC (function)
!
! Return a character string entered by the user through the keyboard.
!
! CHARACTER*(*) CDEF -> character with default value for READC
!               ('@' if there is no default)
! CHARACTER*(*) CVAL -> character string with valid characters
!               ('@' if all characters are valid)
!
!Comment
!------------------------------------------------------------------------------
        CHARACTER*(*) FUNCTION READC(CDEF,CVAL)
        IMPLICIT NONE
        CHARACTER*(*) CDEF
        CHARACTER*(*) CVAL
!
        INTEGER I,L1,L2
        INTEGER TRUEBEG,TRUELEN
        INTEGER NERR
        CHARACTER*255 CADENA
        LOGICAL LRUN,LMANUAL,LHTML
!------------------------------------------------------------------------------
        NERR=0
10      IF(CDEF.NE.'@')THEN
          L1=TRUEBEG(CDEF)
          IF(L1.NE.0)THEN
            L2=TRUELEN(CDEF)
            WRITE(*,100,ADVANCE='NO')'['
            WRITE(*,100,ADVANCE='NO')CDEF(L1:L2)
            WRITE(*,100,ADVANCE='NO')'] ? '
          END IF
        ELSE
          WRITE(*,100,ADVANCE='NO')'? '
          L1=0 !evita un WARNING de compilacion
        END IF
        READ(*,'(A)',ERR=20) CADENA
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
                IF(CDEF.EQ.'@') WRITE(*,100,ADVANCE='NO')'? '
                NERR=NERR+1
                IF(NERR.GT.10) STOP 'FATAL ERROR: too many errors.'
                GOTO 10
              END IF
            END DO
          END IF
        END IF
        READC=CADENA
        CALL LRUNX(LRUN,LMANUAL,LHTML)
        IF(LRUN)THEN
          CALL RMBLANK(CADENA,CADENA,L2)
          !WRITE(*,101)CADENA(1:L2)
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
        RETURN
20      WRITE(*,101)'ERROR: invalid entry. Try again.'
        IF(CDEF.EQ.'@') WRITE(*,100,ADVANCE='NO')'? '
        NERR=NERR+1
        IF(NERR.GT.10) STOP 'FATAL ERROR: too many errors.'
        GOTO 10
100     FORMAT(A)
101     FORMAT(A)
        END
