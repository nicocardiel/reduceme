!------------------------------------------------------------------------------
! Version 20-November-2026                                      File: readf.f90
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
! REAL FUNCTION READF(CDEF)
!
! Input: CDEF
! Output: READF (function)
!
! Return a float number entered by the user through the keyboard.
!
! CHARACTER*(*) CDEF -> character string with default value for READF
!               ('@' if there is no default)
!
!Comment
!------------------------------------------------------------------------------
        REAL FUNCTION READF(CDEF)
        IMPLICIT NONE
        CHARACTER*(*) CDEF
!
        INTEGER I,L1,L2
        INTEGER NERR
        REAL F
        INTEGER TRUEBEG,TRUELEN
        CHARACTER*1 C
        CHARACTER*255 CDUMMY
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
        END IF
        READ(*,'(A)',ERR=20)CADENA
        IF(TRUELEN(CADENA).GE.10) CALL ENDPROGRAM(CADENA)
        IF(TRUELEN(CADENA).EQ.0)THEN
          IF(CDEF.EQ.'@')THEN
            GOTO 10
          END IF
          CADENA=CDEF
        END IF
        DO I=1,TRUELEN(CADENA)
          C=CADENA(I:I)
          IF((INDEX('abcfghijklmnoprstuvwxyz',C).NE.0).OR.&
     &     (INDEX('ABCFGHIJKLMNOPRSTUVWXYZ/',C).NE.0))THEN
            GOTO 20
          END IF
        END DO
        READ(CADENA,*,ERR=20) F
        READF=F
        CALL LRUNX(LRUN,LMANUAL,LHTML)
        IF(LRUN)THEN
          WRITE(CDUMMY,*) F
          CALL RMBLANK(CDUMMY,CDUMMY,L2)
          !WRITE(*,101) CDUMMY(1:L2)
          WRITE(*,101) ACHAR(27)//'[1;46m'//CDUMMY(1:L2)//&
     &     ACHAR(27)//'[1;m'
        END IF
        IF(LMANUAL)THEN
          WRITE(CDUMMY,*) F
          CALL RMBLANK(CDUMMY,CDUMMY,L2)
          WRITE(*,101) '\\ttshade{'//CDUMMY(1:L2)//'}'
        END IF
        IF(LHTML)THEN
          WRITE(CDUMMY,*) F
          CALL RMBLANK(CDUMMY,CDUMMY,L2)
          WRITE(*,101) '<FONT COLOR="#FF0000">'//CDUMMY(1:L2)//'</FONT>'
        END IF
        RETURN
20      WRITE(*,101)'ERROR: invalid character(s) found in '//&
     &   'number. Try again.'
        IF(CDEF.EQ.'@') WRITE(*,100,ADVANCE='NO')'? '
        NERR=NERR+1
        IF(NERR.GT.10) STOP 'FATAL ERROR: too many errors.'
        GOTO 10
100     FORMAT(A)
101     FORMAT(A)
        END
