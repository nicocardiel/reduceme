C------------------------------------------------------------------------------
C Version 21-August-2004                                       File: readilim.f
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
C INTEGER FUNCTION READILIM(CDEF,N1,N2)
C
C Input: CDEF,N1,N2
C Output: READILIM (function)
C
C Return an integer number entered by the user through the keyboard in the
C range from N1 to N2
C
C CHARACTER*(*) CDEF -> character string with default value for READILIM
C               ('@' if there is no default)
C INTEGER       N1 -> first limit for READILIM
C INTEGER       N2 -> second limit for READILIM
C
Comment
C------------------------------------------------------------------------------
        INTEGER FUNCTION READILIM(CDEF,N1,N2)
        IMPLICIT NONE
        CHARACTER*(*) CDEF
        INTEGER N1,N2
C
        INTEGER I,L1,L2
        INTEGER N
        INTEGER NERR
        INTEGER TRUEBEG,TRUELEN
        CHARACTER*1 C
        CHARACTER*255 CDUMMY
        CHARACTER*255 CADENA
        LOGICAL LRUN,LMANUAL,LHTML
C------------------------------------------------------------------------------
        IF(N2.LT.N1)THEN
          WRITE(*,101)'ERROR: N2.LT.N1 in function: READILIM'
          WRITE(*,101)'=> returned value is 0'
          READILIM=0
          RETURN
        END IF
        NERR=0
        WRITE(CDUMMY,'(A1,I10,A5,I10,A1)') '(',N1,',...,',N2,')'
        CALL RMBLANK(CDUMMY,CDUMMY,L2)
        WRITE(*,100) ' '//CDUMMY(1:L2)
10      IF(CDEF.NE.'@')THEN
          L1=TRUEBEG(CDEF)
          IF(L1.NE.0)THEN
            L2=TRUELEN(CDEF)
            WRITE(*,100)' ['
            WRITE(*,100)CDEF(L1:L2)
            WRITE(*,100)'] ? '
          END IF
        ELSE
          WRITE(*,100)'? '
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
          IF((INDEX('abcdefghijklmnopqrstuvwxyz',C).NE.0).OR.
     +     (INDEX('ABCDEFGHIJKLMNOPQRSTUVWXYZ./',C).NE.0))THEN
            GOTO 20
          END IF
        END DO
        READ(CADENA,*,ERR=20) N
        READILIM=N
C
        IF((N.LT.N1).OR.(N.GT.N2)) GOTO 30
        CALL LRUNX(LRUN,LMANUAL,LHTML)
        IF(LRUN)THEN
          WRITE(CDUMMY,*) N
          CALL RMBLANK(CDUMMY,CDUMMY,L2)
          !WRITE(*,101) CDUMMY(1:L2)
          WRITE(*,101) ACHAR(27)//'[1;31m'//CDUMMY(1:L2)//
     +     ACHAR(27)//'[1;m'
        END IF
        IF(LMANUAL)THEN
          WRITE(CDUMMY,*) N
          CALL RMBLANK(CDUMMY,CDUMMY,L2)
          WRITE(*,101) '\\ttshade{'//CDUMMY(1:L2)//'}'
        END IF
        IF(LHTML)THEN
          WRITE(CDUMMY,*) N
          CALL RMBLANK(CDUMMY,CDUMMY,L2)
          WRITE(*,101) '<FONT COLOR="#FF0000">'//CDUMMY(1:L2)//'</FONT>'
        END IF
        RETURN
C------------------------------------------------------------------------------
20      WRITE(*,101)'ERROR: invalid character(s) found in '//
     +   'number. Try again.'
C
        NERR=NERR+1
        IF(NERR.GT.10) STOP 'FATAL ERROR: too many errors.'
        GOTO 10
C------------------------------------------------------------------------------
30      WRITE(*,100)'ERROR: invalid number. Valid range is '
        WRITE(CDUMMY,*)N1
        CALL RMBLANK(CDUMMY,CDUMMY,L2)
        WRITE(*,100) CDUMMY(1:L2)//' to '
        WRITE(CDUMMY,*)N2
        CALL RMBLANK(CDUMMY,CDUMMY,L2)
        WRITE(*,101) CDUMMY(1:L2)//'. Try again.'
C
        NERR=NERR+1
        IF(NERR.GT.10) STOP 'FATAL ERROR: too many errors.'
        GOTO 10
100     FORMAT(A,$)
101     FORMAT(A)
        END
