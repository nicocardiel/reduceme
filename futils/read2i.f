C------------------------------------------------------------------------------
C Version 3-March-2006                                           File: read2i.f
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
C SUBROUTINE READ2I(CDEF,N1,N2)
C
C Input: CDEF
C Output: N1,N2
C
C Return 2 integers N1 and N2 entered by the user through the keyboard.
C
C CHARACTER*(*) CDEF -> character string with default values for N1 and N2 
C               ('@' if there is no default)
C INTEGER N1 -> first integer
C INTEGER N2 -> second integer
C
Comment
C------------------------------------------------------------------------------
        SUBROUTINE READ2I(CDEF,N1,N2)
        IMPLICIT NONE
        INTEGER N1,N2
        CHARACTER*(*) CDEF
C
        INTEGER I,L1,L2
        INTEGER LCADENA
        INTEGER TRUEBEG,TRUELEN
        INTEGER NERR
        CHARACTER*1 C
        CHARACTER*255 CDUMMY
        CHARACTER*255 CADENA
        LOGICAL LRUN,LMANUAL,LHTML
C------------------------------------------------------------------------------
        NERR=0
10      IF(CDEF.NE.'@')THEN
          L1=TRUEBEG(CDEF)
          IF(L1.NE.0)THEN
            L2=TRUELEN(CDEF)
            WRITE(*,100)'['
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
        LCADENA=TRUELEN(CADENA)
        IF((INDEX(CADENA(1:LCADENA),',').EQ.0).AND.
     +   (INDEX(CADENA(1:LCADENA),' ').EQ.0))THEN
          WRITE(*,101) 'ERROR: only one possible number found. '//
     +     'Try again.'
          NERR=NERR+1
          IF(NERR.GT.10) STOP 'FATAL ERROR: too many errors.'
          GOTO 10
        END IF
        READ(CADENA,*,ERR=20) N1,N2
        CALL LRUNX(LRUN,LMANUAL,LHTML)
        IF(LRUN)THEN
          WRITE(CDUMMY,'(I10,A1,I10)')N1,',',N2
          CALL RMBLANK(CDUMMY,CDUMMY,L2)
          WRITE(*,101) CDUMMY(1:L2)
        END IF
        IF(LMANUAL)THEN
          WRITE(CDUMMY,'(I10,A1,I10)')N1,',',N2
          CALL RMBLANK(CDUMMY,CDUMMY,L2)
          WRITE(*,101) '\\ttshade{'//CDUMMY(1:L2)//'}'
        END IF
        IF(LHTML)THEN
          WRITE(CDUMMY,'(I10,A1,I10)')N1,',',N2
          CALL RMBLANK(CDUMMY,CDUMMY,L2)
          WRITE(*,101) '<FONT COLOR="#FF0000">'//CDUMMY(1:L2)//'</FONT>'
        END IF
        RETURN
20      WRITE(*,101)'ERROR: invalid character(s) found in '//
     +   'numbers. Try again.'
        NERR=NERR+1
        IF(NERR.GT.10) STOP 'FATAL ERROR: too many errors.'
        GOTO 10
100     FORMAT(A,$)
101     FORMAT(A)
        END
