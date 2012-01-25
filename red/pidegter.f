C------------------------------------------------------------------------------
C Version 26-November-1996                                     File: pidegter.f
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
C SUBROUTINE PIDEGTER(NTERM,IDN,LCOLOR)
C
C Output: NTERM,IDN,LCOLOR
C
C Open the graphic device(s), detecting whether color is available.
C
C INTEGER NTERM -> No. of opened graphic devices (maximum = MAX_ID_RED)
C INTEGER IDN(MAX_ID_RED) -> logical device number associated to each NTERM
C LOGICAL LCOLOR(MAX_ID_RED) -> .TRUE. if color is available
C
Comment
C------------------------------------------------------------------------------
        SUBROUTINE PIDEGTER(NTERM,IDN,LCOLOR)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'futils.inc'
        INTEGER TRUELEN
        INTEGER NTERM  !number of graphic devices to be employed simultaneously
        INTEGER IDN(MAX_ID_RED)
        LOGICAL LCOLOR(MAX_ID_RED)                  !='y' if color is available
C
        INTEGER ID0
        INTEGER LDEV,LTER
        INTEGER PGOPEN
        INTEGER CI1,CI2
        CHARACTER*255 TERMINAL,DEFDEV
        LOGICAL LOOP
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        CALL GRGENV('DEV',DEFDEV,LDEV)
C
        LOOP=.TRUE.
        NTERM=0
        DO WHILE(LOOP)
          IF(NTERM.GT.0)THEN
            WRITE(*,115)'Graphic device #',NTERM+1,' (NONE=EXIT)'
          ELSE
            WRITE(*,100)'Graphic device #1'
          END IF
          IF(LDEV.EQ.0)THEN
            CALL PGLDEV
            IF(NTERM.GT.0)THEN
              WRITE(*,100)' '
              TERMINAL=READC('NONE','@')
            ELSE
              TERMINAL=READC('@','@')
            END IF
          ELSE
            WRITE(*,100)' (? to see list) '
            IF(NTERM.GT.0)THEN
              TERMINAL=READC('NONE','@')
            ELSE
              TERMINAL=READC(DEFDEV(1:LDEV),'@')
            END IF
            IF(TERMINAL.EQ.'?')THEN
              CALL PGLDEV
              WRITE(*,115)'Graphic device #',NTERM+1,' '
              TERMINAL=READC(DEFDEV(1:LDEV),'@')
            END IF
          END IF
          IF(TERMINAL.NE.'NONE')THEN
            LTER=TRUELEN(TERMINAL)
            ID0=PGOPEN(TERMINAL(1:LTER))
            IF(ID0.LE.0)THEN
              WRITE(*,101)'ERROR: unable to open graphic device. '//
     +         'Please, try again.'
              ID0=PGOPEN('?')
              IF(ID0.LE.0)THEN
                WRITE(*,101)'FATAL ERROR: unable to open graphic '//
     +           'device.'
                STOP
              END IF
            END IF
            NTERM=NTERM+1
            IDN(NTERM)=ID0
            LOOP=(NTERM.LT.MAX_ID_RED)
            CALL PGQCOL(CI1,CI2)
            LCOLOR(NTERM)=(CI2.GE.15)
            CALL PGASK(.FALSE.)
          ELSE
            IF(NTERM.EQ.0)THEN
              WRITE(*,101)'ERROR: you must define at least one '//
     +         'graphic device.'
            ELSE
              LOOP=.FALSE.
            END IF
          END IF
        END DO
C dejamos como ultima ventana activa la primera (si hay mas de una)
        IF(NTERM.GT.1)THEN
          CALL PGSLCT(IDN(1))
        END IF
C------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
115     FORMAT(A,I1,A,$)
        END
