!------------------------------------------------------------------------------
! Version 20-November-2026                                   File: rpgbegin.f90
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
! SUBROUTINE RPGBEGIN(NTERM,IDN,LCOLOR)
!
! Output: NTERM,IDN,LCOLOR
! Output (COMMON): all global variables in button.inc
!
! Open the graphic device(s) and assign the default values to the global
! variables:
!       MAX_XBUTT=6
!       MAX_YBUTT=2
!       PGSCF_BUTT=2
!       PGSCH_BUTT=1.
!       YTEXT_BUTT=0.35
!       X1VPORT=0.1
!       X2VPORT=0.95
!       Y1VPORT=0.1
!       Y2VPORT=0.70
!       X3VPORT=0.05
!       X4VPORT=0.95
!       Y3VPORT=0.80
!       Y4VPORT=0.95
!
! INTEGER NTERM -> number of opened graphic devices to be employed
!         simultaneously
! INTEGER IDN(8) -> identifier of the openned graphic devices
!         (positive values returned by PGOPEN)
! LOGICAL LCOLOR(8) -> determines whether color is available or not
!         in each opened graphic device
!
!
!Comment
!------------------------------------------------------------------------------
        SUBROUTINE RPGBEGIN(NTERM,IDN,LCOLOR)
        IMPLICIT NONE
!------------------------------------------------------------------------------
        INCLUDE 'button.inc'
        INTEGER TLENBUTT
        INCLUDE 'readcbutt.inc'
!------------------------------------------------------------------------------
! subroutine parameters
        INTEGER NTERM  !number of graphic devices to be employed simultaneously
        INTEGER IDN(MAX_ID_BUTT)      !identification number of graphic devices
        LOGICAL LCOLOR(MAX_ID_BUTT)      !determines whether color is available
!
! local variables
        INTEGER ID0
        INTEGER I,LDEV,LTER
        INTEGER PGOPEN
        INTEGER CI1,CI2
        CHARACTER*1 CBUT
        CHARACTER*255 TERMINAL,DEFDEV
        LOGICAL LOOP
        LOGICAL LRUN,LMANUAL,LHTML
!------------------------------------------------------------------------------
        CALL GRGENV('DEV',DEFDEV,LDEV)
! determinamos si los botones se muestran en modo texto (la opción '#' se
! emplea para utilizar los botones en modo texto pero sin preguntas en la
! terminal, para ayudar así en la generación de documentación HTML)
        LOOP=.TRUE.
        DO WHILE(LOOP)
          WRITE(*,100,ADVANCE='NO')&
     &       'Are you using graphic buttons.....(y/n) [y] ? '
          READ(*,101) CBUT
          IF(CBUT.EQ.'') CBUT='y'
          IF((CBUT.NE.'y').AND.(CBUT.NE.'n').AND.(CBUT.NE.'#'))THEN
            WRITE(*,101) 'ERROR: invalid entry. Try again.'
          ELSE
            LOOP=.FALSE.
          END IF
        END DO
!
        IF(CBUT.EQ.'#')THEN
          MODOTEXT_HTML=.TRUE.
          MODOTEXT_BUTT=.FALSE.
          MODOTEXT_PLOTBUTT=.TRUE.
          CALL LRUNXBUTT(LRUN,LMANUAL,LHTML)
          IF(LRUN)THEN
            WRITE(*,101) 'y'
          END IF
          IF(LMANUAL)THEN
            WRITE(*,101) '\\ttshade{y}'
          END IF
          IF(LHTML)THEN
            WRITE(*,101) '<FONT COLOR="#FF0000">y</FONT>'
          END IF
        ELSE
          MODOTEXT_HTML=.FALSE.
          MODOTEXT_BUTT=(CBUT.EQ.'n')
          IF(MODOTEXT_BUTT)THEN
            WRITE(*,*)
            WRITE(*,100,ADVANCE='NO')&
     &         'WARNING: remember to enlarge the terminal'
            WRITE(*,101)' to hold all the text-buttons.'
            WRITE(*,100,ADVANCE='NO')&
     &         'Do you want to plot the buttons anyway (y/n) '
            CBUT(1:1)=READCBUTT('y','yn')
            MODOTEXT_PLOTBUTT=(CBUT.EQ.'y')
            WRITE(*,'(A)')'[1;1f[J'
          ELSE
            MODOTEXT_PLOTBUTT=.FALSE.
          END IF
        END IF
! open graphic output for plots
        LOOP=.TRUE.
        NTERM=0
        DO WHILE(LOOP)
          IF(NTERM.GT.0)THEN
            WRITE(*,115,ADVANCE='NO')&
     &         'Graphic device #',NTERM+1,' (NONE=EXIT)'
          ELSE
            WRITE(*,100,ADVANCE='NO')'Graphic device #1'
          END IF
          IF(LDEV.EQ.0)THEN
            CALL PGLDEV
            IF(NTERM.GT.0)THEN
              WRITE(*,100,ADVANCE='NO')' '
              TERMINAL=READCBUTT('NONE','@')
            ELSE
              TERMINAL=READCBUTT('@','@')
            END IF
          ELSE
            WRITE(*,100,ADVANCE='NO')' (? to see list) '
            IF(NTERM.GT.0)THEN
              TERMINAL=READCBUTT('NONE','@')
            ELSE
              TERMINAL=READCBUTT(DEFDEV(1:LDEV),'@')
            END IF
            IF(TERMINAL.EQ.'?')THEN
              CALL PGLDEV
              WRITE(*,115,ADVANCE='NO')'Graphic device #',NTERM+1,' '
              TERMINAL=READCBUTT(DEFDEV(1:LDEV),'@')
            END IF
          END IF
          IF(TERMINAL.NE.'NONE')THEN
            LTER=TLENBUTT(TERMINAL)
            ID0=PGOPEN(TERMINAL(1:LTER))
            IF(ID0.LE.0)THEN
              WRITE(*,101)'ERROR: unable to open graphic device. '//&
     &         'Please, try again.'
              ID0=PGOPEN('?')
              IF(ID0.LE.0)THEN
                WRITE(*,101)'FATAL ERROR: '//&
     &           'unable to open graphic device.'
                STOP
              END IF
            END IF
            NTERM=NTERM+1
            IDN(NTERM)=ID0
            LOOP=(NTERM.LT.MAX_ID_BUTT)
            CALL PGASK(.FALSE.)
            CALL PGQCOL(CI1,CI2)
            LCOLOR(NTERM)=(CI2.GE.15)
          ELSE
            IF(NTERM.EQ.0)THEN
              WRITE(*,101)'ERROR: you must define at least one '//&
     &         'graphic device.'
            ELSE
              LOOP=.FALSE.
            END IF
          END IF
        END DO
! dejamos como ultima ventana activa la primera (si hay mas de una)
        IF(NTERM.GT.1)THEN
          CALL PGSLCT(IDN(1))
        END IF
! definimos las variables de entorno por defecto
        MAX_XBUTT=6
        MAX_YBUTT=2
        PGSCF_BUTT=2
        PGSCH_BUTT=1.
        YTEXT_BUTT=0.35
        X1VPORT=0.1
        X2VPORT=0.95
        Y1VPORT=0.1
        Y2VPORT=0.70
        X3VPORT=0.05
        X4VPORT=0.95
        Y3VPORT=0.80
        Y4VPORT=0.95
! ningun boton ha sido dibujado todavia
        DO I=1,MAX_NBUTT
          EXIST_BUTT(I)=.FALSE.
        END DO
! los ticks se dibujan hacia dentro
        ITICKS_BUTT=.FALSE.
!
! Redefinimos colores (salida grafica con grises)
!   0: negro
!   1: blanco
!  12: gris intermedio
!  13: gris claro
!  14: gris oscuro
!  15: casi gris intermedio
        CALL PGSCR(12,0.5,0.5,0.5)
        CALL PGSCR(13,0.7,0.7,0.7)
        CALL PGSCR(14,0.3,0.3,0.3)
        CALL PGSCR(15,0.6,0.6,0.6)
! definimos un viewport y un window por defecto
!cc     CALL PGVPORT(X1VPORT,X2VPORT,Y1VPORT,Y2VPORT)
!cc     CALL PGWINDOW(0.,1.,0.,1.)
!cc     CALL PGPAGEB_RED
!cc     CALL PGPAGE_RED
!------------------------------------------------------------------------------
100     FORMAT(A)
101     FORMAT(A)
115     FORMAT(A,I1,A)
        END
