!------------------------------------------------------------------------------
! Version 20-November-2026                                      File: rpgband.f
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
! SUBROUTINE RPGBAND(MODE,POSN,XREF,YREF,XC,YC,CH)
!
! Input: MODE,POSN,XREF,YREF
! Output: XC,YC,CH
!
! This routine is similar to PGBAND, but it also allows the utilization of
! buttons in text mode.
!
! INTEGER     MODE -> display mode (see PGPLOT manual)
! INTEGER     POSN -> if POSN=1, the routine positions the cursor at the
!             position specified by XREF,YREF
! REAL        XREF -> reference position
! REAL        YREF -> reference position
! REAL        XC -> the world x-coordinate of the cursor
! REAL        YC -> the world y-coordinate of the cursor
! CHARACTER*1 CH -> the character typed by the user
!
!Comment
!------------------------------------------------------------------------------
        SUBROUTINE RPGBAND(MODE,POSN,XREF,YREF,XC,YC,CH)
        IMPLICIT NONE
        INTEGER MODE,POSN
        REAL XREF,YREF,XC,YC
        CHARACTER*1 CH
!
! Comando identico a PGCURSE de PGPLOT pero que en caso de utilizar botones en
! modo texto permite la entrada de los parametros de la subrutina desde
! teclado.
!------------------------------------------------------------------------------
        INCLUDE 'button.inc'
        REAL READFBUTT
        INCLUDE 'readcbutt.inc'
!------------------------------------------------------------------------------
! Variables locales
        INTEGER PGBAND,IDUM
        CHARACTER*1 CBUTTON,CPGPAGE
        CHARACTER*255 CADENA
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
        IF(MODOTEXT_BUTT.OR.MODOTEXT_HTML)THEN
          IF(MODOTEXT_BUTT)THEN
            WRITE(*,100,ADVANCE='NO')&
     &         'Are you pressing a button (y/n/#) '
            CBUTTON(1:1)=READCBUTT('y','yn#')
          ELSE
            READ(*,101) CADENA
            CALL ENDPROGRAM(CADENA)
            CBUTTON=CADENA(1:1)
            IF(CBUTTON.EQ.'') CBUTTON='y'
          END IF
          IF(CBUTTON.EQ.'y')THEN
            XC=0.
            YC=0.
            CH='A'         !caracter devuelto por el boton izquierdo del mouse
          ELSEIF(CBUTTON.EQ.'n')THEN
            IF(MODOTEXT_BUTT)THEN
              WRITE(*,100,ADVANCE='NO')'XC'
              XC=READFBUTT('@')
              WRITE(*,100,ADVANCE='NO')'YC'
              YC=READFBUTT('@')
              WRITE(*,100,ADVANCE='NO')'CH '
              CH(1:1)=READCBUTT('A','@')
            ELSE
              READ(*,*) XC
              READ(*,*) YC
              READ(*,101) CH
            END IF
          ELSEIF(CBUTTON.EQ.'#')THEN
            XC=0.
            YC=0.
            CH='A'         !caracter devuelto por el boton izquierdo del mouse
            WRITE(*,100,ADVANCE='NO')&
     &         'Do you want to insert a CALL PGPAGE (y/n) '
            CPGPAGE(1:1)=READCBUTT('y','yn')
            IF(CPGPAGE.EQ.'y')CALL PGPAGE
          END IF
        ELSE
          IDUM=PGBAND(MODE,POSN,XREF,YREF,XC,YC,CH)
          IF(IDUM.EQ.0)THEN
            WRITE(*,101)'ERROR: in subroutine PGBAND.'
            WRITE(*,101)'=> returned XC=YC=0.0, CH="X"'
            XC=0.
            YC=0.
            CH='X'
          END IF
        END IF
!------------------------------------------------------------------------------
100     FORMAT(A)
101     FORMAT(A)
        END
