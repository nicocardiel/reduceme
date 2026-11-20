!------------------------------------------------------------------------------
! Version 20-November-2026                                 File: endprogram.f90
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
! SUBROUTINE ENDPROGRAM(CADENA)
!
! Input: CADENA
!
! Multipurpose routine: stops the program, inserts CALL PGPAGE or 
! captures XServe.
!
! CHARACTER*255 CADENA -> if CADENA='endprogram' CALL PGEND+STOP
!                         if CADENA='newpagenew' CALL PGPAGE
!                         if CADENA='capturegif' capture current X11 image
!
!Comment
!------------------------------------------------------------------------------
        SUBROUTINE ENDPROGRAM(CADENA)
        IMPLICIT NONE
!
        INTEGER TRUELEN
        INTEGER SYSTEMFUNCTION
!
        INTEGER L,L1,L2
        INTEGER ISYSTEM
        CHARACTER*255 CADENA
        CHARACTER*1 COPC
        LOGICAL LRUN,LMANUAL,LHTML
!------------------------------------------------------------------------------
        CALL LRUNX(LRUN,LMANUAL,LHTML)
!------------------------------------------------------------------------------
        IF(CADENA(1:10).EQ.'endprogram')THEN
          IF(LRUN)THEN
            !WRITE(*,101) 'endprogram'
            WRITE(*,101) ACHAR(27)//'[1;31m'//'endprogram'//&
     &       ACHAR(27)//'[1;m'
          END IF
          IF(LMANUAL)THEN
            WRITE(*,101) '\\ttshade{endprogram}'
          END IF
          IF(LHTML)THEN
            WRITE(*,101) '<FONT COLOR="#FF0000">endprogram</FONT>'
          END IF
          WRITE(*,100,ADVANCE='NO')&
     &       'Do you want to STOP the program (y/n) [n] ? '
          READ(*,'(A)') COPC
          IF(LRUN)THEN
            !WRITE(*,101)COPC
            WRITE(*,101) ACHAR(27)//'[1;31m'//COPC//&
     &       ACHAR(27)//'[1;m'
          END IF
          IF(LMANUAL)THEN
            WRITE(*,101) '\\ttshade{'//COPC//'}'
          END IF
          IF(LHTML)THEN
            WRITE(*,101) '<FONT COLOR="#FF0000">'//COPC//'</FONT>'
          END IF
          IF(COPC.EQ.'y')THEN
            CALL PGEND
            STOP
          END IF
!------------------------------------------------------------------------------
        ELSEIF(CADENA(1:10).EQ.'newpagenew')THEN
          IF(LRUN)THEN
            !WRITE(*,101) 'newpagenew'
            WRITE(*,101) ACHAR(27)//'[1;31m'//'newpagenew'//&
     &       ACHAR(27)//'[1;m'
          END IF
          IF(LMANUAL)THEN
            WRITE(*,101) '\\ttshade{newpagenew}'
          END IF
          IF(LHTML)THEN
            WRITE(*,101) '<FONT COLOR="#FF0000">newpagenew</FONT>'
          END IF
          WRITE(*,*)
          WRITE(*,100,ADVANCE='NO')&
     &       'Are you inserting a newpage call (y/n) [n] ? '
          READ(*,'(A)') COPC
          IF(LRUN)THEN
            !WRITE(*,101)COPC
            WRITE(*,101) ACHAR(27)//'[1;31m'//COPC//&
     &       ACHAR(27)//'[1;m'
          END IF
          IF(LMANUAL)THEN
            WRITE(*,101) '\\ttshade{'//COPC//'}'
          END IF
          IF(LHTML)THEN
            WRITE(*,101) '<FONT COLOR="#FF0000">'//COPC//'</FONT>'
          END IF
          IF(COPC.EQ.'y')THEN
            CALL PGPAGE
          END IF
!------------------------------------------------------------------------------
! Capturamos en un fichero GIF la imagen actual en la ventana de PGPLOT.
! La variable cadena debe ser algo como "capturegif#filegif#something_else",
! donde filegif es el nombre del fichero GIF, y something_else es lo que
! finalmente va a retornar esta función a través de la variable CADENA (para
! que la función readi, readf, readc, o la que sea, continue trabajando como
! si no hubiera pasado nada. Si, además, existe el fichero .running_HLPHTML,
! se muestra en pantalla el codigo HTML necesario para insertar la figura
! GIF correspondiente.
        ELSEIF(CADENA(1:10).EQ.'capturegif')THEN
          L=TRUELEN(CADENA)
          IF(L.LE.10) RETURN
          L1=INDEX(CADENA,'#')
          IF(L1.EQ.0) RETURN                           !no encontramos primer #
          IF(L1.EQ.L) RETURN                                     !solo hay un #
          L2=INDEX(CADENA(L1+1:),'#')
          IF(L2.EQ.0) RETURN                          !no encontramos segundo #
          IF(L2.EQ.1) RETURN                      !no hay nombre de fichero gif
          L2=L1+L2
          !capturamos la ventana gráfica y generamos el fichero GIF
          ISYSTEM=SYSTEMFUNCTION('xwd -out '//CADENA(L1+1:L2-1)//'.tmp '&
     &     //'-name "PGPLOT Window 1"\0')
          ISYSTEM=SYSTEMFUNCTION('convert '//CADENA(L1+1:L2-1)//'.tmp '&
     &     //CADENA(L1+1:L2-1)//'.gif\0')
          ISYSTEM=SYSTEMFUNCTION('\\rm '//CADENA(L1+1:L2-1)//'.tmp\0')
          !mostramos en pantalla el código HTML para insertar la figura GIF
          IF(LHTML)THEN
            WRITE(*,101) '<center>'
            WRITE(*,101) '<img SRC="'//CADENA(L1+1:L2-1)//&
     &       '.gif" BORDER="2">'
            WRITE(*,101) '</center>'
          END IF
          CADENA=CADENA(L2+1:)               !devolvemos solo lo que hace falta
!------------------------------------------------------------------------------
        END IF
!------------------------------------------------------------------------------
100     FORMAT(A)
101     FORMAT(A)
        END
