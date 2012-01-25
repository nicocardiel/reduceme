C------------------------------------------------------------------------------
C Version 07-September-2007                                  File: endprogram.f
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
C SUBROUTINE ENDPROGRAM(CADENA)
C
C Input: CADENA
C
C Multipurpose routine: stops the program, inserts CALL PGPAGE or 
C captures XServe.
C
C CHARACTER*255 CADENA -> if CADENA='endprogram' CALL PGEND+STOP
C                         if CADENA='newpagenew' CALL PGPAGE
C                         if CADENA='capturegif' capture current X11 image
C
Comment
C------------------------------------------------------------------------------
        SUBROUTINE ENDPROGRAM(CADENA)
        IMPLICIT NONE
C
        INTEGER TRUELEN
        INTEGER SYSTEMFUNCTION
C
        INTEGER L,L1,L2
        INTEGER ISYSTEM
        CHARACTER*255 CADENA
        CHARACTER*1 COPC
        LOGICAL LRUN,LMANUAL,LHTML
C------------------------------------------------------------------------------
        CALL LRUNX(LRUN,LMANUAL,LHTML)
C------------------------------------------------------------------------------
        IF(CADENA(1:10).EQ.'endprogram')THEN
          IF(LRUN)THEN
            WRITE(*,101) 'endprogram'
          END IF
          IF(LMANUAL)THEN
            WRITE(*,101) '\\ttshade{endprogram}'
          END IF
          IF(LHTML)THEN
            WRITE(*,101) '<FONT COLOR="#FF0000">endprogram</FONT>'
          END IF
          WRITE(*,100) 'Do you want to STOP the program (y/n) [n] ? '
          READ(*,'(A)') COPC
          IF(LRUN)THEN
            WRITE(*,101)COPC
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
C------------------------------------------------------------------------------
        ELSEIF(CADENA(1:10).EQ.'newpagenew')THEN
          IF(LRUN)THEN
            WRITE(*,101) 'newpagenew'
          END IF
          IF(LMANUAL)THEN
            WRITE(*,101) '\\ttshade{newpagenew}'
          END IF
          IF(LHTML)THEN
            WRITE(*,101) '<FONT COLOR="#FF0000">newpagenew</FONT>'
          END IF
          WRITE(*,*)
          WRITE(*,100) 'Are you inserting a newpage call (y/n) [n] ? '
          READ(*,'(A)') COPC
          IF(LRUN)THEN
            WRITE(*,101)COPC
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
C------------------------------------------------------------------------------
C Capturamos en un fichero GIF la imagen actual en la ventana de PGPLOT.
C La variable cadena debe ser algo como "capturegif#filegif#something_else",
C donde filegif es el nombre del fichero GIF, y something_else es lo que
C finalmente va a retornar esta función a través de la variable CADENA (para
C que la función readi, readf, readc, o la que sea, continue trabajando como
C si no hubiera pasado nada. Si, además, existe el fichero .running_HLPHTML,
C se muestra en pantalla el codigo HTML necesario para insertar la figura
C GIF correspondiente.
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
          ISYSTEM=SYSTEMFUNCTION('xwd -out '//CADENA(L1+1:L2-1)//'.tmp '
     +     //'-name "PGPLOT Window 1"\0')
          ISYSTEM=SYSTEMFUNCTION('convert '//CADENA(L1+1:L2-1)//'.tmp '
     +     //CADENA(L1+1:L2-1)//'.gif\0')
          ISYSTEM=SYSTEMFUNCTION('\\rm '//CADENA(L1+1:L2-1)//'.tmp\0')
          !mostramos en pantalla el código HTML para insertar la figura GIF
          IF(LHTML)THEN
            WRITE(*,101) '<center>'
            WRITE(*,101) '<img SRC="'//CADENA(L1+1:L2-1)//
     +       '.gif" BORDER="2">'
            WRITE(*,101) '</center>'
          END IF
          CADENA=CADENA(L2+1:)               !devolvemos solo lo que hace falta
C------------------------------------------------------------------------------
        END IF
C------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
        END
