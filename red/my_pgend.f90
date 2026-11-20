C------------------------------------------------------------------------------
C Version 07-September-2007                                    File: my_pgend.f
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
C SUBROUTINE MY_PGEND
C
C If the file .running_HLPHTML exist, this subroutine allows to capture the 
C last XServe image before calling PGEND.
C
Comment
C------------------------------------------------------------------------------
        SUBROUTINE MY_PGEND
        IMPLICIT NONE
C
        INTEGER TRUELEN
C
        INTEGER SYSTEMFUNCTION
C
        INTEGER L
        INTEGER ISYSTEM
        CHARACTER*255 OUTGIF
        LOGICAL LOK
C------------------------------------------------------------------------------
        INQUIRE(FILE='.running_HLPHTML',EXIST=LOK)
        IF(LOK)THEN
          READ(*,101) OUTGIF
          L=TRUELEN(OUTGIF)
          IF(L.GT.0)THEN
            ISYSTEM=SYSTEMFUNCTION('xwd -out '//OUTGIF(1:L)//
     +       '.tmp -name "PGPLOT Window 1"\0')
            ISYSTEM=SYSTEMFUNCTION('convert '//OUTGIF(1:L)//'.tmp '//
     +       OUTGIF(1:L)//'.gif\0')
            ISYSTEM=SYSTEMFUNCTION('\\rm '//OUTGIF(1:L)//'.tmp \0')
            WRITE(*,101) '<center>'
            WRITE(*,101) '<img SRC="'//OUTGIF(1:L)//'.gif" BORDER="2">'
            WRITE(*,101) '</center>'
          END IF
        END IF
        CALL PGEND
ccc100     FORMAT(A,$)
101     FORMAT(A)
        END
