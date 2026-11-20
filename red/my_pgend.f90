!------------------------------------------------------------------------------
! Version 20-November-2026                                   File: my_pgend.f90
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
! SUBROUTINE MY_PGEND
!
! If the file .running_HLPHTML exist, this subroutine allows to capture the 
! last XServe image before calling PGEND.
!
!Comment
!------------------------------------------------------------------------------
        SUBROUTINE MY_PGEND
        IMPLICIT NONE
!
!
!
        INTEGER L
        INTEGER ISYSTEM
        CHARACTER*255 OUTGIF
        LOGICAL LOK
!------------------------------------------------------------------------------
        INQUIRE(FILE='.running_HLPHTML',EXIST=LOK)
        IF(LOK)THEN
          READ(*,101) OUTGIF
          L=TRUELEN(OUTGIF)
          IF(L.GT.0)THEN
            ISYSTEM=SYSTEMFUNCTION('xwd -out '//OUTGIF(1:L)//&
     &       '.tmp -name "PGPLOT Window 1"\0')
            ISYSTEM=SYSTEMFUNCTION('convert '//OUTGIF(1:L)//'.tmp '//&
     &       OUTGIF(1:L)//'.gif\0')
            ISYSTEM=SYSTEMFUNCTION('\\rm '//OUTGIF(1:L)//'.tmp \0')
            WRITE(*,101) '<center>'
            WRITE(*,101) '<img SRC="'//OUTGIF(1:L)//'.gif" BORDER="2">'
            WRITE(*,101) '</center>'
          END IF
        END IF
        CALL PGEND
!cc100     FORMAT(A,$)
101     FORMAT(A)
        END
