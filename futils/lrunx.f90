!------------------------------------------------------------------------------
! Version 21-August-2004                                          File: lrunx.f
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
! SUBROUTINE LRUNX(LRUN,LMANUAL,LHTML)
!
! Output: LRUN,LMANUAL,LHTML
!
! Determine whether files .running_RUN, .running_MANUAL and .running_HLPHTML
! exist in current the directory.
!
! LOGICAL LRUN -> .TRUE. if file .running_RUN exists (.FALSE. otherwise)
! LOGICAL LMANUAL -> .TRUE. if file .running_MANUAL exists (.FALSE. otherwise)
! LOGICAL LHTML -> .TRUE. if file .running_HLPHTML exists (.FALSE. otherwise)
!
!Comment
!------------------------------------------------------------------------------
        SUBROUTINE LRUNX(LRUN,LMANUAL,LHTML)
        IMPLICIT NONE
        LOGICAL LRUN,LMANUAL,LHTML
!
        INQUIRE(FILE='.running_RUN',EXIST=LRUN)
        INQUIRE(FILE='.running_MANUAL',EXIST=LMANUAL)
        INQUIRE(FILE='.running_HLPHTML',EXIST=LHTML)
        END
