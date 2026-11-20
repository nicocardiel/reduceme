C------------------------------------------------------------------------------
C Version 21-August-2004                                          File: lrunx.f
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
C SUBROUTINE LRUNX(LRUN,LMANUAL,LHTML)
C
C Output: LRUN,LMANUAL,LHTML
C
C Determine whether files .running_RUN, .running_MANUAL and .running_HLPHTML
C exist in current the directory.
C
C LOGICAL LRUN -> .TRUE. if file .running_RUN exists (.FALSE. otherwise)
C LOGICAL LMANUAL -> .TRUE. if file .running_MANUAL exists (.FALSE. otherwise)
C LOGICAL LHTML -> .TRUE. if file .running_HLPHTML exists (.FALSE. otherwise)
C
Comment
C------------------------------------------------------------------------------
        SUBROUTINE LRUNX(LRUN,LMANUAL,LHTML)
        IMPLICIT NONE
        LOGICAL LRUN,LMANUAL,LHTML
C
        INQUIRE(FILE='.running_RUN',EXIST=LRUN)
        INQUIRE(FILE='.running_MANUAL',EXIST=LMANUAL)
        INQUIRE(FILE='.running_HLPHTML',EXIST=LHTML)
        END
