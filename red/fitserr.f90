!------------------------------------------------------------------------------
! Version 24-November-2026                                    file: fitserr.f90
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
! SUBROUTINE FITSERR(ISTATUS,CMSG)
!
! Input: ISTATUS,CMSG
!
! Prints the error message CMSG followed by the CFITSIO description of the
! status ISTATUS and the messages of the CFITSIO error stack (which is
! emptied).
!
! INTEGER       ISTATUS -> CFITSIO status (greater than 0)
! CHARACTER*(*) CMSG -> error message
!
!Comment
!------------------------------------------------------------------------------
        SUBROUTINE FITSERR(ISTATUS,CMSG)
        IMPLICIT NONE
        INTEGER ISTATUS
        CHARACTER*(*) CMSG
!
        CHARACTER*30 CERRTEXT
        CHARACTER*80 CERRMSG
!------------------------------------------------------------------------------
        WRITE(*,101)'ERROR: '//CMSG(1:TRUELEN(CMSG))
        CALL FTGERR(ISTATUS,CERRTEXT)
        WRITE(*,'(A,I4,A)')'FITSIO status ',ISTATUS,': '//&
     &   CERRTEXT(1:TRUELEN(CERRTEXT))
! error stack (FTGMSG returns a blank message when the stack is empty)
        CALL FTGMSG(CERRMSG)
        DO WHILE(TRUELEN(CERRMSG).GT.0)
          WRITE(*,101)'  '//CERRMSG(1:TRUELEN(CERRMSG))
          CALL FTGMSG(CERRMSG)
        END DO
101     FORMAT(A)
        END
