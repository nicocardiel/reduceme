!------------------------------------------------------------------------------
! Version 20-November-2026                                    File: buttspr.f90
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
! SUBROUTINE BUTTSPR(X1,X2,Y1,Y2)
!
! Input: X1,X2,Y1,Y2
! 
! Set the plot region limits.
!
! REAL X1 -> x-coordinate of the left hand edge of the plot region viewport,
!      in normalized device coordinates
! REAL X2 -> x-coordinate of the right hand edge of the plot region viewport,
!      in normalized device coordinates
! REAL Y1 -> y-coordinate of the bottom edge of the plot region viewport,
!      in normalized device coordinates
! REAL Y2 -> y-coordinate of the top edge of the plot region viewport,
!      in normalized device coordinates
!
!Comment
!------------------------------------------------------------------------------
        SUBROUTINE BUTTSPR(X1,X2,Y1,Y2)
        IMPLICIT NONE
        REAL X1,X2,Y1,Y2
!------------------------------------------------------------------------------
        X1VPORT=X1
        X2VPORT=X2
        Y1VPORT=Y1
        Y2VPORT=Y2
        END
