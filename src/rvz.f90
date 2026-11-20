!------------------------------------------------------------------------------
!                                                                    file:rvz.f
!------------------------------------------------------------------------------
! Copyright N. Cardiel & J. Gorgas, Departamento de Astrofisica
! Universidad Complutense de Madrid, 28040-Madrid, Spain
! E-mail: ncl@astrax.fis.ucm.es or fjg@astrax.fis.ucm.es
!------------------------------------------------------------------------------
! This program is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the Free
! Software Foundation; either version 2 of the License, or (at your option) any
! later version. See the file gnu-public-license.txt for details.
!------------------------------------------------------------------------------
!Comment
!
! Program: rvz
! Classification: wavelengths
! Description: Computes the corresponding radial velocity for a given redshift.
!
!Comment
!
        PROGRAM RVZ
        IMPLICIT NONE
!
        REAL C
        PARAMETER (C=299792.46)
!
        REAL Z,X
!
        WRITE(*,100,ADVANCE='NO') 'Redshift? '
        READ(*,*) Z
!
        X=SQRT(1.-(1.+(1.+Z)*(1.+Z))*(1.-(1.+Z)*(1.+Z)))-1.
        X=X/(1.+(1.+Z)*(1.+Z))
        WRITE(*,100,ADVANCE='NO') 'Radial velocity: '
        WRITE(*,*) X*C
!
        STOP
100     FORMAT(A)
        END
