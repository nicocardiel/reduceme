!------------------------------------------------------------------------------
! Version 20-November-2026                                     file:rvrelat.f90
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
! Program: rvrelat
! Classification: miscellany
! Description: Transforms radial velocities computed from classical formulae
! into radial velocities computed with relativistic expressions.
!
!Comment
!
        PROGRAM RVRELAT
        IMPLICIT NONE
        REAL READF
!
        REAL C
        PARAMETER (C=299792.46)
!
        REAL RV1,RV2
        REAL K,Z
!------------------------------------------------------------------------------
        WRITE(*,100,ADVANCE='NO')'Classical Rv (km/sec)'
        RV1=READF('@')
        IF(RV1.LT.0.0) STOP 'ERROR: use only positive velocities.'
        K=1.+RV1/C
        RV2=(K*K-1.)/(K*K+1.)*C
! estas formulas son equivalentes
!       K=RV1/C
!       RV2=C*((1+K)*(1+K)-1)/((1+K)*(1+K)+1)
        WRITE(*,100,ADVANCE='NO')'Relativistic Rv....: '
        WRITE(*,*) RV2
        WRITE(*,100,ADVANCE='NO')'DOUBLE CHECK--> (1+z)='
        Z=RV2/C
        WRITE(*,*) K,(1.+Z)/SQRT(1.-Z*Z)
        STOP
100     FORMAT(A)
        END
