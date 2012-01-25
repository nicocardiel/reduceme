C------------------------------------------------------------------------------
C Version 20-March-1997                                          file:rvrelat.f
C------------------------------------------------------------------------------
C Copyright N. Cardiel & J. Gorgas, Departamento de Astrofisica
C Universidad Complutense de Madrid, 28040-Madrid, Spain
C E-mail: ncl@astrax.fis.ucm.es or fjg@astrax.fis.ucm.es
C------------------------------------------------------------------------------
C This program is free software; you can redistribute it and/or modify it
C under the terms of the GNU General Public License as published by the Free
C Software Foundation; either version 2 of the License, or (at your option) any
C later version. See the file gnu-public-license.txt for details.
C------------------------------------------------------------------------------
Comment
C
C Program: rvrelat
C Classification: miscellany
C Description: Transforms radial velocities computed from classical formulae
C into radial velocities computed with relativistic expressions.
C
Comment
C
        PROGRAM RVRELAT
        IMPLICIT NONE
        REAL READF
C
        REAL C
        PARAMETER (C=299792.46)
C
        REAL RV1,RV2
        REAL K,Z
C------------------------------------------------------------------------------
        WRITE(*,100)'Classical Rv (km/sec)'
        RV1=READF('@')
        IF(RV1.LT.0.0) STOP 'ERROR: use only positive velocities.'
        K=1.+RV1/C
        RV2=(K*K-1.)/(K*K+1.)*C
! estas formulas son equivalentes
!       K=RV1/C
!       RV2=C*((1+K)*(1+K)-1)/((1+K)*(1+K)+1)
        WRITE(*,100)'Relativistic Rv....: '
        WRITE(*,*) RV2
        WRITE(*,100)'DOUBLE CHECK--> (1+z)='
        Z=RV2/C
        WRITE(*,*) K,(1.+Z)/SQRT(1.-Z*Z)
        STOP
100     FORMAT(A,$)
        END
