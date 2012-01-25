C------------------------------------------------------------------------------
C                                                                    file:rvz.f
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
C Program: rvz
C Classification: wavelengths
C Description: Computes the corresponding radial velocity for a given redshift.
C
Comment
C
        PROGRAM RVZ
        IMPLICIT NONE
C
        REAL C
        PARAMETER (C=299792.46)
C
        REAL Z,X
C
        WRITE(*,100) 'Redshift? '
        READ(*,*) Z
C
        X=SQRT(1.-(1.+(1.+Z)*(1.+Z))*(1.-(1.+Z)*(1.+Z)))-1.
        X=X/(1.+(1.+Z)*(1.+Z))
        WRITE(*,100) 'Radial velocity: '
        WRITE(*,*) X*C
C
        STOP
100     FORMAT(A,$)
        END
