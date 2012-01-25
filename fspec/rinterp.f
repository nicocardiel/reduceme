C------------------------------------------------------------------------------
C Version 13-October-2007                                       File: rinterp.f
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
C REAL FUNCTION RINTERP(LAMBDA)
C
C Input: LAMBDA
C Output: RINTERP (function)
C
C Calculate the extinction value A(lambda)/E(B-V) for a given wavelength
C LAMBDA. The tabulated data correspond to Savage & Mathis (1979, Ann.
C Rev. Astron. Astrophys., 17, 13).
C
C REAL LAMBDA -> input wavelength
C
Comment
C------------------------------------------------------------------------------
        REAL FUNCTION RINTERP(LAMBDA)
        IMPLICIT NONE
        INTEGER NPMAX
        PARAMETER(NPMAX=100)
        REAL LAMBDA
        INTEGER I
        REAL X(NPMAX),Y(NPMAX)
        INTEGER NPTOS
        INTEGER N1,N2
C------------------------------------------------------------------------------
        NPTOS=13
        DATA(X(I),Y(I),I=1,13)/2300.000,8.87,
     +                         2400.000,8.00,
     +                         2500.000,7.29,
     +                         2740.000,6.20,
     +                         3440.000,4.90,
     +                         4000.000,4.40,
     +                         4400.000,4.10,
     +                         5500.000,3.10,
     +                         7000.000,2.32,
     +                         9000.000,1.50,
     +                        12500.000,0.87,
     +                        22000.000,0.38,
     +                        34000.000,0.16/
C
        IF((X(1).GT.LAMBDA).OR.(X(13).LT.LAMBDA))THEN
          WRITE(*,101)'ERROR: wavelength out of range in function '//
     +     'RINTERP.'
          WRITE(*,101)'Internal extinction set to 0.0'
          WRITE(*,100)'Press <CR> to continue...'
          READ(*,*)
          RINTERP=0.0
          RETURN
        END IF
C
        N1=0
        N2=0
        DO I=1,NPTOS
          IF(X(I).LE.LAMBDA)N1=I
        END DO
        DO I=NPTOS,1,-1
          IF(X(I).GE.LAMBDA)N2=I
        END DO
        IF((N1.EQ.0).OR.(N2.EQ.0))THEN
          STOP 'FATAL ERROR: N1.EQ.0 or N2.EQ.0 in RINTERP'
        END IF
        IF(N1.EQ.N2)THEN
          RINTERP=Y(N1)
        ELSE
          RINTERP=Y(N1)+((LAMBDA-X(N1))/(X(N2)-X(N1)))*(Y(N2)-Y(N1))
        END IF
C
100     FORMAT(A,$)
101     FORMAT(A)
        END
