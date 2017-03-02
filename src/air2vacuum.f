C------------------------------------------------------------------------------
C                                                             file:air2vacuum.f
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
C Program: air2vacuum
C Classification: wavelengths
C Description: Transforms wavelengths from air to vacuum using the equations
C from Greisen et al. 2006 (Representation of spectral coordinates in FITS,
C A&A, 446, 747-771).
C
Comment
C
        PROGRAM AIR2VACUUM
        IMPLICIT NONE
C
        DOUBLE PRECISION GREISEN2006
        DOUBLE PRECISION CIDDOR1996
C
        INTEGER I
        DOUBLE PRECISION L_AIR,L_VACUUM
        DOUBLE PRECISION WAVE_AIR
        DOUBLE PRECISION N
C------------------------------------------------------------------------------
        WRITE(*,100) 'Wavelength in air or in vacuum (in Angs.)? '
        READ(*,*) L_AIR
C------------------------------------------------------------------------------
C Fórmula antigua, creo que de Cox (2000)
!       WAVE_AIR=10000.0D0/L_AIR !in cm^-1
!       N=1.0D0+6432.8D-8+2949810.0D0/(146.0D8-WAVE_AIR*WAVE_AIR)+
!    +   25540.0D0/(41D8-WAVE_AIR*WAVE_AIR)
!       L_VACUUM=L_AIR*N
!       WRITE(*,100) 'Index of refraction (dry air at standard'//
!    +   ' temperature and pressure): '
!       WRITE(*,*) N
!       WRITE(*,100) 'Vacuum wavelength, air-vacuum (Angs.): '
!       WRITE(*,*) L_VACUUM,L_AIR-L_VACUUM
!       WRITE(*,100) 'Re-predicted air wavelength: '
!       WRITE(*,*) L_AIR/N
C------------------------------------------------------------------------------
        WAVE_AIR=L_AIR/10000.0 !en micras
        DO I=1,2
          IF(I.EQ.1)THEN
            WRITE(*,101) '* Using Greisen et al. (2006):'
            N=GREISEN2006(WAVE_AIR)
          ELSE
            WRITE(*,101) '* Using Ciddor (1996):'
            N=CIDDOR1996(WAVE_AIR)
          END IF
C
          L_VACUUM=L_AIR*N
          WRITE(*,100) 'Refractive index (dry air at standard'//
     +     ' temperature and pressure): '
          WRITE(*,*) N
          WRITE(*,100) '==> Air, Vacuum, Air-Vacuum (Angs.): '
          WRITE(*,'(F10.4,2X,F10.4,2X,F8.4)')
     +     L_AIR,L_VACUUM,L_AIR-L_VACUUM
          WRITE(*,100) '==> Vacuum, Air, Vacuum-Air (Angs.): '
          WRITE(*,'(F10.4,2X,F10.4,2X,F8.4)')
     +     L_AIR,L_AIR/N,L_AIR-L_AIR/N
        END DO
C------------------------------------------------------------------------------
        STOP
100     FORMAT(A,$)
101     FORMAT(A)
        END
C
C******************************************************************************
C Compute refractive index using Eq. 65 in Greisen et al. 2006
C (Representation of spectral coordinates in FITS, A&A, 446, 747-771).
        DOUBLE PRECISION FUNCTION GREISEN2006(WAVE_AIR)
        IMPLICIT NONE
        DOUBLE PRECISION WAVE_AIR  !wavelength in microns
C
        DOUBLE PRECISION N
C
        N=287.6155D0+
     +   1.62887D0/(WAVE_AIR*WAVE_AIR)+
     +   0.01360D0/(WAVE_AIR*WAVE_AIR*WAVE_AIR*WAVE_AIR)
        N=1.D0+N*1.D-6
        GREISEN2006=N
        END
C
C******************************************************************************
C Compute refractive index using Eq. 1 in Prieto 2011
C (see http://www.as.utexas.edu/~hebe/apogee/docs/air_vacuum.pdf)
C which provides the expresion given by Ciddor 1996.
        DOUBLE PRECISION FUNCTION CIDDOR1996(WAVE_AIR)
        IMPLICIT NONE
        DOUBLE PRECISION WAVE_AIR  !wavelength in microns
C
        DOUBLE PRECISION N
        DOUBLE PRECISION A,B1,B2,C1,C2
C
        A=0.D0
        B1=5.792105D-2
        B2=1.67917D-3
        C1=238.0185D0
        C2=57.362D0
        N=A+B1/(C1-1.D0/(WAVE_AIR**2))+B2/(C2-1.D0/(WAVE_AIR**2))
        N=1.D0+N
        CIDDOR1996=N
        END
