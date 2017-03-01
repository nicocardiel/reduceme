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
C Fórmulas nuevas (Eq. 65 in Greisen et al. 2006)
        WAVE_AIR=L_AIR/10000.0 !en micras
        N=287.6155D0+
     +   1.62887D0/(WAVE_AIR*WAVE_AIR)+
     +   0.01360D0/(WAVE_AIR*WAVE_AIR*WAVE_AIR*WAVE_AIR)
        N=1.D0+N*1.D-6
        L_VACUUM=L_AIR*N
        WRITE(*,100) 'Refraction index (dry air at standard'//
     +   ' temperature and pressure): '
        WRITE(*,*) N
        WRITE(*,100) '==> Air, Vacuum, Air-Vacuum (Angs.): '
        WRITE(*,'(F10.4,2X,F10.4,2X,F8.4)')L_AIR,L_VACUUM,L_AIR-L_VACUUM
        WRITE(*,100) '==> Vacuum, Air, Vacuum-Air (Angs.): '
        WRITE(*,'(F10.4,2X,F10.4,2X,F8.4)')L_AIR,L_AIR/N,L_AIR-L_AIR/N
C------------------------------------------------------------------------------
        STOP
100     FORMAT(A,$)
        END
