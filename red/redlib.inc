!------------------------------------------------------------------------------
! Version 20-November-2026                                     File: redlib.inc
!------------------------------------------------------------------------------
! Copyright N. Cardiel & J. Gorgas, Departamento de Astrofisica
! Universidad Complutense de Madrid, 28040-Madrid, Spain
! E-mail: ncl@astrax.fis.ucm.es or fjg@astrax.fis.ucm.es
!------------------------------------------------------------------------------
! This file is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the Free
! Software Foundation; either version 2 of the License, or (at your option) any
! later version. See the file gnu-public-license.txt for details.
!------------------------------------------------------------------------------
! Image dimensions and wavelength calibration: local variables of every
! program unit that includes this file. The global parameters and image
! header variables are in the module REDLIB_M (futils/redlib_m.f90.in).
!------------------------------------------------------------------------------
        INTEGER NSCAN                 !image dimension in the spatial direction
        INTEGER NCHAN              !image dimension in the wavelength direction
        REAL    STWV                     !central wavelength of the first pixel
        REAL    DISP       !dispersion (Angs/pixel) in the wavelength direction
!------------------------------------------------------------------------------
