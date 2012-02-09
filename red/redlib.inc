C------------------------------------------------------------------------------
C Version 07-September-2007                                    File: redlib.inc
C------------------------------------------------------------------------------
C Copyright N. Cardiel & J. Gorgas, Departamento de Astrofisica
C Universidad Complutense de Madrid, 28040-Madrid, Spain
C E-mail: ncl@astrax.fis.ucm.es or fjg@astrax.fis.ucm.es
C------------------------------------------------------------------------------
C This file is free software; you can redistribute it and/or modify it
C under the terms of the GNU General Public License as published by the Free
C Software Foundation; either version 2 of the License, or (at your option) any
C later version. See the file gnu-public-license.txt for details.
C------------------------------------------------------------------------------
        INTEGER NSMAX                                            !maximum NSCAN
        PARAMETER (NSMAX=2700)
C
        INTEGER NCMAX                                            !maximum NCHAN
        PARAMETER (NCMAX=4400)
C
        INTEGER NMAXFFT
        PARAMETER (NMAXFFT=2**16)     !2**16 maximum spectrum size for FFT work
C------------------------------------------------------------------------------
        INTEGER NMAXBUFF                   !maximum no. of Buffers in plotsplus
        PARAMETER (NMAXBUFF=6)
C
        CHARACTER*12 CLAVE_RED         !fixed identification for REDUCEME files
        PARAMETER (CLAVE_RED='abcdefghijkl')
C
        INTEGER NBDMAX        !maximum no. of bands allowed for a generic index
        PARAMETER (NBDMAX=198)                       !note that NBDMAX=NWVMAX/2
C
        INTEGER NWVMAX      !2*maximum no. of bands allowed for a generic index
        PARAMETER (NWVMAX=396)                       !see subroutine selindex.f
C
        INTEGER NINDMAX                         !maximum no. of indices defined
        PARAMETER (NINDMAX=100)                      !see subroutine selindex.f
C
        INTEGER NLINMAX                           !maximum no. of lines defined
        PARAMETER (NLINMAX=40)                       !see subroutine sellines.f
C
        INTEGER NPBANDMAX           !maximum no. of points to define a bandpass
        PARAMETER (NPBANDMAX=40)                      !see subroutine selband.f
C
        INTEGER MAX_ID_RED      !maximum number of simultaneous graphic devices
        PARAMETER (MAX_ID_RED=8)
C
        CHARACTER*12 CREDUCEVERSION                            !current version
        PARAMETER (CREDUCEVERSION='REDUCEMEv4.1')
C------------------------------------------------------------------------------
        INTEGER NSCAN                 !image dimension in the spatial direction
        INTEGER NCHAN              !image dimension in the wavelength direction
        REAL    STWV                     !central wavelength of the first pixel
        REAL    DISP       !dispersion (Angs/pixel) in the wavelength direction
C------------------------------------------------------------------------------
        REAL AIRMASS                                                   !airmass
        REAL TIMEXPOS                                            !exposure time
        CHARACTER*255 OBJECT                       !name of the observed object
        CHARACTER*255 FITSFILE             !file name of the original FITS file
        CHARACTER*255 COMMENT       !comment to be included in the image header
        CHARACTER*20 THISPROGRAM                !current program being executed
C NOTE: if OBJECT, FITSFILE or COMMENT are redimensioned, the changes must
C also be performed in basicred.f and imath.f
C------------------------------------------------------------------------------
        COMMON/BLKRED01/AIRMASS
        COMMON/BLKRED02/TIMEXPOS
        COMMON/BLKRED03/OBJECT
        COMMON/BLKRED04/FITSFILE
        COMMON/BLKRED05/COMMENT
        COMMON/BLKRED06/THISPROGRAM
C------------------------------------------------------------------------------
!100     FORMAT(A,$)
!101     FORMAT(A)
!110     FORMAT(A,I6)
C------------------------------------------------------------------------------
