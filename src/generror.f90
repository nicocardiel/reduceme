!------------------------------------------------------------------------------
! Version 20-November-2026                                   file: generror.f90
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
! Program: generror
! Classification: error handling
! Description: Creates an error file from an initial image, taking into
! account the readout noise and gain of the detector.
!
!Comment
! 
! Genera un fichero de errores a partir de una imagen inicial y 
! teniendo en cuenta el ruido de lectura y la ganancia del detector
!
        PROGRAM GENERROR
        USE REDLIB_M
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        REAL READF
!
        INTEGER I,J
        REAL ERR(NCMAX,NSMAX)
        REAL RNOISE,GAIN
        REAL RNOISE2
        CHARACTER*255 INFILE,OUTFILE
        LOGICAL LNEGATIVE
!------------------------------------------------------------------------------
        THISPROGRAM='generror'
        CALL WELCOME('20-November-2026')
! leemos la imagen inicial
        WRITE(*,100,ADVANCE='NO')'Input file name'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        DO I=1,NSCAN
          READ(20) (ERR(J,I),J=1,NCHAN)
        END DO
        CLOSE(20)
! entramos el ruido de lectura y la ganancia
        WRITE(*,100,ADVANCE='NO')&
     &     'Read-out Noise (counts=Digital Units, DU)'
        RNOISE=READF('@')
        RNOISE2=RNOISE*RNOISE
        WRITE(*,100,ADVANCE='NO')&
     &     'Gain (elect/DU)                          '
        GAIN=READF('@')
! calculamos el error en cada pixel
        LNEGATIVE=.FALSE.              !la imagen no contiene valores negativos
        DO I=1,NSCAN
          DO J=1,NCHAN
            IF(ERR(J,I).GE.0.)THEN
              ERR(J,I)=SQRT(ERR(J,I)/GAIN+RNOISE2)
            ELSE
              LNEGATIVE=.TRUE.
              ERR(J,I)=SQRT(-ERR(J,I)/GAIN+RNOISE2)
            END IF
          END DO
        END DO
        IF(LNEGATIVE)THEN
          WRITE(*,101)'WARNING; the image contains negative '//&
     &     'pixel value(s)'
        END IF
! salvamos el fichero de errores
        WRITE(*,100,ADVANCE='NO')'Output file name '
        CALL GUESSEF(INFILE,OUTFILE)
        OUTFILE=OUTFILEX(30,OUTFILE,NSCAN,NCHAN,STWV,DISP,1,.TRUE.)
        DO I=1,NSCAN
          WRITE(30) (ERR(J,I),J=1,NCHAN)
        END DO
        CLOSE(30)
!
        STOP
100     FORMAT(A)
101     FORMAT(A)
        END
