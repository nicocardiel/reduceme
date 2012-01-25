C------------------------------------------------------------------------------
C Version 6-December-1996                                      file: generror.f
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
C Program: generror
C Classification: error handling
C Description: Creates an error file from an initial image, taking into
C account the readout noise and gain of the detector.
C
Comment
C 
C Genera un fichero de errores a partir de una imagen inicial y 
C teniendo en cuenta el ruido de lectura y la ganancia del detector
C
        PROGRAM GENERROR
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        REAL READF
C
        INTEGER I,J
        REAL ERR(NCMAX,NSMAX)
        REAL RNOISE,GAIN
        REAL RNOISE2
        CHARACTER*75 INFILE,OUTFILE
        LOGICAL LNEGATIVE
C------------------------------------------------------------------------------
        THISPROGRAM='generror'
        CALL WELCOME('6-December-1996')
C leemos la imagen inicial
        WRITE(*,100)'Input file name'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        DO I=1,NSCAN
          READ(20) (ERR(J,I),J=1,NCHAN)
        END DO
        CLOSE(20)
C entramos el ruido de lectura y la ganancia
        WRITE(*,100)'Read-out Noise (counts=Digital Units, DU)'
        RNOISE=READF('@')
        RNOISE2=RNOISE*RNOISE
        WRITE(*,100)'Gain (elect/DU)                          '
        GAIN=READF('@')
C calculamos el error en cada pixel
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
          WRITE(*,101)'WARNING; the image contains negative '//
     +     'pixel value(s)'
        END IF
C salvamos el fichero de errores
        WRITE(*,100)'Output file name '
        CALL GUESSEF(INFILE,OUTFILE)
        OUTFILE=OUTFILEX(30,OUTFILE,NSCAN,NCHAN,STWV,DISP,1,.TRUE.)
        DO I=1,NSCAN
          WRITE(30) (ERR(J,I),J=1,NCHAN)
        END DO
        CLOSE(30)
C
        STOP
100     FORMAT(A,$)
101     FORMAT(A)
        END
