C------------------------------------------------------------------------------
C Version 6-June-2003                                          file: resample.f
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
C Program: resample
C Classification: arithmetic & manipulations
C Description: Transforms an image with an initial STWV and DISP into another 
C image with a different STWV and DISP.
C
Comment
C
C Este programa transforma una imagen con una determinada STWV y DISP a otra
C imagen con una nueva STWV y DISP.
C
        PROGRAM RESAMPLE
        IMPLICIT NONE
C
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER READILIM
        REAL READF
C
        INTEGER I,J
        INTEGER NCHAN2
        REAL A(NCMAX,NSMAX),ERR(NCMAX,NSMAX)
        REAL STWV2,DISP2
        REAL S(NCMAX),SS(NCMAX),X(NCMAX),XX(NCMAX)
        REAL FACTOR
        CHARACTER*1 CERR,CFLUX
        CHARACTER*50 CDUMMY
        CHARACTER*75 INFILE,ERRFILE,OUTFILE
C------------------------------------------------------------------------------
        THISPROGRAM='resample'
        CALL WELCOME('6-June-2003')
C imagen de entrada (y errores)
        WRITE(*,100)'Work with error images (y/n) '
        CERR(1:1)=READC('n','yn')
        WRITE(*,100)'Preserve flux (y/n) '
        CFLUX(1:1)=READC('n','yn')
C
        WRITE(*,100)'Input file name'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        DO I=1,NSCAN
          READ(20) (A(J,I),J=1,NCHAN)
        END DO
        CLOSE(20)
        IF(CERR.EQ.'y')THEN
          CALL GUESSEF(INFILE,ERRFILE)
          WRITE(*,100)'Error file name '
          ERRFILE=INFILEX(21,ERRFILE,NSCAN,NCHAN,STWV,DISP,21,.TRUE.) !match
          DO I=1,NSCAN
            READ(21) (ERR(J,I),J=1,NCHAN)
          END DO
          CLOSE(21)
        END IF
C------------------------------------------------------------------------------
C nuevos parametros de la imagen de salida
        WRITE(CDUMMY,*) STWV
        WRITE(*,100)'New STWV '
        STWV2=READF(CDUMMY)
C
        WRITE(CDUMMY,*) DISP
        WRITE(*,100)'New DISP '
        DISP2=READF(CDUMMY)
C
        WRITE(CDUMMY,*) NCHAN
        WRITE(*,100)'New NCHAN '
        NCHAN2=READILIM(CDUMMY,1,NCMAX)
C------------------------------------------------------------------------------
        IF(CFLUX.EQ.'y')THEN
          FACTOR=DISP2/DISP
        ELSE
          FACTOR=1.0
        END IF
C------------------------------------------------------------------------------
C creamos los nuevos espectros
        DO J=1,NCHAN
          X(J)=STWV+REAL(J-1)*DISP
        END DO
C
        DO I=1,NSCAN
          DO J=1,NCHAN
            S(J)=A(J,I)
          END DO
          CALL REBINING(X,S,NCHAN,XX,SS,NCHAN2,STWV2,DISP2)
          DO J=1,NCHAN2
            A(J,I)=SS(J)*FACTOR
          END DO
        END DO
C
        IF(CERR.EQ.'y')THEN
          DO I=1,NSCAN
            DO J=1,NCHAN
              S(J)=ERR(J,I)
            END DO
            CALL REBINING(X,S,NCHAN,XX,SS,NCHAN2,STWV2,DISP2)
            DO J=1,NCHAN2
              ERR(J,I)=SS(J)*FACTOR
            END DO
          END DO
        END IF
C------------------------------------------------------------------------------
C salvamos la nueva imagen
        WRITE(*,100)'Output file name'
        OUTFILE=OUTFILEX(30,'@',NSCAN,NCHAN2,STWV2,DISP2,1,.FALSE.)
        DO I=1,NSCAN
          WRITE(30) (A(J,I),J=1,NCHAN2)
        END DO
        CLOSE(30)
        IF(CERR.EQ.'y')THEN
          CALL GUESSEF(OUTFILE,ERRFILE)
          WRITE(*,100)'Output error file name '
          OUTFILE=OUTFILEX(31,ERRFILE,NSCAN,NCHAN2,STWV2,DISP2,1,.TRUE.)
          DO I=1,NSCAN
            WRITE(31) (ERR(J,I),J=1,NCHAN2)
          END DO
          CLOSE(31)
        END IF
C------------------------------------------------------------------------------
        STOP
100     FORMAT(A,$)
        END
