C------------------------------------------------------------------------------
C Version 10-November-1997                                      file: rvshift.f
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
C Program: rvshift
C Classification: wavelengths
C Description: Applies a radial velocity shift to an image.
C
Comment
C
C Traslada un espectro una cierta velocidad, manteniendolo en escala lineal
C con la misma STWV y DISP que el espectro de entrada.
C
        PROGRAM RVSHIFT
C
        IMPLICIT NONE
C
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        REAL READF
C
        INTEGER I,J
        REAL S(NCMAX),SS(NCMAX)
        REAL RADVEL
        CHARACTER*1 CERR
        CHARACTER*75 INFILE,ERRINFILE,OUTFILE,ERROUTFILE
C------------------------------------------------------------------------------
        THISPROGRAM='rvshift'
        CALL WELCOME('10-November-1997')
C------------------------------------------------------------------------------
        WRITE(*,100)'Work with error images (y/n) '
        CERR(1:1)=READC('y','yn')
C
        WRITE(*,100)'Input file name'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        IF((STWV.EQ.0.0).AND.(DISP.EQ.0.0))THEN
          WRITE(*,101)'FATAL ERROR: STWV=0.0 and DISP=0.0!'
          CLOSE(20)
          STOP
        END IF
C
        IF(CERR.EQ.'y')THEN
          WRITE(*,100)'Input error file name '
          CALL GUESSEF(INFILE,ERRINFILE)
          ERRINFILE=
     +     INFILEX(21,ERRINFILE,NSCAN,NCHAN,STWV,DISP,21,.TRUE.) !........match
        END IF
C
        WRITE(*,100)'Output file name'
        OUTFILE=OUTFILEX(30,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        IF(CERR.EQ.'y')THEN
          WRITE(*,100)'Output error file name '
          CALL GUESSEF(OUTFILE,ERROUTFILE)
          ERROUTFILE=
     +     OUTFILEX(31,ERROUTFILE,NSCAN,NCHAN,STWV,DISP,1,.TRUE.) !.......match
        END IF
        
C
        WRITE(*,100)'Radial velocity (km/sec)'
        RADVEL=READF('@')
C
        WRITE(*,100)'Wait...'
        DO I=1,NSCAN
          READ(20) (S(J),J=1,NCHAN)
          CALL RVREBIN(RADVEL,NCHAN,S,SS,STWV,DISP)
          WRITE(30) (SS(J),J=1,NCHAN)
          IF(CERR.EQ.'y')THEN
            READ(21) (S(J),J=1,NCHAN)
            CALL RVREBIN(RADVEL,NCHAN,S,SS,STWV,DISP)
            WRITE(31) (SS(J),J=1,NCHAN)
          END IF
        END DO
        WRITE(*,101)'   ...OK!'
C------------------------------------------------------------------------------
        CLOSE(20)
        CLOSE(30)
        IF(CERR.EQ.'y')THEN
          CLOSE(21)
          CLOSE(31)
        END IF
C------------------------------------------------------------------------------
        STOP
100     FORMAT(A,$)
101     FORMAT(A)
C
        END
