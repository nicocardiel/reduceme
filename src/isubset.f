C------------------------------------------------------------------------------
C Version 13-December-1996                                      file: isubset.f
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
C Program: isubset
C Classification: arithmetic & manipulations
C Description: Produces a subset of an image.
C
Comment
C
C Extrae un SUBSET de una imagen.
C
        PROGRAM ISUBSET
C
        IMPLICIT NONE
C Include NSMAX,NCMAX,NSCAN,NCHAN
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
C
        INTEGER I,J
        INTEGER NS1,NS2,NC1,NC2
        REAL A(NCMAX,NSMAX),ERR(NCMAX,NSMAX)
        REAL NEWSTWV
        CHARACTER*1 CERR,CUPDATE
        CHARACTER*75 INFILE,ERRFILE,OUTFILE
C------------------------------------------------------------------------------
        THISPROGRAM='isubset'
        CALL WELCOME('13-December-1996')
C
        WRITE(*,100)'Work with error images (y/n) '
        CERR(1:1)=READC('y','yn')
C
        WRITE(*,100)'Input file name'
        INFILE=INFILEX(30,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        DO I=1,NSCAN
          READ(30) (A(J,I),J=1,NCHAN)
        END DO
        CLOSE(30)
C
        IF(CERR.EQ.'y')THEN
          WRITE(*,100)'Input error file name '
          CALL GUESSEF(INFILE,ERRFILE)
          INFILE=INFILEX(32,ERRFILE,NSCAN,NCHAN,STWV,DISP,21,.TRUE.) !match
          DO I=1,NSCAN
            READ(32) (ERR(J,I),J=1,NCHAN)
          END DO
          CLOSE(32)
        END IF
C------------------------------------------------------------------------------
10      WRITE(*,100)'First and last scan   '
        CALL READ2I('@',NS1,NS2)
        IF((NS1.LT.1).OR.(NS2.GT.NSCAN))THEN
          WRITE(*,101)'ERROR: number(s) out of range. Try again.'
          GOTO 10
        END IF
        IF(NS1.GT.NS2)THEN
          WRITE(*,101)'ERROR: wrong order. Try again.'
          GOTO 10
        END IF
20      WRITE(*,100)'First and last channel'
        CALL READ2I('@',NC1,NC2)
        IF((NC1.LT.1).OR.(NC2.GT.NCHAN))THEN
          WRITE(*,101)'ERROR: number(s) out of range. Try again.'
          GOTO 20
        END IF
        IF(NC1.GT.NC2)THEN
          WRITE(*,101)'ERROR: wrong order. Try again.'
          GOTO 20
        END IF
C------------------------------------------------------------------------------
        IF((STWV.NE.0.0).AND.(DISP.NE.0.0))THEN
          IF(NC1.NE.1)THEN
            NEWSTWV=STWV+REAL(NC1-1)*DISP
            WRITE(*,100)'>>> Calculated STWV for subimage: '
            WRITE(*,*)NEWSTWV
            WRITE(*,100)'Update STWV (y/n) '
            CUPDATE(1:1)=READC('y','yn')
            IF(CUPDATE.EQ.'y') STWV=NEWSTWV
          END IF
        END IF
C------------------------------------------------------------------------------
        WRITE(*,100)'Output file name'
        OUTFILE=OUTFILEX(40,'@',NS2-NS1+1,NC2-NC1+1,STWV,DISP,1,.FALSE.)
        DO I=NS1,NS2
          WRITE(40) (A(J,I),J=NC1,NC2)
        END DO
        CLOSE(40)
C
        IF(CERR.EQ.'y')THEN
          WRITE(*,100)'Output error file name '
          CALL GUESSEF(OUTFILE,ERRFILE)
          OUTFILE=OUTFILEX(42,ERRFILE,NS2-NS1+1,NC2-NC1+1,
     +     STWV,DISP,1,.TRUE.)
          DO I=NS1,NS2
            WRITE(42) (ERR(J,I),J=NC1,NC2)
          END DO
          CLOSE(42)
        END IF
C------------------------------------------------------------------------------
        STOP
100     FORMAT(A,$)
101     FORMAT(A)
        END
