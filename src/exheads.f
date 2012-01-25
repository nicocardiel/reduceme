C------------------------------------------------------------------------------
C Version 6-December-1996                                       file: exheads.f
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
C Program: exheads
C Classification: examination & statistics
C Description: Shows image header information without prompting.
C
Comment
C
C Lee la cabecera de una imagen para obtener la informacion de su
C tamanho. Este programa no realiza ninguna modificacion sobre el
C fichero.
C
        PROGRAM EXHEADS
        IMPLICIT NONE
C Include NSMAX,NCMAX,NSCAN,NCHAN
        INCLUDE 'redlib.inc'
        INCLUDE 'futils.inc'
        INTEGER TRUELEN
C
        INTEGER L
        CHARACTER*12 CLAVE
        CHARACTER*50 CDUMMY
        CHARACTER*80 INFILE
        LOGICAL LOGFILE
C------------------------------------------------------------------------------
ccc        THISPROGRAM='exheads'
ccc     CALL WELCOME('6-December-1996')
C
ccc        WRITE(*,100)'File name to be examined'
        INFILE(1:80)=READC('@','@')
        INQUIRE(FILE=INFILE,EXIST=LOGFILE)
        IF(LOGFILE)THEN
          OPEN(20,FILE=INFILE,STATUS='OLD',FORM='UNFORMATTED')
        ELSE
          WRITE(*,101)'ERROR: this file does not exist.'
          STOP
        END IF
        READ(20,ERR=20)CLAVE
        IF(CLAVE.EQ.CLAVE_RED)THEN
          READ(20)NSCAN,NCHAN
          WRITE(CDUMMY,'(I10,A1,I10)')NSCAN,',',NCHAN
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          WRITE(*,101)'Image size (NSCAN,NCHAN): '//CDUMMY(1:L)
          READ(20) STWV,DISP
ccc          WRITE(CDUMMY,*)STWV
ccc          CALL RMBLANK(CDUMMY,CDUMMY,L)
ccc          WRITE(*,101)'STWV    : '//CDUMMY(1:L)
ccc          WRITE(CDUMMY,*)DISP
ccc          CALL RMBLANK(CDUMMY,CDUMMY,L)
ccc          WRITE(*,101)'DISP    : '//CDUMMY(1:L)
          READ(20)AIRMASS
ccc          WRITE(CDUMMY,*)AIRMASS
ccc          CALL RMBLANK(CDUMMY,CDUMMY,L)
ccc          WRITE(*,101)'Airmass : '//CDUMMY(1:L)
          READ(20)TIMEXPOS
ccc          WRITE(CDUMMY,*)TIMEXPOS
ccc          CALL RMBLANK(CDUMMY,CDUMMY,L)
ccc          WRITE(*,101)'Timexpos: '//CDUMMY(1:L)
          READ(20)L
          IF(L.GT.0)THEN
            READ(20)OBJECT(1:L)
            WRITE(*,101)'Object  : '//OBJECT(1:TRUELEN(OBJECT))
          ELSE
            WRITE(*,101)'Object  : [not found]'
          END IF
          READ(20)L 
          IF(L.GT.0)THEN
            READ(20)FITSFILE(1:L)
ccc            WRITE(*,101)'FITSfile: '//FITSFILE(1:TRUELEN(FITSFILE))
          ELSE
ccc            WRITE(*,101)'FITSfile: [not found]'
          END IF
          READ(20)L
          IF(L.GT.0)THEN
            READ(20)COMMENT(1:L)
ccc            WRITE(*,101)'Comment : '//COMMENT(1:TRUELEN(COMMENT))
          ELSE
ccc            WRITE(*,101)'Comment : [not found]'
          END IF
          GOTO 30
        END IF
C
20      CLOSE(20)
        WRITE(*,101)'>>> WARNING <<<'
        WRITE(*,101)'This file does not contain head information.'
        STOP
C
30      IF(NSCAN.GT.NSMAX)THEN
          WRITE(*,101)'ERROR: NSCAN.GT.NSMAX'
        END IF
        IF(NCHAN.GT.NCMAX)THEN
          WRITE(*,101)'ERROR: NCHAN.GT.NCMAX'
        END IF
        CLOSE(20)
        STOP
C
ccc100     FORMAT(A,$)
101     FORMAT(A)
        END
