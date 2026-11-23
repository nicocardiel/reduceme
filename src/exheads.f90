!------------------------------------------------------------------------------
! Version 24-November-2026                                    file: exheads.f90
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
! Program: exheads
! Classification: examination & statistics
! Description: Shows image header information without prompting.
!
!Comment
!
! Lee la cabecera de una imagen para obtener la informacion de su
! tamanho. Este programa no realiza ninguna modificacion sobre el
! fichero.
!
        PROGRAM EXHEADS
        USE RED_M
        USE FUTILS_M
        USE REDLIB_M
        IMPLICIT NONE
! Include NSMAX,NCMAX,NSCAN,NCHAN
        INCLUDE 'redlib.inc'
!
        INTEGER L,LBRACKET
        CHARACTER*12 CLAVE
        CHARACTER*50 CDUMMY
        CHARACTER*255 INFILE
        LOGICAL LOGFILE,LFITSOK
!------------------------------------------------------------------------------
!cc        THISPROGRAM='exheads'
!cc     CALL WELCOME('6-December-1996')
!
!cc        WRITE(*,100)'File name to be examined'
        INFILE=READC_FILE('@','@')
! FITS HDU specification (file.fits[n]): check the file without it
        L=TRUELEN(INFILE)
        LBRACKET=INDEX(INFILE(1:MAX(L,1)),'[')
        IF((L.GT.0).AND.(LBRACKET.GT.1).AND.&
     &   (INFILE(MAX(L,1):MAX(L,1)).EQ.']'))THEN
          INQUIRE(FILE=INFILE(1:LBRACKET-1),EXIST=LOGFILE)
        ELSE
          INQUIRE(FILE=INFILE,EXIST=LOGFILE)
        END IF
        IF(LOGFILE)THEN
          IF(ISFITSFILE(INFILE))THEN
! FITS file: converted into a scratch file in REDUCEME format
            CALL FITS2SCRATCH(20,INFILE,.FALSE.,LFITSOK)
            IF(.NOT.LFITSOK) STOP
          ELSE
            OPEN(20,FILE=INFILE,STATUS='OLD',FORM='UNFORMATTED')
          END IF
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
!cc          WRITE(CDUMMY,*)STWV
!cc          CALL RMBLANK(CDUMMY,CDUMMY,L)
!cc          WRITE(*,101)'STWV    : '//CDUMMY(1:L)
!cc          WRITE(CDUMMY,*)DISP
!cc          CALL RMBLANK(CDUMMY,CDUMMY,L)
!cc          WRITE(*,101)'DISP    : '//CDUMMY(1:L)
          READ(20)AIRMASS
!cc          WRITE(CDUMMY,*)AIRMASS
!cc          CALL RMBLANK(CDUMMY,CDUMMY,L)
!cc          WRITE(*,101)'Airmass : '//CDUMMY(1:L)
          READ(20)TIMEXPOS
!cc          WRITE(CDUMMY,*)TIMEXPOS
!cc          CALL RMBLANK(CDUMMY,CDUMMY,L)
!cc          WRITE(*,101)'Timexpos: '//CDUMMY(1:L)
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
!cc            WRITE(*,101)'FITSfile: '//FITSFILE(1:TRUELEN(FITSFILE))
          ELSE
!cc            WRITE(*,101)'FITSfile: [not found]'
          END IF
          READ(20)L
          IF(L.GT.0)THEN
            READ(20)COMMENT(1:L)
!cc            WRITE(*,101)'Comment : '//COMMENT(1:TRUELEN(COMMENT))
          ELSE
!cc            WRITE(*,101)'Comment : [not found]'
          END IF
          GOTO 30
        END IF
!
20      CLOSE(20)
        WRITE(*,101)'>>> WARNING <<<'
        WRITE(*,101)'This file does not contain head information.'
        STOP
!
30      IF(NSCAN.GT.NSMAX)THEN
          WRITE(*,101)'ERROR: NSCAN.GT.NSMAX'
        END IF
        IF(NCHAN.GT.NCMAX)THEN
          WRITE(*,101)'ERROR: NCHAN.GT.NCMAX'
        END IF
        CLOSE(20)
        STOP
!
!cc100     FORMAT(A,$)
101     FORMAT(A)
        END
