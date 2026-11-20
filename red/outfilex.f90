!------------------------------------------------------------------------------
! Version 20-November-2026                                   File: outfilex.f90
!------------------------------------------------------------------------------
! Copyright N. Cardiel & J. Gorgas, Departamento de Astrofisica
! Universidad Complutense de Madrid, 28040-Madrid, Spain
! E-mail: ncl@astrax.fis.ucm.es or fjg@astrax.fis.ucm.es
!------------------------------------------------------------------------------
! This routine is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the Free
! Software Foundation; either version 2 of the License, or (at your option) any
! later version. See the file gnu-public-license.txt for details.
!------------------------------------------------------------------------------
!Comment
!
! CHARACTER*255 FUNCTION OUTFILEX(NF,FILENAME,NSCAN,NCHAN,STWV,DISP,MODE,LERR)
!
! Input: NF,FILENAME,NSCAN,NCHAN,STWV,DISP,MODE,LERR
! Input (COMMON): AIRMASS,TIMEXPOS,OBJECT,FITSFILE,COMMENT
! Output: OUTFILEX
!
! This function opens an output file (if the file does NOT exist), writing 
! the header keywords (usually declared in a previous call to INFILEX). 
! The function returns the name of the file opened. This function does NOT 
! write the data records (this action must be performed after a call to this 
! function).
!
! INTEGER       NF -> logical unit number of the file to be opened
! CHARACTER*(*) FILENAME -> default file name ('@' means there is not default)
! INTEGER       NSCAN -> no. of scans (pixels in the spatial direction)
! INTEGER       NCHAN -> no. of channels (pixels in the wavelength direction)
! REAL          STWV -> central wavelength of the first pixel
! REAL          DISP -> dispersion (Angstroms/pixel)
! INTEGER       MODE -> indicates the expected file format:
!               MODE=1: unformatted (with full header), i.e. REDUCEME format
!               MODE=2: unformatted (without header)
!               MODE=3: formatted (ascii file)
!               MODE=11,12,13: like MODE=1,2,3 but, if the file exist, it
!               is opened directly (without prompting)
! LOGICAL       LERR -> if .TRUE. the input file corresponds to an error frame;
!               if .FALSE. the input file does not corresponds to an error
!               frame (only when MODE=1, 11 or 21; otherwise it has no effect)
!
! Apart from NSCAN, NCHAN, STWV, and DISP, other global variables (declared
! through COMMON blocks in redlib.inc), are also saved: AIRMASS,
! TIMEXPOS, OBJECT, FITSFILE and COMMENT.
!
!Comment
!
!
! LERR=.TRUE.: el fichero a salvar es una imagen de errores. En este caso
!              anhadimos al nombre del objeto la terminacion " @ERROR@"
! LERR=.FALSE.: el fichero a salvar no es una imagen de errores
!
! El nombre del fichero se lee con READC_FILE: si el paquete se ha compilado
! con readline, TAB autocompleta el nombre del fichero. Si el nombre contiene
! "*" o "?", READC_FILE muestra los ficheros que coinciden (llamada a "ls")
! y vuelve a preguntar.
!------------------------------------------------------------------------------
        CHARACTER*255 FUNCTION&
     &   OUTFILEX(NF,FILENAME,NSCAN,NCHAN,STWV,DISP,MODE,LERR)
        USE REDLIB_M
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INTEGER TRUEBEG,TRUELEN
        CHARACTER*255 READC_FILE
!
        INTEGER NF
        INTEGER MODE
        INTEGER L1,L2
        INTEGER NERR,NMODE
        CHARACTER*(*) FILENAME
        LOGICAL LPROMPT,LOGFILE,LERR,LDIRECT
!------------------------------------------------------------------------------
! comprobamos si la imagen tiene que ser abierta directamente
        NMODE=MODE
        LDIRECT=.FALSE.
        IF(NMODE.GT.10)THEN
          NMODE=NMODE-10
          LDIRECT=.TRUE.
        END IF
!
        IF((NMODE.GT.3).OR.(NMODE.LT.1))THEN
          WRITE(*,*)
          WRITE(*,110)'MODE=',NMODE
          WRITE(*,101)'FATAL ERROR: invalid MODE value in function '//&
     &     'OUTFILEX.'
          STOP
        END IF
!
        LPROMPT=.FALSE.
        IF(FILENAME(1:1).EQ.'@')THEN
          L1=0
          L2=0
        ELSE
          L1=TRUEBEG(FILENAME)
          IF(L1.NE.0)THEN
            L2=TRUELEN(FILENAME)
            LPROMPT=.TRUE.
          ELSE
            L2=0
          END IF
        END IF
        NERR=0
!
        IF(LDIRECT)THEN
          IF(.NOT.LPROMPT)THEN
            WRITE(*,101)'FATAL ERROR: no file name given.'
            STOP
          END IF
          OUTFILEX=FILENAME(L1:L2)
          GOTO 11
        END IF
10      IF(LPROMPT)THEN
          OUTFILEX(1:255)=READC_FILE(FILENAME(L1:L2),'@')
        ELSE
            OUTFILEX(1:255)=READC_FILE('@','@')
        END IF
11      INQUIRE(FILE=OUTFILEX,EXIST=LOGFILE)
        IF(LOGFILE)THEN
          WRITE(*,101)'ERROR: This file already exist. You cannot '//&
     &     'overwrite it. Try again.'
!cc          IF(.NOT.LPROMPT) WRITE(*,100)'? '
          NERR=NERR+1
          IF(NERR.GT.10) STOP 'FATAL ERROR: too many errors.'
          GOTO 10
        ELSE
          IF(NMODE.EQ.1)THEN
            OPEN(NF,FILE=OUTFILEX,STATUS='NEW',FORM='UNFORMATTED')
            WRITE(NF)'abcdefghijkl'
            WRITE(NF) NSCAN,NCHAN
            WRITE(NF) STWV,DISP
            WRITE(NF) AIRMASS
            WRITE(NF) TIMEXPOS
            L1=TRUELEN(OBJECT)
            IF(LERR)THEN                     !el fichero a salvar es de errores
              WRITE(NF) L1+8
              IF(L1.GT.0)THEN
                WRITE(NF) OBJECT(1:L1)//' @ERROR@'
              ELSE
                WRITE(NF) ' @ERROR@'
              END IF
            ELSE                          !el fichero a salvar no es de errores
              WRITE(NF) L1
              IF(L1.GT.0)THEN
                WRITE(NF) OBJECT(1:L1)
              END IF
            END IF
            L1=TRUELEN(FITSFILE)
            WRITE(NF) L1
            IF(L1.GT.0)THEN
              WRITE(NF) FITSFILE(1:L1)
            END IF
            L1=TRUELEN(COMMENT)
            WRITE(NF) L1
            IF(L1.GT.0)THEN
              WRITE(NF) COMMENT(1:L1)
            END IF
          ELSEIF(NMODE.EQ.2)THEN
            OPEN(NF,FILE=OUTFILEX,STATUS='NEW',FORM='UNFORMATTED')
          ELSEIF(NMODE.EQ.3)THEN
            OPEN(NF,FILE=OUTFILEX,STATUS='NEW',FORM='FORMATTED')
          END IF
        END IF
!cc100     FORMAT(A,$)
101     FORMAT(A)
110     FORMAT(A,I6)
        END
