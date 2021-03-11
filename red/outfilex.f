C------------------------------------------------------------------------------
C Version 07-September-2007                                    File: outfilex.f
C------------------------------------------------------------------------------
C Copyright N. Cardiel & J. Gorgas, Departamento de Astrofisica
C Universidad Complutense de Madrid, 28040-Madrid, Spain
C E-mail: ncl@astrax.fis.ucm.es or fjg@astrax.fis.ucm.es
C------------------------------------------------------------------------------
C This routine is free software; you can redistribute it and/or modify it
C under the terms of the GNU General Public License as published by the Free
C Software Foundation; either version 2 of the License, or (at your option) any
C later version. See the file gnu-public-license.txt for details.
C------------------------------------------------------------------------------
Comment
C
C CHARACTER*255 FUNCTION OUTFILEX(NF,FILENAME,NSCAN,NCHAN,STWV,DISP,MODE,LERR)
C
C Input: NF,FILENAME,NSCAN,NCHAN,STWV,DISP,MODE,LERR
C Input (COMMON): AIRMASS,TIMEXPOS,OBJECT,FITSFILE,COMMENT
C Output: OUTFILEX
C
C This function opens an output file (if the file does NOT exist), writing 
C the header keywords (usually declared in a previous call to INFILEX). 
C The function returns the name of the file opened. This function does NOT 
C write the data records (this action must be performed after a call to this 
C function).
C
C INTEGER       NF -> logical unit number of the file to be opened
C CHARACTER*(*) FILENAME -> default file name ('@' means there is not default)
C INTEGER       NSCAN -> no. of scans (pixels in the spatial direction)
C INTEGER       NCHAN -> no. of channels (pixels in the wavelength direction)
C REAL          STWV -> central wavelength of the first pixel
C REAL          DISP -> dispersion (Angstroms/pixel)
C INTEGER       MODE -> indicates the expected file format:
C               MODE=1: unformatted (with full header), i.e. REDUCEME format
C               MODE=2: unformatted (without header)
C               MODE=3: formatted (ascii file)
C               MODE=11,12,13: like MODE=1,2,3 but, if the file exist, it
C               is opened directly (without prompting)
C LOGICAL       LERR -> if .TRUE. the input file corresponds to an error frame;
C               if .FALSE. the input file does not corresponds to an error
C               frame (only when MODE=1, 11 or 21; otherwise it has no effect)
C
C Apart from NSCAN, NCHAN, STWV, and DISP, other global variables (declared
C through COMMON blocks in redlib.inc), are also saved: AIRMASS,
C TIMEXPOS, OBJECT, FITSFILE and COMMENT.
C
Comment
C
C
C LERR=.TRUE.: el fichero a salvar es una imagen de errores. En este caso
C              anhadimos al nombre del objeto la terminacion " @ERROR@"
C LERR=.FALSE.: el fichero a salvar no es una imagen de errores
C
C Si el nombre del fichero contiene "*" o "?", la rutina realiza una llamada
C a la funcion SYSTEM para obtener un directorio con el nombre introducido.
C------------------------------------------------------------------------------
        CHARACTER*255 FUNCTION
     +   OUTFILEX(NF,FILENAME,NSCAN,NCHAN,STWV,DISP,MODE,LERR)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'futils.inc'
        INTEGER TRUEBEG,TRUELEN
C
        INTEGER SYSTEMFUNCTION
C
        INTEGER NF
        INTEGER MODE
        INTEGER L1,L2
        INTEGER NERR,NMODE
        INTEGER ISYSTEM
        CHARACTER*(*) FILENAME
        LOGICAL LPROMPT,LOGFILE,LERR,LDIRECT
        LOGICAL LOG1,LOG2
C------------------------------------------------------------------------------
C comprobamos si la imagen tiene que ser abierta directamente
        NMODE=MODE
        LDIRECT=.FALSE.
        IF(NMODE.GT.10)THEN
          NMODE=NMODE-10
          LDIRECT=.TRUE.
        END IF
C
        IF((NMODE.GT.3).OR.(NMODE.LT.1))THEN
          WRITE(*,*)
          WRITE(*,110)'MODE=',NMODE
          WRITE(*,101)'FATAL ERROR: invalid MODE value in function '//
     +     'OUTFILEX.'
          STOP
        END IF
C
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
C
        IF(LDIRECT)THEN
          IF(.NOT.LPROMPT)THEN
            WRITE(*,101)'FATAL ERROR: no file name given.'
            STOP
          END IF
          OUTFILEX=FILENAME(L1:L2)
          GOTO 11
        END IF
10      IF(LPROMPT)THEN
          OUTFILEX(1:255)=READC(FILENAME(L1:L2),'@')
        ELSE
            OUTFILEX(1:255)=READC('@','@')
        END IF
        LOG1=(INDEX(OUTFILEX,'*').NE.0)
        LOG2=(INDEX(OUTFILEX,'?').NE.0)
        IF(LOG1.OR.LOG2)THEN
           ISYSTEM=SYSTEMFUNCTION('ls '//OUTFILEX)
           GOTO 10
        END IF
11      INQUIRE(FILE=OUTFILEX,EXIST=LOGFILE)
        IF(LOGFILE)THEN
          WRITE(*,101)'ERROR: This file already exist. You cannot '//
     +     'overwrite it. Try again.'
ccc          IF(.NOT.LPROMPT) WRITE(*,100)'? '
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
ccc100     FORMAT(A,$)
101     FORMAT(A)
110     FORMAT(A,I6)
        END
