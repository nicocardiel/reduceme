!------------------------------------------------------------------------------
! Version 24-November-2026                                    File: infilex.f90
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
! CHARACTER*255 FUNCTION INFILEX(NF,FILENAME,NSCAN,NCHAN,STWV,DISP,MODE,LERR)
!
! Input: NF,FILENAME,MODE,LERR  (and NSCAN,NCHAN,STWV,DISP if MODE=21)
! Output: NSCAN,NCHAN,STWV,DISP
! Output (COMMON): AIRMASS,TIMEXPOS,OBJECT,FITSFILE,COMMENT
!
! This function opens an input file (if the file exist), reading the header 
! keywords, and verifying whether the required file has REDUCEME format and 
! that the image dimensions do not exceed the maximum values (NCMAX,NSMAX,
! both 65536, defined in the module REDLIB_M). The function returns the name
! of the file opened. This function does NOT read the data records (this
! action must be performed after a call to this function).
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
!               MODE=4: formatted, with RECL=2880 to read FITS files
!               MODE=11,12,13,14: like MODE=1,2,3,4 but, if the file exist, it
!               is opened directly (without prompting)
!               MODE=21: like MODE=1 but the function verifies whether the
!               input values of NSCAN,NCHAN,STWV, and DISP are identical with
!               those in the header of the file.
! LOGICAL       LERR -> if .TRUE. the input file corresponds to an error frame;
!               if .FALSE. the input file does not corresponds to an error
!               frame (only when MODE=1, 11 or 21; otherwise it has no effect)
!
! Apart from NSCAN, NCHAN, STWV, and DISP, other global variables (declared
! through COMMON blocks in redlib.inc), are also (re)declared: AIRMASS,
! TIMEXPOS, OBJECT, FITSFILE and COMMENT.
!
! FITS files (MODE=1, 11 or 21): if the file is a FITS file (see ISFITSFILE;
! the CFITSIO syntax file.fits[n] selects an HDU), FITS2SCRATCH converts it
! into a scratch file in REDUCEME format connected to the unit NF, which is
! then read as any REDUCEME file: the calling program reads the data records
! as usual.
!
!Comment
!
!------------------------------------------------------------------------------
! Si el fichero tiene formato libre conteniendo varias columnas, MODE retorna
! negativo con el numero de columnas.
!
! LERR no tiene ningun efecto para MODE=2,3,4  (12,13,14)
!
! Para distinguir si una imagen es de errores o no, anhadimos al nombre del
! objeto 8 caracteres: ' @ERROR@'
!                       12345678
!
! IMPORTANTE:
! Si LERR=.TRUE. y el programa no encuentra ' @ERROR@', aparece un aviso y
! se pide al usuario que confirme continuar.
! Si LERR=.TRUE. y el programa encuentra ' @ERROR@', se elimina esta tira
! de caracteres de la variable global OBJECT (dado que luego se anhadira
! a la hora de salvar la imagen con OUTFILEX).
! Si LERR=.FALSE. y el program encuentra ' @ERROR@', aparece un aviso y
! se pide al usuario que confirme continuar. Asimismo, esta cadena NO es
! borrada en este caso de la variable global OBJECT, de modo que al ser
! salvada con OUTFILEX, lo normal es que se salve como si no fuera imagen
! de errores (por lo que no se le anhadiria la cadena ' @ERROR@' y se 
! perderia el caracter de imagen de errores), pero al contener la cadena
! de control en la variable OBJECT, el fichero salvado sigue siendo de 
! errores.
! El nombre del fichero se lee con READC_FILE: si el paquete se ha compilado
! con readline, TAB autocompleta el nombre del fichero. Si el nombre contiene
! "*" o "?", READC_FILE muestra los ficheros que coinciden (llamada a "ls")
! y vuelve a preguntar.
!------------------------------------------------------------------------------
        CHARACTER*255 FUNCTION& 
     &   INFILEX(NF,FILENAME,NSCAN,NCHAN,STWV,DISP,MODE,LERR)
        USE REDLIB_M
        IMPLICIT NONE
! Include NSMAX,NCMAX,NSCAN,NCHAN
        INCLUDE 'redlib.inc'
!
        INTEGER NF
        INTEGER MODE
        CHARACTER*(*) FILENAME
        LOGICAL LERR
!
        INTEGER I
        INTEGER L1,L2,L0,LCOMA
        INTEGER NMODE
        INTEGER NERR
        INTEGER NSKIP
        INTEGER NSCAN2,NCHAN2
        REAL STWV2,DISP2
        CHARACTER*1 COPC,CFOR,CNM,CCONT
        CHARACTER*12 CLAVE
        CHARACTER*255 CDUMMY
        INTEGER LBRACKET
        LOGICAL LPROMPT,LOGFILE,LDIRECT,LMATCH,LASK,LFITSOK
        CHARACTER*255 OBJECT2
!------------------------------------------------------------------------------
! determinamos si hay que buscar coincidencia en los parametros fundamentales
! de la imagen que se quiere abrir
        LMATCH=.FALSE.
        NMODE=MODE
        IF(MODE.EQ.21)THEN
          NMODE=1
          LMATCH=.TRUE.
        END IF
!
        LDIRECT=.FALSE.                          !abrir el fichero directamente
        IF(NMODE.GT.10)THEN
          NMODE=MODE-10
          LDIRECT=.TRUE.
        END IF
!
        IF((NMODE.GT.4).OR.(NMODE.LT.1))THEN
          WRITE(*,*)
          WRITE(*,110)'MODE=',NMODE
          WRITE(*,101)'FATAL ERROR: invalid MODE value in '//&
     &     'function INFILEX.'
          STOP
        END IF
!
        LPROMPT=.FALSE.
        IF(ICHAR(FILENAME(1:1)).EQ.64)THEN
          L1=0
          L2=0
        ELSE
          L1=TRUEBEG(FILENAME)
          IF(L1.NE.0)THEN
            LCOMA=INDEX(FILENAME,',')
            IF(LCOMA.NE.0)THEN
              L2=LCOMA-1
            ELSE
              L2=TRUELEN(FILENAME)
            END IF
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
          INFILEX=FILENAME(L1:L2)
          GOTO 11
        END IF
10      IF(LPROMPT)THEN
          INFILEX(1:255)=READC_FILE(FILENAME(L1:L2),'@')
        ELSE
          INFILEX(1:255)=READC_FILE('@','@')
        END IF
11      LCOMA=INDEX(INFILEX,',') !permitimos poner una coma detras del nombre
        IF(LCOMA.GT.1)THEN
          INQUIRE(FILE=INFILEX(1:LCOMA-1),EXIST=LOGFILE)
        ELSE
! FITS HDU specification (file.fits[n]): check the file without it
          L0=TRUELEN(INFILEX)
          LBRACKET=INDEX(INFILEX(1:MAX(L0,1)),'[')
          IF((L0.GT.0).AND.(LBRACKET.GT.1).AND.&
     &     (INFILEX(MAX(L0,1):MAX(L0,1)).EQ.']'))THEN
            INQUIRE(FILE=INFILEX(1:LBRACKET-1),EXIST=LOGFILE)
          ELSE
            INQUIRE(FILE=INFILEX,EXIST=LOGFILE)
          END IF
        END IF
        IF(.NOT.LOGFILE)THEN
          WRITE(*,101)'ERROR: This file does not exist. Try again.'
          NERR=NERR+1
          IF(NERR.GT.10) STOP 'FATAL ERROR: too many errors.'
          GOTO 10
        END IF
!
        IF(NMODE.GT.1) GOTO 40
!
        LCOMA=INDEX(INFILEX,',')
        IF(LCOMA.GT.1)THEN
          CDUMMY=INFILEX(1:LCOMA-1)
        ELSE
          CDUMMY=INFILEX
        END IF
        IF(ISFITSFILE(CDUMMY))THEN
! FITS file: converted into a scratch file in REDUCEME format (unit NF)
          CALL FITS2SCRATCH(NF,CDUMMY,LERR,LFITSOK)
          IF(.NOT.LFITSOK)THEN
            IF(LDIRECT) STOP 'FATAL ERROR: the FITS file cannot be read.'
            NERR=NERR+1
            IF(NERR.GT.10) STOP 'FATAL ERROR: too many errors.'
            WRITE(*,101)'Enter new file name:'
            GOTO 10
          END IF
        ELSE
          OPEN(NF,FILE=CDUMMY,STATUS='OLD',FORM='UNFORMATTED',ERR=999)
        END IF
        READ(NF,ERR=20)CLAVE
        IF(CLAVE.EQ.CLAVE_RED)THEN
!..............................................................................
          READ(NF)NSCAN2,NCHAN2
          WRITE(*,110)'>>> NSCAN : ',NSCAN2
          WRITE(*,110)'>>> NCHAN : ',NCHAN2
!..............................................................................
          READ(NF)STWV2,DISP2
          WRITE(CDUMMY,*)STWV2
          L1=TRUEBEG(CDUMMY)
          L2=TRUELEN(CDUMMY)
          WRITE(*,101)'>>> STWV  : '//CDUMMY(L1:L2)
          WRITE(CDUMMY,*)DISP2
          L1=TRUEBEG(CDUMMY)
          L2=TRUELEN(CDUMMY)
          WRITE(*,101)'>>> DISP  : '//CDUMMY(L1:L2)
!..............................................................................
          READ(NF)AIRMASS
          READ(NF)TIMEXPOS
!..............................................................................
          READ(NF)L1
          IF(L1.GT.0)THEN
            READ(NF)OBJECT2(1:L1)
            WRITE(*,101)'>>> OBJECT: '//OBJECT2(1:L1)
          ELSE
            WRITE(*,101)'>>> OBJECT: [not found]'
          END IF
          IF(LERR)THEN                             !esperamos imagen de errores
            LASK=.FALSE.
            IF(L1.LT.7)THEN
              LASK=.TRUE.
            ELSE
              LASK=(OBJECT2(L1-6:L1).NE.'@ERROR@')
            END IF
            IF(LASK)THEN
              WRITE(*,101)'WARNING: Expected error file has not '//&
     &         'been found.'
              WRITE(*,101)'Current file name: '//&
     &         INFILEX(1:TRUELEN(INFILEX))
              WRITE(*,100,ADVANCE='NO')&
     &           'Do you want to continue with this '//&
     &         'file (y/n) '
              CCONT(1:1)=READC('n','yn')
              IF(CCONT.EQ.'n')THEN
                CLOSE(NF)
                WRITE(*,101)'Enter new file name:'
                GOTO 10
              END IF
            ELSE                      !si la imagen de errores se ha encontrado
              L1=L1-7                         !eliminamos ' @ERROR@' del nombre
            END IF
          ELSE                                  !no esperamos imagen de errores
            LASK=.FALSE.
            IF(L1.GE.7)THEN
              LASK=(OBJECT2(L1-6:L1).EQ.'@ERROR@')
            END IF
            IF(LASK)THEN
              WRITE(*,101)'WARNING: Unexpected error file has '//&
     &         'been found.'
              WRITE(*,101)'Current file name: '//&
     &         INFILEX(1:TRUELEN(INFILEX))
              WRITE(*,100,ADVANCE='NO')&
     &           'Do you want to continue with this '//&
     &         'file (y/n) '
              CCONT(1:1)=READC('n','yn')
              IF(CCONT.EQ.'n')THEN
                CLOSE(NF)
                WRITE(*,101)'Enter new file name:'
                GOTO 10
              END IF
            END IF
          END IF
          IF(L1+1.LE.255)THEN
            DO I=L1+1,255
              OBJECT2(I:I)=CHAR(32)
            END DO
          END IF
!..............................................................................
          READ(NF)L0
          IF(L0.GT.0) READ(NF)FITSFILE(1:L0)
          IF(L0+1.LE.255)THEN
            DO I=L0+1,255
              FITSFILE(I:I)=CHAR(32)
            END DO
          END IF
!..............................................................................
          READ(NF)L0
          IF(L0.GT.0) READ(NF)COMMENT(1:L0)
          IF(L0+1.LE.255)THEN
            DO I=L0+1,255
              COMMENT(I:I)=CHAR(32)
            END DO
          END IF
!..............................................................................
          IF(LMATCH)THEN                 !comparamos si los valores son iguales
            IF(NSCAN.NE.NSCAN2)THEN
              WRITE(*,100,ADVANCE='NO')&
     &           'WARNING: NSCAN does not match requested '//&
     &         'value: '
              WRITE(*,*)NSCAN
              WRITE(*,100,ADVANCE='NO')'Continue (y/n) '
              CNM(1:1)=READC('n','yn')
              IF(CNM.EQ.'n')THEN
                CLOSE(NF)
                STOP
              END IF
              NSCAN=NSCAN2
            END IF
            IF(NCHAN.NE.NCHAN2)THEN
              WRITE(*,100,ADVANCE='NO')&
     &           'WARNING: NCHAN does not match requested '//&
     &         'value: '
              WRITE(*,*)NCHAN
              WRITE(*,100,ADVANCE='NO')'Continue (y/n) '
              CNM(1:1)=READC('n','yn')
              IF(CNM.EQ.'n')THEN
                CLOSE(NF)
                STOP
              END IF
              NCHAN=NCHAN2
            END IF
            IF(STWV.NE.STWV2)THEN
              WRITE(*,100,ADVANCE='NO')&
     &           'WARNING: STWV does not match requested '//&
     &         'value: '
              WRITE(*,*)STWV
              WRITE(*,100,ADVANCE='NO')'Continue (y/n) '
              CNM(1:1)=READC('n','yn')
              IF(CNM.EQ.'n')THEN
                CLOSE(NF)
                STOP
              END IF
              STWV=STWV2
            END IF
            IF(DISP.NE.DISP2)THEN
              WRITE(*,100,ADVANCE='NO')&
     &           'WARNING: DISP does not match requested '//&
     &         'value: '
              WRITE(*,*)DISP
              WRITE(*,100,ADVANCE='NO')'Continue (y/n) '
              CNM(1:1)=READC('n','yn')
              IF(CNM.EQ.'n')THEN
                CLOSE(NF)
                STOP
              END IF
              DISP=DISP2
            END IF
            IF(OBJECT2(1:L1).NE.OBJECT(1:TRUELEN(OBJECT)))THEN
              WRITE(*,101)'WARNING: OBJECT does not match requested '//&
     &         'value:'
              WRITE(*,101)'>>> '//OBJECT(1:TRUELEN(OBJECT))
!cc           WRITE(*,100)'Continue (y/n) '
!cc           CNM(1:1)=READC('n','yn')
              CNM='y'
              IF(CNM.EQ.'n')THEN
                CLOSE(NF)
                STOP
              END IF
              OBJECT=OBJECT2
            END IF
          ELSE
            NSCAN=NSCAN2
            NCHAN=NCHAN2
            STWV=STWV2
            DISP=DISP2
            OBJECT=OBJECT2
          END IF
          GOTO 30
!..............................................................................
        END IF
!
20      CLOSE(NF)
        WRITE(*,101)'>>> WARNING <<<'
        WRITE(*,101)'This file does not contain head information.'
        WRITE(*,100,ADVANCE='NO')'Do you want to read it anyway (y/n) '
        COPC(1:1)=READC('y','yn')
        IF(COPC.EQ.'y')THEN
          WRITE(*,101)'(1) unformatted file'
          WRITE(*,101)'(2) formatted file (multiple column)'
          WRITE(*,100,ADVANCE='NO')'Option'
          CFOR(1:1)=READC('@','12')
          IF(CFOR.EQ.'1')THEN
            OPEN(NF,FILE=INFILEX,STATUS='OLD',FORM='UNFORMATTED')
            WRITE(*,100,ADVANCE='NO')'NSCAN'
            NSCAN=READI('@')
            WRITE(*,100,ADVANCE='NO')'NCHAN'
            NCHAN=READI('@')
            WRITE(*,100,ADVANCE='NO')'STWV'
            STWV=READF('0.0')
            WRITE(*,100,ADVANCE='NO')'DISP'
            DISP=READF('0.0')
            AIRMASS=0.
            TIMEXPOS=0.
            OBJECT=CHAR(32)
            FITSFILE=CHAR(32)
            COMMENT=CHAR(32)
          ELSEIF(CFOR.EQ.'2')THEN
            WRITE(*,100,ADVANCE='NO')'No. of lines to be skipped'
            NSKIP=READI('0')
            OPEN(NF,FILE=INFILEX,STATUS='OLD',FORM='FORMATTED')
            DO I=1,NSKIP
              READ(NF,*)
            END DO
            WRITE(*,100,ADVANCE='NO')'No. of columns'
            MODE=-READI('2')
          ELSE
            WRITE(*,101)'FATAL ERROR in INFILEX!!!'
            STOP
          END IF
        ELSE
          GOTO 10
        END IF
!
30      IF(NSCAN.GT.NSMAX)THEN
          WRITE(*,*)
          WRITE(*,100,ADVANCE='NO')'NSCAN, NSMAX: '
          WRITE(*,*) NSCAN,NSMAX
          WRITE(*,101)'FATAL ERROR: NSCAN.GT.NSMAX'
          CLOSE(NF)
          STOP
        END IF
        IF(NCHAN.GT.NCMAX)THEN
          WRITE(*,*)
          WRITE(*,100,ADVANCE='NO')'NCHAN, NCMAX: '
          WRITE(*,*) NCHAN,NCMAX
          WRITE(*,101)'FATAL ERROR: NCHAN.GT.NCMAX'
          CLOSE(NF)  
          STOP
        END IF
        RETURN
!
40      IF(NMODE.EQ.2)THEN
          OPEN(NF,FILE=INFILEX,STATUS='OLD',FORM='UNFORMATTED',ERR=999)
        ELSEIF(NMODE.EQ.3)THEN
          OPEN(NF,FILE=INFILEX,STATUS='OLD',FORM='FORMATTED',ERR=999)
        ELSEIF(NMODE.EQ.4)THEN
!cc          OPEN(NF,FILE=INFILEX,STATUS='OLD',READONLY,FORM='FORMATTED',
!cc     +     RECORDTYPE='FIXED',RECL=2880,ACCESS='DIRECT',ERR=999)
          OPEN(NF,FILE=INFILEX,STATUS='OLD',FORM='FORMATTED',&
     &     RECL=2880,ACCESS='DIRECT',ERR=999)
        END IF
!
        RETURN
!
999     WRITE(*,101)'FATAL ERROR: opening the file: '//&
     &   INFILEX(1:TRUELEN(INFILEX))
        STOP
100     FORMAT(A)
101     FORMAT(A)
110     FORMAT(A,I6)
        END
