C------------------------------------------------------------------------------
C Version 04-October-2007                                       File: infilex.f
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
C CHARACTER*75 FUNCTION INFILEX(NF,FILENAME,NSCAN,NCHAN,STWV,DISP,MODE,LERR)
C
C Input: NF,FILENAME,MODE,LERR  (and NSCAN,NCHAN,STWV,DISP if MODE=21)
C Output: NSCAN,NCHAN,STWV,DISP
C Output (COMMON): AIRMASS,TIMEXPOS,OBJECT,FITSFILE,COMMENT
C
C This function opens an input file (if the file exist), reading the header 
C keywords, and verifying whether the required file has REDUCEME format and 
C that the image dimensions do not exceed the maximum expected values (defined 
C in NCMAX,NSMAX --see the file redlib.inc--). The function returns the name
C of the file opened. This function does NOT read the data records (this
C action must be performed after a call to this function).
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
C               MODE=4: formatted, with RECL=2880 to read FITS files
C               MODE=11,12,13,14: like MODE=1,2,3,4 but, if the file exist, it
C               is opened directly (without prompting)
C               MODE=21: like MODE=1 but the function verifies whether the
C               input values of NSCAN,NCHAN,STWV, and DISP are identical with
C               those in the header of the file.
C LOGICAL       LERR -> if .TRUE. the input file corresponds to an error frame;
C               if .FALSE. the input file does not corresponds to an error
C               frame (only when MODE=1, 11 or 21; otherwise it has no effect)
C
C Apart from NSCAN, NCHAN, STWV, and DISP, other global variables (declared
C through COMMON blocks in redlib.inc), are also (re)declared: AIRMASS,
C TIMEXPOS, OBJECT, FITSFILE and COMMENT.
C
Comment
C
C------------------------------------------------------------------------------
C Si el fichero tiene formato libre conteniendo varias columnas, MODE retorna
C negativo con el numero de columnas.
C
C LERR no tiene ningun efecto para MODE=2,3,4  (12,13,14)
C
C Para distinguir si una imagen es de errores o no, anhadimos al nombre del
C objeto 8 caracteres: ' @ERROR@'
C                       12345678
C
C IMPORTANTE:
C Si LERR=.TRUE. y el programa no encuentra ' @ERROR@', aparece un aviso y
C se pide al usuario que confirme continuar.
C Si LERR=.TRUE. y el programa encuentra ' @ERROR@', se elimina esta tira
C de caracteres de la variable global OBJECT (dado que luego se anhadira
C a la hora de salvar la imagen con OUTFILEX).
C Si LERR=.FALSE. y el program encuentra ' @ERROR@', aparece un aviso y
C se pide al usuario que confirme continuar. Asimismo, esta cadena NO es
C borrada en este caso de la variable global OBJECT, de modo que al ser
C salvada con OUTFILEX, lo normal es que se salve como si no fuera imagen
C de errores (por lo que no se le anhadiria la cadena ' @ERROR@' y se 
C perderia el caracter de imagen de errores), pero al contener la cadena
C de control en la variable OBJECT, el fichero salvado sigue siendo de 
C errores.
C Si el nombre del fichero contiene "*" o "?", la rutina realiza una llamada
C a la funcion SYSTEM para obtener un directorio con el nombre introducido.
C------------------------------------------------------------------------------
        CHARACTER*75 FUNCTION 
     +   INFILEX(NF,FILENAME,NSCAN,NCHAN,STWV,DISP,MODE,LERR)
        IMPLICIT NONE
C Include NSMAX,NCMAX,NSCAN,NCHAN
        INCLUDE 'redlib.inc'
        INCLUDE 'futils.inc'
        INTEGER TRUEBEG,TRUELEN
        INTEGER READI
        REAL READF
C
        INTEGER NF
        INTEGER MODE
        CHARACTER*(*) FILENAME
        LOGICAL LERR
C
        INTEGER SYSTEMFUNCTION
C
        INTEGER I
        INTEGER L1,L2,L0,LCOMA
        INTEGER NMODE
        INTEGER NERR
        INTEGER NSKIP
        INTEGER NSCAN2,NCHAN2
        INTEGER ISYSTEM
        REAL STWV2,DISP2
        CHARACTER*1 COPC,CFOR,CNM,CCONT
        CHARACTER*12 CLAVE
        CHARACTER*255 CDUMMY
        LOGICAL LPROMPT,LOGFILE,LDIRECT,LMATCH,LASK
        LOGICAL LOG1,LOG2
        CHARACTER*255 OBJECT2
c------------------------------------------------------------------------------
C determinamos si hay que buscar coincidencia en los parametros fundamentales
C de la imagen que se quiere abrir
        LMATCH=.FALSE.
        NMODE=MODE
        IF(MODE.EQ.21)THEN
          NMODE=1
          LMATCH=.TRUE.
        END IF
C
        LDIRECT=.FALSE.                          !abrir el fichero directamente
        IF(NMODE.GT.10)THEN
          NMODE=MODE-10
          LDIRECT=.TRUE.
        END IF
C
        IF((NMODE.GT.4).OR.(NMODE.LT.1))THEN
          WRITE(*,*)
          WRITE(*,110)'MODE=',NMODE
          WRITE(*,101)'FATAL ERROR: invalid MODE value in '//
     +     'function INFILEX.'
          STOP
        END IF
C
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
C
        IF(LDIRECT)THEN
          IF(.NOT.LPROMPT)THEN
            WRITE(*,101)'FATAL ERROR: no file name given.'
            STOP
          END IF
          INFILEX=FILENAME(L1:L2)
          GOTO 11
        END IF
10      IF(LPROMPT)THEN
          INFILEX(1:75)=READC(FILENAME(L1:L2),'@')
        ELSE
          INFILEX(1:75)=READC('@','@')
        END IF
        LOG1=(INDEX(INFILEX,'*').NE.0)
        LOG2=(INDEX(INFILEX,'?').NE.0)
        IF(LOG1.OR.LOG2)THEN
           ISYSTEM=SYSTEMFUNCTION('ls '//INFILEX)
           GOTO 10
        END IF
11      LCOMA=INDEX(INFILEX,',') !permitimos poner una coma detras del nombre
        IF(LCOMA.GT.1)THEN
          INQUIRE(FILE=INFILEX(1:LCOMA-1),EXIST=LOGFILE)
        ELSE
          INQUIRE(FILE=INFILEX,EXIST=LOGFILE)
        END IF
        IF(.NOT.LOGFILE)THEN
          WRITE(*,101)'ERROR: This file does not exist. Try again.'
          NERR=NERR+1
          IF(NERR.GT.10) STOP 'FATAL ERROR: too many errors.'
          GOTO 10
        END IF
C
        IF(NMODE.GT.1) GOTO 40
C
        LCOMA=INDEX(INFILEX,',')
        IF(LCOMA.GT.1)THEN
          OPEN(NF,FILE=INFILEX(1:LCOMA-1),STATUS='OLD',
     +     FORM='UNFORMATTED',ERR=999)
        ELSE
          OPEN(NF,FILE=INFILEX,STATUS='OLD',FORM='UNFORMATTED',ERR=999)
        END IF
        READ(NF,ERR=20)CLAVE
        IF(CLAVE.EQ.CLAVE_RED)THEN
C..............................................................................
          READ(NF)NSCAN2,NCHAN2
          WRITE(*,110)'>>> NSCAN : ',NSCAN2
          WRITE(*,110)'>>> NCHAN : ',NCHAN2
C..............................................................................
          READ(NF)STWV2,DISP2
          WRITE(CDUMMY,*)STWV2
          L1=TRUEBEG(CDUMMY)
          L2=TRUELEN(CDUMMY)
          WRITE(*,101)'>>> STWV  : '//CDUMMY(L1:L2)
          WRITE(CDUMMY,*)DISP2
          L1=TRUEBEG(CDUMMY)
          L2=TRUELEN(CDUMMY)
          WRITE(*,101)'>>> DISP  : '//CDUMMY(L1:L2)
C..............................................................................
          READ(NF)AIRMASS
          READ(NF)TIMEXPOS
C..............................................................................
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
              WRITE(*,101)'WARNING: Expected error file has not '//
     +         'been found.'
              WRITE(*,101)'Current file name: '//
     +         INFILEX(1:TRUELEN(INFILEX))
              WRITE(*,100)'Do you want to continue with this '//
     +         'file (y/n) '
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
              WRITE(*,101)'WARNING: Unexpected error file has '//
     +         'been found.'
              WRITE(*,101)'Current file name: '//
     +         INFILEX(1:TRUELEN(INFILEX))
              WRITE(*,100)'Do you want to continue with this '//
     +         'file (y/n) '
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
C..............................................................................
          READ(NF)L0
          IF(L0.GT.0) READ(NF)FITSFILE(1:L0)
          IF(L0+1.LE.255)THEN
            DO I=L0+1,255
              FITSFILE(I:I)=CHAR(32)
            END DO
          END IF
C..............................................................................
          READ(NF)L0
          IF(L0.GT.0) READ(NF)COMMENT(1:L0)
          IF(L0+1.LE.255)THEN
            DO I=L0+1,255
              COMMENT(I:I)=CHAR(32)
            END DO
          END IF
C..............................................................................
          IF(LMATCH)THEN                 !comparamos si los valores son iguales
            IF(NSCAN.NE.NSCAN2)THEN
              WRITE(*,100)'WARNING: NSCAN does not match requested '//
     +         'value: '
              WRITE(*,*)NSCAN
              WRITE(*,100)'Continue (y/n) '
              CNM(1:1)=READC('n','yn')
              IF(CNM.EQ.'n')THEN
                CLOSE(NF)
                STOP
              END IF
              NSCAN=NSCAN2
            END IF
            IF(NCHAN.NE.NCHAN2)THEN
              WRITE(*,100)'WARNING: NCHAN does not match requested '//
     +         'value: '
              WRITE(*,*)NCHAN
              WRITE(*,100)'Continue (y/n) '
              CNM(1:1)=READC('n','yn')
              IF(CNM.EQ.'n')THEN
                CLOSE(NF)
                STOP
              END IF
              NCHAN=NCHAN2
            END IF
            IF(STWV.NE.STWV2)THEN
              WRITE(*,100)'WARNING: STWV does not match requested '//
     +         'value: '
              WRITE(*,*)STWV
              WRITE(*,100)'Continue (y/n) '
              CNM(1:1)=READC('n','yn')
              IF(CNM.EQ.'n')THEN
                CLOSE(NF)
                STOP
              END IF
              STWV=STWV2
            END IF
            IF(DISP.NE.DISP2)THEN
              WRITE(*,100)'WARNING: DISP does not match requested '//
     +         'value: '
              WRITE(*,*)DISP
              WRITE(*,100)'Continue (y/n) '
              CNM(1:1)=READC('n','yn')
              IF(CNM.EQ.'n')THEN
                CLOSE(NF)
                STOP
              END IF
              DISP=DISP2
            END IF
            IF(OBJECT2(1:L1).NE.OBJECT(1:TRUELEN(OBJECT)))THEN
              WRITE(*,101)'WARNING: OBJECT does not match requested '//
     +         'value:'
              WRITE(*,101)'>>> '//OBJECT(1:TRUELEN(OBJECT))
ccc           WRITE(*,100)'Continue (y/n) '
ccc           CNM(1:1)=READC('n','yn')
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
C..............................................................................
        END IF
C
20      CLOSE(NF)
        WRITE(*,101)'>>> WARNING <<<'
        WRITE(*,101)'This file does not contain head information.'
        WRITE(*,100)'Do you want to read it anyway (y/n) '
        COPC(1:1)=READC('y','yn')
        IF(COPC.EQ.'y')THEN
          WRITE(*,101)'(1) unformatted file'
          WRITE(*,101)'(2) formatted file (multiple column)'
          WRITE(*,100)'Option'
          CFOR(1:1)=READC('@','12')
          IF(CFOR.EQ.'1')THEN
            OPEN(NF,FILE=INFILEX,STATUS='OLD',FORM='UNFORMATTED')
            WRITE(*,100)'NSCAN'
            NSCAN=READI('@')
            WRITE(*,100)'NCHAN'
            NCHAN=READI('@')
            WRITE(*,100)'STWV'
            STWV=READF('0.0')
            WRITE(*,100)'DISP'
            DISP=READF('0.0')
            AIRMASS=0.
            TIMEXPOS=0.
            OBJECT=CHAR(32)
            FITSFILE=CHAR(32)
            COMMENT=CHAR(32)
          ELSEIF(CFOR.EQ.'2')THEN
            WRITE(*,100)'No. of lines to be skipped'
            NSKIP=READI('0')
            OPEN(NF,FILE=INFILEX,STATUS='OLD',FORM='FORMATTED')
            DO I=1,NSKIP
              READ(NF,*)
            END DO
            WRITE(*,100)'No. of columns'
            MODE=-READI('2')
          ELSE
            WRITE(*,101)'FATAL ERROR in INFILEX!!!'
            STOP
          END IF
        ELSE
          GOTO 10
        END IF
C
30      IF(NSCAN.GT.NSMAX)THEN
          WRITE(*,*)
          WRITE(*,100)'NSCAN, NSMAX: '
          WRITE(*,*) NSCAN,NSMAX
          WRITE(*,101)'FATAL ERROR: NSCAN.GT.NSMAX'
          CLOSE(NF)
          STOP
        END IF
        IF(NCHAN.GT.NCMAX)THEN
          WRITE(*,*)
          WRITE(*,100)'NCHAN, NCMAX: '
          WRITE(*,*) NCHAN,NCMAX
          WRITE(*,101)'FATAL ERROR: NCHAN.GT.NCMAX'
          CLOSE(NF)  
          STOP
        END IF
        RETURN
C
40      IF(NMODE.EQ.2)THEN
          OPEN(NF,FILE=INFILEX,STATUS='OLD',FORM='UNFORMATTED',ERR=999)
        ELSEIF(NMODE.EQ.3)THEN
          OPEN(NF,FILE=INFILEX,STATUS='OLD',FORM='FORMATTED',ERR=999)
        ELSEIF(NMODE.EQ.4)THEN
ccc          OPEN(NF,FILE=INFILEX,STATUS='OLD',READONLY,FORM='FORMATTED',
ccc     +     RECORDTYPE='FIXED',RECL=2880,ACCESS='DIRECT',ERR=999)
          OPEN(NF,FILE=INFILEX,STATUS='OLD',FORM='FORMATTED',
     +     RECL=2880,ACCESS='DIRECT',ERR=999)
        END IF
C
        RETURN
C
999     WRITE(*,101)'FATAL ERROR: opening the file: '//
     +   INFILEX(1:TRUELEN(INFILEX))
        STOP
100     FORMAT(A,$)
101     FORMAT(A)
110     FORMAT(A,I6)
        END
