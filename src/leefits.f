C------------------------------------------------------------------------------
C Version 19-July-2001                                          file: leefits.f
C @ ncl & fjg
C------------------------------------------------------------------------------
Comment
C
C Program: leefits
C Classification: input/output
C Description: Reads a FITS file and creates a new file with REDUCEME format.
C
Comment
C
C Este programa lee un fichero en formato FITS normal y lo reescribe en otro
C fichero con el formato para la reduccion con nuestros programas.
C
        PROGRAM LEEFITS
        IMPLICIT NONE
C
        INTEGER MAXNAXIS
        PARAMETER(MAXNAXIS=3)               !valor maximo admisible para NAXIS
        INTEGER NRESERVK
        PARAMETER(NRESERVK=3)                   !numero de "reserved keywords"
        INTEGER NEXTRAK
        PARAMETER(NEXTRAK=3)  !numero de "keywords" extra (no imprescindibles)
C
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER TRUELEN
        INTEGER READILIM
        REAL READF
C
        INTEGER BITPIX
        INTEGER NAXIS(0:MAXNAXIS)                               !NAXIS=NAXIS(0)
        INTEGER ICHANGE
        INTEGER I,J,I0,L0,K,LLL
        INTEGER NS1,NS2
        INTEGER IHEADER
        REAL CRPIX1,CRVAL1,CDELT1
        REAL FNULL
        REAL IMAGEN(NCMAX,NSMAX)
        REAL DATAMIN,DATAMAX
        CHARACTER*1 CHIST,CWHOLE
        CHARACTER*8 RESERVK(NRESERVK),EXTRAK(NEXTRAK)
        CHARACTER*8 RESERVK_TEMPORAL
        CHARACTER*50 CDUMMY
        CHARACTER*80 INFILE,OUTFILE
        CHARACTER*2880 HDU
        LOGICAL LHIST,LEXTENSION
        LOGICAL LRESERVK(NRESERVK),LEXTRAK(NEXTRAK)
        LOGICAL LRUN,LMANUAL,LHTML
        LOGICAL LNULL(NCMAX,NSMAX)
        LOGICAL ANYNULL
        LOGICAL LOOP
C
        COMMON/BLKIMAGEN/IMAGEN !imagen FITS leida en formato REAL
        COMMON/BLKLNULL/LNULL,ANYNULL   !mascara que indica si existen NaN, etc.
        COMMON/BLKDATMINMAX/DATAMIN,DATAMAX !maximo y minimo en la imagen
        COMMON/BLKNAXIS/NAXIS   !dimensiones
        COMMON/BLKFITSFILE/INFILE !nombre de la imagen
C------------------------------------------------------------------------------
C Introducimos en el DATA las "reserved keywords" que seran buscadas en el
C formato FITS normal. NO CAMBIAR el orden de las palabras; el programa
C asume el orden establecido. Si se deseen buscar nuevas "keywords" basta
C con introducirlas al final del DATA. 
C IMPORTANTE: rellenar las "keywords" con blancos si su longitud es inferior
C a 8 caracteres.
        DATA (RESERVK(I),I=1,NRESERVK)/
     +   'OBJECT  ','EXPTIME ','AIRMASS '/
        DATA (EXTRAK(I),I=1,NEXTRAK)/
     +   'CRPIX1  ','CRVAL1  ','CDELT1  '/
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        THISPROGRAM='leefits'
        CALL WELCOME('9-April-1999')
C------------------------------------------------------------------------------
5       WRITE(*,101) '* The following reserved Keywords will be '//
     +   'required: '
        DO I=1,NRESERVK
          WRITE(*,'(I1,1X,A1,1X,A)') I,'-',RESERVK(I)
        END DO
        WRITE(*,100) 'Keyword to be changed (0=NONE) '
        ICHANGE=READILIM('0',0,NRESERVK)
        IF(ICHANGE.NE.0)THEN
          WRITE(*,'(A)') 'The following keyword will be changed: '//
     +     RESERVK(ICHANGE)
          WRITE(*,100) 'New keyword (max. 8 characters)'
          RESERVK_TEMPORAL(1:8)=READC('@','@')
          RESERVK(ICHANGE)=RESERVK_TEMPORAL
          GOTO 5
        END IF
C
        DO I=1,NRESERVK
          LRESERVK(I)=.FALSE.  !inicializamos: ninguna "reserved keyword found" 
        END DO
        DO I=1,NEXTRAK
          LEXTRAK(I)=.FALSE.                !ninguna extra "keyword" encontrada
        END DO
C
        WRITE(*,*)
        WRITE(*,100) 'FITS file name'
        INFILE=INFILEX(20,'@',0,0,0.,0.,4,.FALSE.)
        FITSFILE=INFILE
C
        WRITE(*,100) 'Run silently (y/n) '
        CHIST(1:1)=READC('y','yn')
        LHIST=(CHIST.EQ.'n')
        WRITE(*,*)
C
C Leemos el primer HDU (Header Data Unit)
        IHEADER=1
        READ(20,'(A2880)',REC=IHEADER) HDU
        WRITE(*,101) '---> REQUIRED KEYWORDS <---'
C
C Comprobamos el tipo de formato buscando las "required keywords"
C----------------------->123456789012345678901234567890
        IF(HDU(1:30).NE.'SIMPLE  =                    T')THEN
          WRITE(*,101) 'FATAL ERROR: the file does not conform to '//
     +     'FITS standards'
          CLOSE(20)
          STOP
        ELSE
          WRITE(*,101) '> File seems to conform FITS standards'
        END IF
C
        NAXIS(0)=0
        LEXTENSION=.FALSE.      !suponemos que inicialmente no es una EXTENSION
10      READ(HDU(91:110),*)BITPIX
        READ(HDU(171:190),*)NAXIS(0)
        WRITE(*,110) '> BITPIX  = ',BITPIX
        WRITE(*,110) '> NAXIS   = ',NAXIS(0)
        IF(NAXIS(0).GT.MAXNAXIS)THEN
          WRITE(*,101) 'FATAL ERROR: NAXIS.GT.MAXNAXIS'
          CLOSE(20)
          STOP
        END IF
C
C Si NAXIS=0 entonces buscamos la confirmacion de "standard extension".
C En ese caso, buscamos el siguiente HDU, y verificamos que la palabra
C clave XTENSION esta presente.
        IF(NAXIS(0).EQ.0)THEN
C---------------------------->123456789012345678901234567890
          IF(HDU(241:270).NE.'EXTEND  =                    T')THEN
            WRITE(*,101)
     +       '> FATAL ERROR: Standard extension not confirmed.'
            CLOSE(20)
            STOP
          END IF
          WRITE(*,101) '> There are standard extensions (?)'
          IHEADER=IHEADER+1
          READ(20,'(A2880)',REC=IHEADER) HDU
          IF(HDU(1:9).NE.'XTENSION=')THEN
            WRITE(*,101) 'FATAL ERROR: XTENSION has not been found.'
            CLOSE(20)
            STOP
          END IF
          WRITE(*,101) '> '//HDU(1:30)
          IF(HDU(12:19).NE.'BINTABLE')THEN
            WRITE(*,101)
     +       'FATAL ERROR: This extension cannot be handled.'
            CLOSE(20)
            STOP
          END IF
          LEXTENSION=.TRUE.
          GOTO 10
        END IF
C
C leemos NAXIS1, NAXIS2,..., NAXISn (con n=NAXIS)
        DO I=1,MAXNAXIS
          NAXIS(I)=0
        END DO
        DO I=1,NAXIS(0)
          READ(HDU((I+2)*80+11:(I+2)*80+30),*)NAXIS(I)
          WRITE(*,'(A,I1,A,I6)') '> NAXIS',I,'  = ',NAXIS(I)
        END DO
C
C Si es una EXTENSION, tratamos el fichero de forma diferente
        IF(LEXTENSION) GOTO 200
C
        L0=NAXIS(0)+4             !numero de linea (de 80 char.) del HDU
        WRITE(*,101) '---> RESERVED KEYWORDS <---'
C Buscamos "reserved keywords"
20      DO I=L0,36
          I0=(I-1)*80
          IF(HDU(I0+1:I0+3).EQ.'END') GOTO 22                     !fin de HDU's
          DO K=1,NEXTRAK
            IF(.NOT.LEXTRAK(K))THEN
              IF(HDU(I0+1:I0+8).EQ.EXTRAK(K))THEN
                LEXTRAK(K)=.TRUE.
                WRITE(*,101) '> '//HDU(I0+1:I0+30)
                IF(K.EQ.1)THEN
                  IF(TRUELEN(HDU(I0+11:I0+30)).GT.0)THEN
                    READ(HDU(I0+11:I0+30),*) CRPIX1
                  ELSE
                    WRITE(*,101) 'ERROR: '//RESERVK(K)//
     +               ' --> invalid field'
                    LRESERVK(K)=.FALSE.
                  END IF
                ELSEIF(K.EQ.2)THEN
                  IF(TRUELEN(HDU(I0+11:I0+30)).GT.0)THEN
                    READ(HDU(I0+11:I0+30),*) CRVAL1
                  ELSE
                    WRITE(*,101) 'ERROR: '//RESERVK(K)//
     +               ' --> invalid field'
                    LRESERVK(K)=.FALSE.
                  END IF
                ELSEIF(K.EQ.3)THEN
                  IF(TRUELEN(HDU(I0+11:I0+30)).GT.0)THEN
                    READ(HDU(I0+11:I0+30),*) CDELT1
                  ELSE
                    WRITE(*,101) 'ERROR: '//RESERVK(K)//
     +               ' --> invalid field'
                    LRESERVK(K)=.FALSE.
                  END IF
                ELSE
                  WRITE(*,101)
     +             'FATAL ERROR: Invalid Required Keyword number.'
                  CLOSE(20)
                  STOP
                END IF
                GOTO 21                         !provocamos la salida del bucle
              END IF
            END IF
          END DO
          DO K=1,NRESERVK                     !buscamos cada "reserved keyword"
            IF(.NOT.LRESERVK(K))THEN             !solo leemos primera aparicion 
              IF(HDU(I0+1:I0+8).EQ.RESERVK(K))THEN
                LRESERVK(K)=.TRUE.
                WRITE(*,101) '> '//HDU(I0+1:I0+30)
                IF(K.EQ.1)THEN
                  LLL=INDEX(HDU(I0+12:I0+30),CHAR(39))-1
                  OBJECT=HDU(I0+12:I0+11+LLL)
                ELSEIF(K.EQ.2)THEN
                  IF(TRUELEN(HDU(I0+11:I0+30)).GT.0)THEN
                    READ(HDU(I0+11:I0+30),*) TIMEXPOS
                  ELSE
                    WRITE(*,101) 'ERROR: '//RESERVK(K)//
     +               ' --> invalid field'
                    LRESERVK(K)=.FALSE.
                  END IF
                ELSEIF(K.EQ.3)THEN
                  IF(TRUELEN(HDU(I0+11:I0+30)).GT.0)THEN
                    READ(HDU(I0+11:I0+30),*) AIRMASS
                  ELSE
                    WRITE(*,101) 'ERROR: '//RESERVK(K)//
     +               ' --> invalid field'
                    LRESERVK(K)=.FALSE.
                  END IF
                ELSE
                  WRITE(*,101)
     +             'FATAL ERROR: Invalid Required Keyword number.'
                  CLOSE(20)
                  STOP
                END IF
                GOTO 21                         !provocamos la salida del bucle
              END IF
            END IF
          END DO
21        CONTINUE
        END DO
        L0=1
        IHEADER=IHEADER+1
        READ(20,'(A2880)',REC=IHEADER) HDU
        GOTO 20
C
22      WRITE(*,110) '> No. of HDUs read: ',IHEADER
        DO I=1,NRESERVK
          IF(.NOT.LRESERVK(I))THEN
            WRITE(*,101) 'WARNING: '//RESERVK(I)//' has not been found.'
            WRITE(*,100) 'New value for '//RESERVK(I)//' '
            IF(I.EQ.1)THEN
              WRITE(*,100) '(CHARACTER*20) '
              OBJECT=READC(INFILE,'@')
            ELSEIF(I.EQ.2)THEN
              WRITE(*,100) '(REAL) '
              TIMEXPOS=READF('-999')
            ELSEIF(I.EQ.3)THEN
              WRITE(*,100) '(REAL) '
              AIRMASS=READF('-999')
            END IF
            LRESERVK(I)=.TRUE.
          END IF
        END DO
C
C Si NAXIS > 3 el programa hay que modificarlo
        IF(NAXIS(0).EQ.1)THEN
          NAXIS(2)=1
          NAXIS(3)=1
        ELSEIF(NAXIS(0).EQ.2)THEN
          NAXIS(3)=1
        ELSEIF(NAXIS(0).GT.3)THEN
          WRITE(*,101) 'FATAL ERROR: NAXIS.GT.3'
          WRITE(*,101) 
     +     'This situation has not been implemented. Sorry.'
          CLOSE(20)
          STOP
        END IF
C
        WRITE(*,*)
        WRITE(*,101) 'Comments (max. 255 characters, '//
     +   '<ENTER>=No-comment)? '
        READ(*,'(A)')COMMENT
        CALL LRUNX(LRUN,LMANUAL,LHTML)
        IF(LRUN)THEN
          WRITE(*,101)COMMENT(1:TRUELEN(COMMENT))
        END IF
        IF(LMANUAL)THEN
          WRITE(*,101)COMMENT(1:TRUELEN(COMMENT))
          WRITE(*,101) '\\ttshade{'//COMMENT(1:TRUELEN(COMMENT))//'}'
        END IF
        IF(LHTML)THEN
          WRITE(*,101) '<FONT COLOR="#FF0000">'//
     +     COMMENT(1:TRUELEN(COMMENT))//'</FONT>'
        END IF
C
        IF(LEXTRAK(1).AND.LEXTRAK(2).AND.LEXTRAK(3))THEN
          IF(CRPIX1.EQ.1.0)THEN
            STWV=CRVAL1
            DISP=CDELT1
          ELSE
            STWV=CRVAL1+CDELT1*(1.-CRPIX1)
            DISP=CDELT1
          END IF
        ELSEIF(LEXTRAK(2).AND.LEXTRAK(3))THEN
          STWV=CRVAL1
          DISP=CDELT1
        ELSE
          STWV=0.
          DISP=0.
        END IF
        WRITE(*,*)
        WRITE(*,100) '> STWV: '
        WRITE(*,*) STWV
        WRITE(*,100) '> DISP: '
        WRITE(*,*) DISP
        WRITE(*,*)
C
        IF(NAXIS(2).GT.1)THEN
          WRITE(*,100) 'Save whole frame (y/n) '
          CWHOLE(1:1)=READC('y','yn')
          IF(CWHOLE.EQ.'n')THEN
            WRITE(*,100) 'First scan: '
            NS1=READILIM('1',1,NAXIS(2))
            WRITE(*,100) 'Last  scan: '
            WRITE(CDUMMY,*) NAXIS(2)
            NS2=READILIM(CDUMMY,NS1,NAXIS(2))
          ELSE
            NS1=1
            NS2=NAXIS(2)
          END IF
        ELSE
          NS1=1
          NS2=1
        END IF
        NSCAN=NS2-NS1+1
        CLOSE(20) !cerramos el fichero: SLEEFITS lo vuelve a abrir y cerrar
C
        LOOP=.TRUE.
        K=0
        DO WHILE(LOOP)
          IF(NAXIS(3).EQ.1)THEN
            LOOP=.FALSE.
            K=1
          ELSE
            WRITE(*,100) '>>> NAXIS3='
            WRITE(*,*) NAXIS(3)
            WRITE(*,100) 'Which NAXIS3 do you want to save (0=QUIT)'
            K=K+1
            IF(K.GT.NAXIS(3)) K=0
            WRITE(CDUMMY,*) K
            K=READILIM(CDUMMY,0,NAXIS(3))
            IF(K.EQ.0) LOOP=.FALSE.
          END IF
          IF(K.GT.0)THEN
            WRITE(*,100) 'Output file name'
            OUTFILE=OUTFILEX(30,'@',NSCAN,NAXIS(1),STWV,DISP,1,.FALSE.)
C leemos la matriz de datos
            CALL SLEEFITS(LHIST,K)
            WRITE(*,100) '>>> DATAMIN: '
            WRITE(*,*) DATAMIN
            WRITE(*,100) '>>> DATAMAX: '
            WRITE(*,*) DATAMAX
            IF(ANYNULL)THEN
              WRITE(*,*)
              WRITE(*,101) '****************************************'
              WRITE(*,101) 'WARNING: image contains undefined pixels'
              WRITE(*,101) '****************************************'
              WRITE(*,*)
              WRITE(*,100) 'Real value for undefined pixels '
              FNULL=READF('@')
              DO I=NS1,NS2
                DO J=1,NAXIS(1)
                  IF(LNULL(J,I)) IMAGEN(J,I)=FNULL
                END DO
              END DO
            END IF
            WRITE(*,100) 'Please wait (saving file)...'
            DO I=NS1,NS2
              WRITE(30) (IMAGEN(J,I),J=1,NAXIS(1))
            END DO
            WRITE(*,101) '  ..OK!'
            CLOSE(30)
          END IF
        END DO
        STOP
C------------------------------------------------------------------------------
200     CONTINUE
        WRITE(*,101) 'This program is unable to reed an extension '//
     +   'Binary Table'
220     DO I=1,36
          I0=(I-1)*80
          IF(HDU(I0+1:I0+3).EQ.'END') GOTO 222                    !fin de HDU's
        END DO
        IHEADER=IHEADER+1
        READ(20,'(A2880)',REC=IHEADER) HDU
        GOTO 220
222     WRITE(*,110) '> No. of HDUs read: ',IHEADER
C
ccc230     CLOSE(20)
        CLOSE(20)
C
        STOP
C
100     FORMAT(A,$)
101     FORMAT(A)
110     FORMAT(A,I6)
        END
C
C******************************************************************************
C Subrutina para leer una image FITS (si LHIST=.TRUE., muestra la informacion
C de la cabecera FITS).
        SUBROUTINE SLEEFITS(LHIST,NAXIS3)
        IMPLICIT NONE
        LOGICAL LHIST
        INTEGER NAXIS3
C
        INCLUDE 'redlib.inc'
        INTEGER TRUELEN
        INTEGER READI
C
        INTEGER MAXNAXIS
        PARAMETER(MAXNAXIS=3)                !valor maximo admisible para NAXIS
C
        INTEGER JROW(NCMAX)
        INTEGER I,J,L
        INTEGER FIRSTPIX
        INTEGER BITPIX,NAXIS(0:MAXNAXIS)
        INTEGER ISTATUS,IREADWRITE,IUNIT
        INTEGER BLOCKSIZE,NULLVAL
        INTEGER NKEYS,NSPACE,NFOUND
        REAL IMAGEN(NCMAX,NSMAX),FROW(NCMAX)
        REAL DATAMIN,DATAMAX
        DOUBLE PRECISION DROW(NCMAX)
        CHARACTER*50 MYCOMMENT,CDUMMY
        CHARACTER*80 INFILE,CLINEA
        LOGICAL ANYNULL,LANYNULL
        LOGICAL LROW(NCMAX),LNULL(NCMAX,NSMAX)
C
        COMMON/BLKIMAGEN/IMAGEN !imagen FITS leida en formato REAL
        COMMON/BLKLNULL/LNULL,ANYNULL   !mascara que indica si existen NaN, etc.
        COMMON/BLKDATMINMAX/DATAMIN,DATAMAX !maximo y minimo en la imagen
        COMMON/BLKNAXIS/NAXIS   !dimensiones
        COMMON/BLKFITSFILE/INFILE !nombre de la imagen
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C inicializamos variables
        ISTATUS=0               !controla posibles errores durante la ejecucion
        IREADWRITE=0                      !la imagen se abrira en modo READONLY
        NULLVAL=-999
        LANYNULL=.FALSE.
C------------------------------------------------------------------------------
C localizamos un numero de unidad de fichero no utilizada
ccc     CALL FTGIOU(IUNIT,ISTATUS)
        IUNIT=80 !ojo, IUNIT=99 entra en conflicto con el Postcript y es
                 !precisamente este numero el que toma esta funcion
C abrimos el fichero
        CALL FTOPEN(IUNIT,INFILE,IREADWRITE,BLOCKSIZE,ISTATUS)
C determinamos el numero de keywords en la cabecera y las mostramos
        IF(LHIST)THEN
          CALL FTGHSP(IUNIT,NKEYS,NSPACE,ISTATUS)
          DO I=1,NKEYS
            CALL FTGREC(IUNIT,I,CLINEA,ISTATUS)
            L=TRUELEN(CLINEA)
            WRITE(*,101)CLINEA(1:L)
          END DO
          IF(ISTATUS.EQ.0)THEN                                !todo ha ido bien
            WRITE(*,101)'END'
            WRITE(*,*)
          END IF
        END IF
C leemos BITPIX
        CALL FTGKYJ(IUNIT,'BITPIX',BITPIX,MYCOMMENT,ISTATUS)
        WRITE(CDUMMY,*) BITPIX
        WRITE(*,100) 'BITPIX (double check) '
        BITPIX=READI(CDUMMY)
C comprobamos que NAXIS=2
        CALL FTGKYJ(IUNIT,'NAXIS',NAXIS(0),MYCOMMENT,ISTATUS)
        IF(NAXIS(0).GT.MAXNAXIS)THEN
          WRITE(*,100)'FATAL ERROR: NAXIS >'
          WRITE(*,*)MAXNAXIS
          CALL FTCLOS(IUNIT,ISTATUS)
          STOP
        END IF
C leemos NAXIS1...NAXIS3 [notar que el quinto parametro es NAXIS(1) en lugar
C de NAXIS para asi recuperar NAXIS(1)...NAXIS(3)]
        CALL FTGKNJ(IUNIT,'NAXIS',1,NAXIS(0),NAXIS(1),NFOUND,ISTATUS)
        IF(NAXIS(1).GT.NCMAX)THEN
          WRITE(*,101)'* FATAL ERROR in subroutine SLEEFITS:'
          WRITE(*,101)'NAXIS(1) > NCMAX'
          STOP
        END IF
        IF(NAXIS(2).GT.NSMAX)THEN
          WRITE(*,101)'* FATAL ERROR in subroutine SLEEFITS:'
          WRITE(*,101)'NAXIS(2) > NSMAX'
          STOP
        END IF
C leemos la imagen
        WRITE(*,100) 'Please wait (reading FITS file)...'
        IF(BITPIX.EQ.16)THEN
          DO I=1,NAXIS(2)
            FIRSTPIX=(I-1)*NAXIS(1)+1+(NAXIS3-1)*(NAXIS(1)*NAXIS(2))
            CALL FTGPFJ(IUNIT,1,FIRSTPIX,NAXIS(1),JROW(1),LROW(1),
     +       ANYNULL,ISTATUS)
            DO J=1,NAXIS(1)
              IMAGEN(J,I)=REAL(JROW(J))
            END DO
            IF(ANYNULL)THEN
              DO J=1,NAXIS(1)
                LNULL(J,I)=LROW(J)
              END DO
              LANYNULL=.TRUE.
            END IF
          END DO
        ELSEIF(BITPIX.EQ.32)THEN
          DO I=1,NAXIS(2)
            FIRSTPIX=(I-1)*NAXIS(1)+1+(NAXIS3-1)*(NAXIS(1)*NAXIS(2))
            CALL FTGPFE(IUNIT,1,FIRSTPIX,NAXIS(1),FROW(1),LROW(1),
     +       ANYNULL,ISTATUS)
            DO J=1,NAXIS(1)
              IMAGEN(J,I)=FROW(J)
            END DO
            IF(ANYNULL)THEN
              DO J=1,NAXIS(1)
                LNULL(J,I)=LROW(J)
              END DO
              LANYNULL=.TRUE.
            END IF
          END DO
ccc       CALL FTG2DE(IUNIT,1,NULLVAL,NCMAX,NAXIS(1),NAXIS(2),
ccc     +     IMAGEN,ANYNULL,ISTATUS)
        ELSEIF(BITPIX.EQ.-32)THEN
          DO I=1,NAXIS(2)
            FIRSTPIX=(I-1)*NAXIS(1)+1+(NAXIS3-1)*(NAXIS(1)*NAXIS(2))
            CALL FTGPFE(IUNIT,1,FIRSTPIX,NAXIS(1),FROW(1),LROW(1),
     +       ANYNULL,ISTATUS)
            DO J=1,NAXIS(1)
              IMAGEN(J,I)=FROW(J)
            END DO
            IF(ANYNULL)THEN
              DO J=1,NAXIS(1)
                LNULL(J,I)=LROW(J)
              END DO
              LANYNULL=.TRUE.
            END IF
          END DO
ccc       CALL FTG2DE(IUNIT,1,NULLVAL,NCMAX,NAXIS(1),NAXIS(2),
ccc     +     IMAGEN,ANYNULL,ISTATUS)
        ELSEIF(BITPIX.EQ.-64)THEN
          DO I=1,NAXIS(2)
            FIRSTPIX=(I-1)*NAXIS(1)+1+(NAXIS3-1)*(NAXIS(1)*NAXIS(2))
            CALL FTGPFD(IUNIT,1,FIRSTPIX,NAXIS(1),DROW(1),LROW(1),
     +       ANYNULL,ISTATUS)
            DO J=1,NAXIS(1)
              IMAGEN(J,I)=REAL(DROW(J))
            END DO
            IF(ANYNULL)THEN
              DO J=1,NAXIS(1)
                LNULL(J,I)=LROW(J)
              END DO
              LANYNULL=.TRUE.
            END IF
          END DO
        ELSE
          WRITE(*,100)'FATAL ERROR in subroutine SLEEFITS: BITPIX ='
          WRITE(*,*) BITPIX
          CALL FTCLOS(IUNIT,ISTATUS)
          STOP
        END IF
        WRITE(*,101) '  ..OK!'
C cerramos el fichero
        CALL FTCLOS(IUNIT,ISTATUS)
C liberamos el numero de unidad del fichero utilizado
        CALL FTFIOU(IUNIT,ISTATUS)
C calculamos maximo y minimo
        DATAMIN=1.E30
        DATAMAX=-1.E30
        DO I=1,NAXIS(2)
          DO J=1,NAXIS(1)
            IF(.NOT.LNULL(J,I))THEN
              IF(IMAGEN(J,I).LT.DATAMIN) DATAMIN=IMAGEN(J,I)
              IF(IMAGEN(J,I).GT.DATAMAX) DATAMAX=IMAGEN(J,I)
            END IF
          END DO
        END DO
C chequeamos si se ha producido algun error
        IF(ISTATUS.GT.0)THEN
          CALL PRINTERROR(ISTATUS)
        END IF
        ANYNULL=LANYNULL                  !basta que haya ocurrido una sola vez
C------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
        END
C
C******************************************************************************
C
        SUBROUTINE PRINTERROR(ISTATUS)
C Print out the FITSIO error messages to the user
        INTEGER ISTATUS
        CHARACTER ERRTEXT*30,ERRMESSAGE*80
C Check if status is OK (no error); if so, simply return
        IF(ISTATUS.LE.0) RETURN
C Get the text string which describes the error
        CALL FTGERR(ISTATUS,ERRTEXT)
        WRITE(*,'(A,$)')'FITSIO Error Status = '
        WRITE(*,*)ISTATUS
        WRITE(*,'(A)')ERRTEXT
C Read and print out all the error messages on the FITSIO stack
        CALL FTGMSG(ERRMESSAGE)
        DO WHILE(ERRMESSAGE.NE.' ')
          WRITE(*,'(A)') ERRMESSAGE
          CALL FTGMSG(ERRMESSAGE)
        END DO
        END
