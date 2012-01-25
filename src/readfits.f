C------------------------------------------------------------------------------
C Version 24-March-1998                                        file: readfits.f
C @ ncl & fjg
C------------------------------------------------------------------------------
Comment
C
C Program: readfits
C Classification: input/output
C Description: Reads a FITS file and creates a new file with REDUCEME format
C (this program is obsolete; use leefits).
C
Comment
C
C Este programa lee un fichero en formato FITS normal y lo reescribe en otro
C fichero con el formato para la reduccion con nuestros programas.
C
        PROGRAM READFITS
        IMPLICIT NONE
C
        INTEGER MAXNAXIS
        PARAMETER(MAXNAXIS=9)               !valor maximo admisible para NAXIS
        INTEGER NRESERVK
        PARAMETER(NRESERVK=5)                   !numero de "reserved keywords"
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
        REAL BSCALE
        REAL BZERO
C
        INTEGER I,J,I0,L0,K,IX,IY,LLL
        INTEGER I1,I2
        INTEGER BYTEPIX,NIBL,IBL,NBREAD,N
        INTEGER NS1,NS2,NS0
        INTEGER IHEADER
        REAL ROW(NCMAX)
        REAL CRPIX1,CRVAL1,CDELT1
        CHARACTER*1 CHIST,CSWAP,CWHOLE,CWARNING
        CHARACTER*8 RESERVK(NRESERVK),EXTRAK(NEXTRAK)
        CHARACTER*8 RESERVK_TEMPORAL
        CHARACTER*50 CDUMMY
        CHARACTER*75 INFILE,OUTFILE
        CHARACTER*2880 HDU
        LOGICAL LSWAP
        LOGICAL LHIST,LEXTENSION
        LOGICAL LRESERVK(NRESERVK),LEXTRAK(NEXTRAK)
        LOGICAL LRUN,LMANUAL,LHTML
C------------------------------------------------------------------------------
C Introducimos en el DATA las "reserved keywords" que seran buscadas en el
C formato FITS normal. NO CAMBIAR el orden de las palabras; el programa
C asume el orden establecido. Si se deseen buscar nuevas "keywords" basta
C con introducirlas al final del DATA. 
C IMPORTANTE: rellenar las "keywords" con blancos si su longitud es inferior
C a 8 caracteres.
        DATA (RESERVK(I),I=1,NRESERVK)/
     +   'BSCALE  ','BZERO   ','OBJECT  ','EXPTIME ','AIRMASS '/
        DATA (EXTRAK(I),I=1,NEXTRAK)/
     +   'CRPIX1  ','CRVAL1  ','CDELT1  '/
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        THISPROGRAM='readfits'
        CALL WELCOME('24-March-1998')
C------------------------------------------------------------------------------
5       WRITE(*,101)'* The following reserved Keywords will be '//
     +   'required: '
        DO I=1,NRESERVK
          WRITE(*,'(I1,1X,A1,1X,A)')I,'-',RESERVK(I)
        END DO
        WRITE(*,100)'Keyword to be changed (0=NONE) '
        ICHANGE=READILIM('0',0,NRESERVK)
        IF(ICHANGE.NE.0)THEN
          WRITE(*,'(A)')'The following keyword will be changed: '//
     +     RESERVK(ICHANGE)
          WRITE(*,100)'New keyword (max. 8 characters)'
          RESERVK_TEMPORAL(1:1)=READC('@','@')
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
        WRITE(*,100)'FITS file name'
        INFILE=INFILEX(20,'@',0,0,0.,0.,4,.FALSE.)
        FITSFILE=INFILE
C
        WRITE(*,100)'Run silently (y/n) '
        CHIST(1:1)=READC('n','yn')
        LHIST=(CHIST.EQ.'n')
        WRITE(*,*)
C
C Leemos el primer HDU (Header Data Unit)
        IHEADER=1
        READ(20,'(A2880)',REC=IHEADER) HDU
        IF(LHIST) WRITE(*,101)'---> REQUIRED KEYWORDS <---'
C
C Comprobamos el tipo de formato buscando las "required keywords"
C----------------------->123456789012345678901234567890
        IF(HDU(1:30).NE.'SIMPLE  =                    T')THEN
          WRITE(*,101)'FATAL ERROR: the file does not conform to FITS'//
     +     ' standards'
          CLOSE(20)
          STOP
        ELSE
          IF(LHIST) WRITE(*,101)'> File seems to conform FITS standards'
        END IF
C
        NAXIS(0)=0
        LEXTENSION=.FALSE.      !suponemos que inicialmente no es una EXTENSION
10      READ(HDU(91:110),*)BITPIX
        READ(HDU(171:190),*)NAXIS(0)
        IF(LHIST)THEN
          WRITE(*,110)'> BITPIX  = ',BITPIX
          WRITE(*,110)'> NAXIS   = ',NAXIS(0)
        END IF
        IF(NAXIS(0).GT.MAXNAXIS)THEN
          WRITE(*,101)'FATAL ERROR: NAXIS.GT.MAXNAXIS'
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
            WRITE(*,101)'> FATAL ERROR: Standard extension'//
     +       ' not confirmed.'
            CLOSE(20)
            STOP
          END IF
          IF(LHIST) WRITE(*,101)'> There are standard extensions (?)'
          IHEADER=IHEADER+1
          READ(20,'(A2880)',REC=IHEADER) HDU
          IF(HDU(1:9).NE.'XTENSION=')THEN
            WRITE(*,101)'FATAL ERROR: XTENSION has not been found.'
            CLOSE(20)
            STOP
          END IF
          IF(LHIST) WRITE(*,101) '> '//HDU(1:30)
          IF(HDU(12:19).NE.'BINTABLE')THEN
            WRITE(*,101)'FATAL ERROR: This extension cannot be handled.'
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
          IF(LHIST) WRITE(*,'(A,I1,A,I6)')'> NAXIS',I,'  = ',NAXIS(I)
        END DO
        IF(NAXIS(0).EQ.1)THEN
          NAXIS(2)=1
        END IF
C
C Si es una EXTENSION, tratamos el fichero de forma diferente
        IF(LEXTENSION) GOTO 200
C
        L0=NAXIS(0)+4             !numero de linea (de 80 char.) del HDU
        IF(LHIST) WRITE(*,101)'---> RESERVED KEYWORDS <---'
C Buscamos "reserved keywords"
20      DO I=L0,36
          I0=(I-1)*80
          IF(HDU(I0+1:I0+3).EQ.'END') GOTO 22                     !fin de HDU's
          DO K=1,NEXTRAK
            IF(.NOT.LEXTRAK(K))THEN
              IF(HDU(I0+1:I0+8).EQ.EXTRAK(K))THEN
                LEXTRAK(K)=.TRUE.
                WRITE(*,101)'> '//HDU(I0+1:I0+30)
                IF(K.EQ.1)THEN
                  IF(TRUELEN(HDU(I0+11:I0+30)).GT.0)THEN
                    READ(HDU(I0+11:I0+30),*) CRPIX1
                  ELSE
                    WRITE(*,101)'ERROR: '//RESERVK(K)//
     +               ' --> invalid field'
                    LRESERVK(K)=.FALSE.
                  END IF
                ELSEIF(K.EQ.2)THEN
                  IF(TRUELEN(HDU(I0+11:I0+30)).GT.0)THEN
                    READ(HDU(I0+11:I0+30),*) CRVAL1
                  ELSE
                    WRITE(*,101)'ERROR: '//RESERVK(K)//
     +               ' --> invalid field'
                    LRESERVK(K)=.FALSE.
                  END IF
                ELSEIF(K.EQ.3)THEN
                  IF(TRUELEN(HDU(I0+11:I0+30)).GT.0)THEN
                    READ(HDU(I0+11:I0+30),*) CDELT1
                  ELSE
                    WRITE(*,101)'ERROR: '//RESERVK(K)//
     +               ' --> invalid field'
                    LRESERVK(K)=.FALSE.
                  END IF
                ELSE
                  WRITE(*,101)'FATAL ERROR: Invalid Required '//
     +             'Keyword number.'
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
                WRITE(*,101)'> '//HDU(I0+1:I0+30)
                IF(K.EQ.1)THEN
                  IF(TRUELEN(HDU(I0+11:I0+30)).GT.0)THEN
                    READ(HDU(I0+11:I0+30),*) BSCALE
                  ELSE
                    WRITE(*,101)'ERROR: '//RESERVK(K)//
     +               ' --> invalid field'
                    LRESERVK(K)=.FALSE.
                  END IF
                ELSEIF(K.EQ.2)THEN
                  IF(TRUELEN(HDU(I0+11:I0+30)).GT.0)THEN
                    READ(HDU(I0+11:I0+30),*) BZERO
                  ELSE
                    WRITE(*,101)'ERROR: '//RESERVK(K)//
     +               ' --> invalid field'
                    LRESERVK(K)=.FALSE.
                  END IF
                ELSEIF(K.EQ.3)THEN
                  LLL=INDEX(HDU(I0+12:I0+30),CHAR(39))-1
                  OBJECT=HDU(I0+12:I0+11+LLL)
                ELSEIF(K.EQ.4)THEN
                  IF(TRUELEN(HDU(I0+11:I0+30)).GT.0)THEN
                    READ(HDU(I0+11:I0+30),*) TIMEXPOS
                  ELSE
                    WRITE(*,101)'ERROR: '//RESERVK(K)//
     +               ' --> invalid field'
                    LRESERVK(K)=.FALSE.
                  END IF
                ELSEIF(K.EQ.5)THEN
                  IF(TRUELEN(HDU(I0+11:I0+30)).GT.0)THEN
                    READ(HDU(I0+11:I0+30),*) AIRMASS
                  ELSE
                    WRITE(*,101)'ERROR: '//RESERVK(K)//
     +               ' --> invalid field'
                    LRESERVK(K)=.FALSE.
                  END IF
                ELSE
                  WRITE(*,101) 'FATAL ERROR: Invalid Required'//
     +             ' Keyword number.'
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
22      WRITE(*,110)'> No. of HDUs read: ',IHEADER
        DO I=1,NRESERVK
          IF(.NOT.LRESERVK(I))THEN
            WRITE(*,101)'WARNING: '//RESERVK(I)//' has not been found.'
            WRITE(*,100)'New value for '//RESERVK(I)//' '
            IF(I.EQ.1)THEN
              WRITE(*,100)'(REAL) '
              BSCALE=READF('1.0')
            ELSEIF(I.EQ.2)THEN
              WRITE(*,100)'(REAL) '
              BZERO=READF('0.0')
            ELSEIF(I.EQ.3)THEN
              WRITE(*,100)'(CHARACTER*20) '
              OBJECT=READC(INFILE,'@')
            ELSEIF(I.EQ.4)THEN
              WRITE(*,100)'(REAL) '
              TIMEXPOS=READF('-999')
            ELSEIF(I.EQ.5)THEN
              WRITE(*,100)'(REAL) '
              AIRMASS=READF('-999')
            END IF
            LRESERVK(I)=.TRUE.
          END IF
        END DO
C
C Si NAXIS > 2 el programa hay que modificarlo
        IF(NAXIS(0).GT.2)THEN
          CWARNING='n'
          IF((NAXIS(0).EQ.3).AND.(NAXIS(3).EQ.1))THEN
            WRITE(*,101)'WARNING: NAXIS.GT.2 but NAXIS(3)=1'
            WRITE(*,100)'Do you want to continue (y/n) [y] '
            CWARNING(1:1)=READC('y','yn')
          END IF
          IF(CWARNING.EQ.'n')THEN
            WRITE(*,101)'FATAL ERROR: NAXIS.GT.2'
            WRITE(*,101) 'This situation has not been implemented.'//
     +       ' Sorry.'
            CLOSE(20)
            STOP
          END IF
        END IF
C
C Si BITPIX no es 16, el programa se detiene
        IF(BITPIX.NE.16)THEN
          WRITE(*,101)'FATAL ERROR: BITPIX value cannot be handled.'
          CLOSE(20)
          STOP
        END IF
C
        WRITE(*,100)'Swap first bit (y/n) '
        CSWAP(1:1)=READC('y','yn')
        LSWAP=(CSWAP.EQ.'y')
C
        WRITE(*,101)'Comments (max. 255 characters, '//
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
            STWV=0.
            DISP=0.
          END IF
        ELSE
          STWV=0.
          DISP=0.
        END IF
        WRITE(*,*)
        WRITE(*,100)'> STWV: '
        WRITE(*,*) STWV
        WRITE(*,100)'> DISP: '
        WRITE(*,*) DISP
        WRITE(*,*)
C
        WRITE(*,100)'Save whole frame (y/n) '
        CWHOLE(1:1)=READC('y','yn')
        IF(CWHOLE.EQ.'n')THEN
          IF(NAXIS(2).GT.1)THEN
            WRITE(*,100)'First scan: '
            NS1=READILIM('1',1,NAXIS(2))
            WRITE(*,100)'Last  scan: '
            WRITE(CDUMMY,*) NAXIS(2)
            NS2=READILIM(CDUMMY,NS1,NAXIS(2))
          ELSE
            NS1=1
            NS2=1
          END IF
        ELSE
          NS1=1
          NS2=NAXIS(2)
        END IF
        NSCAN=NS2-NS1+1
C
        WRITE(*,100)'Output file name'
        OUTFILE=OUTFILEX(30,'@',NSCAN,NAXIS(1),STWV,DISP,1,.FALSE.)
        WRITE(*,100)'Please wait...'
C leemos la matriz de datos
        BYTEPIX=2
        NBREAD=NAXIS(1)*NAXIS(2)*BYTEPIX/2880
        NS0=0
        DO IBL=0,NBREAD
          IHEADER=IHEADER+1
          READ(20,'(A2880)',REC=IHEADER) HDU
          NIBL=IBL*1440
          DO J=0,1439
            I1=ICHAR(HDU(J*2+1:J*2+1))
            I2=ICHAR(HDU(J*2+2:J*2+2))
            N=NIBL+J
            IY=N/NAXIS(1)
            IX=N-NAXIS(1)*IY
            IF(LSWAP)THEN
              IF(I1.GT.127) I1=I1-256
            END IF
            ROW(IX+1)=REAL(I1*256+I2)
            ROW(IX+1)=ROW(IX+1)*BSCALE+BZERO
            IF(IX+1.EQ.NAXIS(1))THEN
              NS0=NS0+1
              IF((NS0.GE.NS1).AND.(NS0.LE.NS2))THEN
                WRITE(30) (ROW(K),K=1,NAXIS(1))
              END IF
              IF(IY+1.EQ.NAXIS(2)) GOTO 30
            END IF
          END DO
        END DO
C
30      WRITE(*,101)'  ..OK!'
        CLOSE(20)
        CLOSE(30)
        STOP
C------------------------------------------------------------------------------
200     CONTINUE
        WRITE(*,101)'Este programa no puede leer una extension '//
     +   'Binary Table'
220     DO I=1,36
          I0=(I-1)*80
          IF(HDU(I0+1:I0+3).EQ.'END') GOTO 222                    !fin de HDU's
        END DO
        IHEADER=IHEADER+1
        READ(20,'(A2880)',REC=IHEADER) HDU
        GOTO 220
222     WRITE(*,110)'> No. of HDUs read: ',IHEADER
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
