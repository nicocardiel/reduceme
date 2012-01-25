C------------------------------------------------------------------------------
C Version 24-March-1998                                        file: fitshead.f
C @ ncl & fjg
C------------------------------------------------------------------------------
Comment
C
C Program: fitshead
C Classification: miscellany
C Description: Reads a FITS file and creates an output file with keyword 
C information.
C
Comment
C
C Este programa lee un fichero en formato FITS normal y escribe algunos
C parametros fundamentales de la cabecera.
C
        PROGRAM FITSHEAD
        IMPLICIT NONE
C
        INTEGER MAXNAXIS
        PARAMETER(MAXNAXIS=9)               !valor maximo admisible para NAXIS
        INTEGER NRESERVK
        PARAMETER(NRESERVK=16)                  !numero de "reserved keywords"
C
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER TRUEBEG,TRUELEN
        INTEGER READILIM
C
        INTEGER NCMAX_LOCAL
        PARAMETER (NCMAX_LOCAL=2048)
C
        INTEGER BITPIX
        INTEGER NAXIS(0:MAXNAXIS)                               !NAXIS=NAXIS(0)
        INTEGER BLANK
        REAL BSCALE
        REAL BZERO
        REAL GAIN                                                     !ganancia
        REAL READNOIS                                         !ruido de lectura
        CHARACTER*20 BUNIT
        CHARACTER*20 RA,DEC,EQUINOX
        CHARACTER*20 PACKDATE
        CHARACTER*20 TIMSTART,TIMEND
C
        INTEGER I,J,I0,L0,K,IX,IY,LLL
        INTEGER I1,I2
        INTEGER L1,L2
        INTEGER BYTEPIX,NIBL,IBL,NBREAD,N
        INTEGER IHEADER
        INTEGER NLINES
        INTEGER ICHANGE
        REAL ROW(NCMAX_LOCAL)
        REAL DATAMIN,DATAMAX
        REAL ACIMUT,ALTURA
        REAL ZENDIST
        CHARACTER*1 CHIST,CSWAP,CLIST,CPACKDATE
        CHARACTER*8 RESERVK(NRESERVK),RESERVK_TEMPORAL
        CHARACTER*75 INFILE,OUTFILE
        CHARACTER*255 CLINEA
        CHARACTER*2880 HDU
        LOGICAL LSWAP
        LOGICAL LHIST,LEXTENSION
        LOGICAL LRESERVK(NRESERVK)
C
        COMMON/BLKCOMPUTE0/CPACKDATE
        COMMON/BLKCOMPUTE1/RA,DEC,EQUINOX,PACKDATE,TIMSTART,TIMEND
        COMMON/BLKCOMPUTE2/ACIMUT,ALTURA
        COMMON/BLKCOMPUTE4/RESERVK
C------------------------------------------------------------------------------
C Introducimos en el DATA las "reserved keywords" que seran buscadas en el
C formato FITS normal. NO CAMBIAR el orden de las palabras; el programa
C asume el orden establecido. Si se deseen buscar nuevas "keywords" basta
C con introducirlas al final del DATA. 
C IMPORTANTE: rellenar las "keywords" con blancos si su longitud es inferior
C a 8 caracteres.
        DATA (RESERVK(I),I=1,NRESERVK)/
     +   'BSCALE  ','BZERO   ','BUNIT   ','BLANK   ','OBJECT  ',
     +   'EXPTIME ','AIRMASS ','RA      ','DEC     ','GAIN    ',
     +   'READNOIS','PACKDATE','TIMSTART','TIMEND  ','EQUINOX ',
     +   'ZENDIST '/
C
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
        OUTFILEX=OUTFILEX
C
        THISPROGRAM='fitshead'
        CALL WELCOME('24-March-1998')
C
4       WRITE(*,101)'* The following reserved Keywords will be '//
     +   'required: '
        DO I=1,NRESERVK
          WRITE(*,'(I2,1X,A1,1X,A)')I,'-',RESERVK(I)
        END DO
        WRITE(*,100)'Keyword to be changed (0=NONE) '
        ICHANGE=READILIM('0',0,NRESERVK)
        IF(ICHANGE.NE.0)THEN
          WRITE(*,'(A)')'The following keyword will be changed: '//
     +     RESERVK(ICHANGE)
          WRITE(*,100)'New keyword'
          RESERVK_TEMPORAL(1:8)=READC('@','@')
          RESERVK(ICHANGE)=RESERVK_TEMPORAL
          GOTO 4
        END IF
C
        DO I=1,NRESERVK
          LRESERVK(I)=.FALSE.  !inicializamos: ninguna "reserved keyword found" 
        END DO
C
        BSCALE=-99.
        BZERO=-99.
        BUNIT='XXXXXXXXXXXXXXXXXXXX'
        BLANK=99
        OBJECT='XXXXXXXXXXXXXXXXXXXX'
        TIMEXPOS=-99.
        AIRMASS=-99.
        RA='XX:XX:XX.XX'
        DEC='XXX:XX:XX.X'
        GAIN=-99.
        READNOIS=-99.
        PACKDATE='XX/XX/XX'
        TIMSTART='XX:XX:XX.XXX'
        TIMEND='XX:XX:XX.XXX'
        EQUINOX='XXXXX.X'
        ZENDIST=-99
C
        WRITE(*,100)'FITS file name'
        INFILE=INFILEX(20,'@',0,0,0.,0.,4,.FALSE.)
        FITSFILE=INFILE
C
        WRITE(*,100)'Output file name (header information)'
        OUTFILE(1:75)=READC('@','@')
        OPEN(30,FILE=OUTFILE,STATUS='UNKNOWN',FORM='FORMATTED')
        WRITE(*,100)'[l]ong, [s]hort list or s[p]ecial (l/s/p) '
        CLIST(1:1)=READC('l','lsp')
        IF(CLIST.EQ.'p')THEN
          WRITE(*,*)
          WRITE(*,100)'Date will be obtained from: '
          WRITE(*,100)RESERVK(12)
          WRITE(*,*)
          WRITE(*,101)'(1) format YY/MM/DD'
          WRITE(*,101)'(2) format DD/MM/YY'
          WRITE(*,100)'Option (1/2)'
          CPACKDATE(1:1)=READC('@','12')
        END IF
        NLINES=0
5       READ(30,'(A)',END=6) CLINEA
        NLINES=NLINES+1
        GOTO 5
6       IF(NLINES.EQ.0)THEN
          IF(CLIST.EQ.'l')THEN
            WRITE(30,101)'FITS--file BPIX  NAX    N1    N2 OBJECT'//
     +       '                   R.A.         DEC.           MIN '//
     +       '     MAX   BLANK   BSCALE   BZERO       AIRMASS    '//
     +       ' TIMEXPOS BUNIT  GAIN READNOIS'
          ELSEIF(CLIST.EQ.'s')THEN
            WRITE(30,101)'FITS--file    N1    N2 OBJECT          '//
     +       '   MIN      MAX      AIRMASS    TIMEXPOS'
          ELSE
            WRITE(30,101)'FITS--file BPIX  NAX    N1    N2 OBJECT'//
     +       '                   R.A.         DEC.     '//
     +       '       AIRMASS    TIMEXPOS   DATE       UT1'//
     +       '          UT2         A       h       z'
          END IF
        END IF
C
        WRITE(*,100)'Run silently (y/n) '
        CHIST(1:1)=READC('n','yn')
        LHIST=(CHIST.EQ.'n')
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
        LEXTENSION=.FALSE.      !suponemos que inicialmente no es una EXTENSION
10      READ(HDU(91:110),*)BITPIX
        READ(HDU(171:190),*)NAXIS(0)
        IF(LHIST)THEN
          WRITE(*,110)'> BITPIX  = ',BITPIX
          WRITE(*,110)'> NAXIS   = ',NAXIS(0)
        END IF
        L1=TRUEBEG(FITSFILE)
        L2=TRUELEN(FITSFILE)
        IF(L2-L1+1.GT.10) L1=L2-9
        WRITE(30,'(A10,$)')FITSFILE(L1:L2)
        IF((CLIST.EQ.'l').OR.(CLIST.EQ.'p'))
     +   WRITE(30,'(1X,I4,1X,I4,$)') BITPIX,NAXIS(0)
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
!---------------------------->123456789012345678901234567890
          IF(HDU(241:270).NE.'EXTEND  =                    T')THEN
            WRITE(*,101)'> FATAL ERROR: Standard extension not '//
     +       'confirmed.'
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
        DO I=1,NAXIS(0)
          READ(HDU((I+2)*80+11:(I+2)*80+30),*)NAXIS(I)
          IF(LHIST) WRITE(*,'(A,I1,A,I6)')'> NAXIS',I,'  = ',NAXIS(I)
          WRITE(30,'(1X,I5,$)') NAXIS(I)
        END DO
        IF(NAXIS(1).GT.NCMAX_LOCAL)THEN
          WRITE(*,101)'FATAL ERROR: NAXIS1.GT.NCMAX'
          CLOSE(20)
          STOP
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
          DO K=1,NRESERVK                     !buscamos cada "reserved keyword"
            IF(.NOT.LRESERVK(K))THEN             !solo leemos primera aparicion 
              IF(HDU(I0+1:I0+8).EQ.RESERVK(K))THEN
                LRESERVK(K)=.TRUE.
                IF(LHIST) WRITE(*,101)'> '//HDU(I0+1:I0+30)
                IF(K.EQ.1)THEN
                  IF(TRUELEN(HDU(I0+11:I0+30)).GT.0)THEN
                    READ(HDU(I0+11:I0+30),*) BSCALE
                  ELSE
                    WRITE(*,101)'ERROR: '//RESERVK(K)//
     +               ' --> invalid field'
                  END IF
                ELSEIF(K.EQ.2)THEN
                  IF(TRUELEN(HDU(I0+11:I0+30)).GT.0)THEN
                    READ(HDU(I0+11:I0+30),*) BZERO
                  ELSE
                    WRITE(*,101)'ERROR: '//RESERVK(K)//
     +               ' --> invalid field'
                  END IF
                ELSEIF(K.EQ.3)THEN
                  LLL=INDEX(HDU(I0+12:I0+30),CHAR(39))-1
                  BUNIT=HDU(I0+12:I0+11+LLL)
                ELSEIF(K.EQ.4)THEN
                  IF(TRUELEN(HDU(I0+11:I0+30)).GT.0)THEN
                    READ(HDU(I0+11:I0+30),*) BLANK
                  ELSE
                    WRITE(*,101)'ERROR: '//RESERVK(K)//
     +               ' --> invalid field'
                  END IF
                ELSEIF(K.EQ.5)THEN
                  LLL=INDEX(HDU(I0+12:I0+30),CHAR(39))-1
                  OBJECT=HDU(I0+12:I0+11+LLL)
                ELSEIF(K.EQ.6)THEN
                  IF(TRUELEN(HDU(I0+11:I0+30)).GT.0)THEN
                    READ(HDU(I0+11:I0+30),*) TIMEXPOS
                  ELSE
                    WRITE(*,101)'ERROR: '//RESERVK(K)//
     +               ' --> invalid field'
                  END IF
                ELSEIF(K.EQ.7)THEN
                  IF(TRUELEN(HDU(I0+11:I0+30)).GT.0)THEN
                    READ(HDU(I0+11:I0+30),*) AIRMASS
                  ELSE
                    WRITE(*,101)'ERROR: '//RESERVK(K)//
     +               ' --> invalid field'
                  END IF
                ELSEIF(K.EQ.8)THEN
                  LLL=INDEX(HDU(I0+12:I0+30),CHAR(39))-1
                  RA=HDU(I0+12:I0+11+LLL)
                ELSEIF(K.EQ.9)THEN
                  LLL=INDEX(HDU(I0+12:I0+30),CHAR(39))-1
                  DEC=HDU(I0+12:I0+11+LLL)
                ELSEIF(K.EQ.10)THEN
                  IF(TRUELEN(HDU(I0+11:I0+30)).GT.0)THEN
                    READ(HDU(I0+11:I0+30),*) GAIN
                  ELSE
                    WRITE(*,101)'ERROR: '//RESERVK(K)//
     +               ' --> invalid field'
                  END IF
                ELSEIF(K.EQ.11)THEN
                  IF(TRUELEN(HDU(I0+11:I0+30)).GT.0)THEN
                    READ(HDU(I0+11:I0+30),*) READNOIS
                  ELSE
                    WRITE(*,101)'ERROR: '//RESERVK(K)//
     +               ' --> invalid field'
                  END IF
                ELSEIF(K.EQ.12)THEN
                  LLL=INDEX(HDU(I0+12:I0+30),CHAR(39))-1
                  PACKDATE=HDU(I0+12:I0+11+LLL)
                ELSEIF(K.EQ.13)THEN
                  LLL=INDEX(HDU(I0+12:I0+30),CHAR(39))-1
                  TIMSTART=HDU(I0+12:I0+11+LLL)
                ELSEIF(K.EQ.14)THEN
                  LLL=INDEX(HDU(I0+12:I0+30),CHAR(39))-1
                  TIMEND=HDU(I0+12:I0+11+LLL)
                ELSEIF(K.EQ.15)THEN
                  LLL=INDEX(HDU(I0+12:I0+30),CHAR(39))-1
                  EQUINOX=HDU(I0+12:I0+11+LLL)
                ELSEIF(K.EQ.16)THEN
                  IF(TRUELEN(HDU(I0+11:I0+30)).GT.0)THEN
                    READ(HDU(I0+11:I0+30),*) ZENDIST
                  ELSE
                    WRITE(*,101)'ERROR: '//RESERVK(K)//
     +               ' --> invalid field'
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
21        CONTINUE
        END DO
        L0=1
        IHEADER=IHEADER+1
        READ(20,'(A2880)',REC=IHEADER) HDU
        GOTO 20
C
22      IF(LHIST) WRITE(*,110)'> No. of HDUs read: ',IHEADER
        DO I=1,NRESERVK
          IF(.NOT.LRESERVK(I))THEN
            WRITE(*,101) 'WARNING: '//RESERVK(I)//
     +       ' has not been found.'
          END IF
        END DO
C
C Si NAXIS > 2 el programa hay que modificarlo
        IF(NAXIS(0).GT.2)THEN
          IF((NAXIS(0).EQ.3).AND.(NAXIS(3).EQ.1))THEN
            WRITE(*,101)'WARNING: NAXIS.GT.2 but NAXIS(3)=1'
          ELSE
            WRITE(*,101)'FATAL ERROR: NAXIS.GT.2'
            WRITE(*,101)'This situation has not been implemented.'//
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
        IF(CLIST.EQ.'p') GOTO 30
C
        WRITE(*,100)'Swap first bit (y/n) '
        CSWAP(1:1)=READC('y','yn')
        LSWAP=(CSWAP.EQ.'y')
C
        DATAMIN=1.E20
        DATAMAX=-1.E20
C leemos la matriz de datos
        BYTEPIX=2
        NBREAD=NAXIS(1)*NAXIS(2)*BYTEPIX/2880
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
              DO K=1,NAXIS(1)
                IF(ROW(K).LT.DATAMIN) DATAMIN=ROW(K)
                IF(ROW(K).GT.DATAMAX) DATAMAX=ROW(K)
              END DO
              IF(IY+1.EQ.NAXIS(2)) GOTO 30
            END IF
          END DO
        END DO
C 
30      CLOSE(20)
C
        IF(CLIST.EQ.'l')THEN
          WRITE(30,'(1X,A20,2(1X,A13),2(1X,F8.0),1X,I4,1X,F10.4,1X,
     +     F10.2,1X,F12.7,1X,F10.2,1X,A,1X,F7.3,2X,F7.3)')
     +     OBJECT,RA,DEC,DATAMIN,DATAMAX,BLANK,BSCALE,BZERO,AIRMASS,
     +     TIMEXPOS,BUNIT(1:TRUELEN(BUNIT)),GAIN,READNOIS
        ELSEIF(CLIST.EQ.'s')THEN
          WRITE(30,'(1X,A14,2(1X,F8.0),1X,F12.6,1X,F10.2)')
     +     OBJECT,DATAMIN,DATAMAX,AIRMASS,
     +     TIMEXPOS
        ELSEIF(CLIST.EQ.'p')THEN
          IF(TIMEND(1:2).EQ.'XX')THEN
            IF(TIMSTART(1:2).NE.'XX')THEN
              IF(TIMEXPOS.NE.-99.)THEN
                CALL NEWTIMEND(TIMSTART,TIMEXPOS,TIMEND)
              END IF
            END IF
          END IF
          CALL COMPUTE
          WRITE(30,'(1X,A20,2(1X,A13),1X,F12.7,1X,
     +     F10.2,3(1X,A),4(1X,F7.3))')
     +     OBJECT,RA,DEC,AIRMASS,
     +     TIMEXPOS,
     +     PACKDATE(1:8),
     +     TIMSTART(1:12),
     +     TIMEND(1:12),
ccc     +     PACKDATE(1:TRUELEN(PACKDATE)),
ccc     +     TIMSTART(1:TRUELEN(TIMSTART)),
ccc     +     TIMEND(1:TRUELEN(TIMEND)),
     +     ACIMUT,ALTURA,ZENDIST,90.-ALTURA-ZENDIST
        ELSE
          WRITE(*,101)'FATAL ERROR: invalid CLIST'
          STOP
        END IF
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
C
C******************************************************************************
C
        SUBROUTINE COMPUTE
        IMPLICIT NONE
        INCLUDE 'futils.inc'
        INTEGER TRUELEN
C
        INTEGER NRESERVK
        PARAMETER(NRESERVK=16)                  !numero de "reserved keywords"
C
        CHARACTER*8 RESERVK(NRESERVK)
        CHARACTER*20 RA,DEC,EQUINOX,PACKDATE,TIMSTART,TIMEND
        CHARACTER*1 CPACKDATE
C
        REAL PI
        PARAMETER(PI=3.141593)
C
        INTEGER ANO,MES,DIA
        INTEGER IXX1,IXX2
        INTEGER K,KK
        INTEGER L1,L2
        INTEGER LCOLON
        REAL MONTH
        REAL TI,TF,TII,TFF
        REAL X(3),X0(3)
        REAL GIO,ZETA,TETA
        REAL M(3,3)
        REAL ARINI,DECINI,ARFIN,DECFIN,HOR
        REAL ACIMUT,ALTURA
        REAL FXX3
        REAL HORA,TS0,TU,TU1,TU2,TSG,TSL
        REAL LAT,LONG,ALTOBS
        CHARACTER*1 CKEYB
        CHARACTER*15 OBS
        DOUBLE PRECISION FJ
C
        COMMON/BLKCOMPUTE0/CPACKDATE
        COMMON/BLKCOMPUTE1/RA,DEC,EQUINOX,PACKDATE,TIMSTART,TIMEND
        COMMON/BLKCOMPUTE2/ACIMUT,ALTURA
        COMMON/BLKCOMPUTE4/RESERVK
C------------------------------------------------------------------------------
        IF((RA(1:2).EQ.'XX').OR.
     +     (DEC(1:2).EQ.'XX').OR.
     +     (EQUINOX(1:2).EQ.'XX').OR.
     +     (PACKDATE(1:2).EQ.'XX').OR.
     +     (TIMSTART(1:2).EQ.'XX').OR.
     +     (TIMEND(1:2).EQ.'XX'))THEN
          WRITE(*,101)'WARNING: insufficient Keyword information to '//
     +     'compute telescope position.'
          WRITE(*,100)'Do you want to introduce the missing parameters '
          WRITE(*,100)'through keyboard (y/n) '
          CKEYB(1:1)=READC('y','yn')
          IF(CKEYB.EQ.'n') RETURN
          IF(RA(1:2).EQ.'XX')THEN
            WRITE(*,100)'Enter '
            WRITE(*,100) RESERVK(8)(1:TRUELEN(RESERVK(8)))
            WRITE(*,100)' (HH:MM:SS.ddd) '
            RA(1:20)=READC('00:00:00.000','@')
          END IF
          IF(DEC(1:2).EQ.'XX')THEN
            WRITE(*,100)'Enter '
            WRITE(*,100) RESERVK(9)(1:TRUELEN(RESERVK(9)))
            WRITE(*,100)' (+DD:MM:SS.dd) '
            DEC(1:20)=READC('+00:00:00.00','@')
          END IF
          IF(EQUINOX(1:2).EQ.'XX')THEN
            WRITE(*,100)'Enter '
            WRITE(*,100) RESERVK(15)(1:TRUELEN(RESERVK(15)))
            WRITE(*,100)' (B1950.0 or J2000.0) '
            EQUINOX(1:20)=READC('J2000.0','@')
          END IF
          IF(PACKDATE(1:2).EQ.'XX')THEN
            WRITE(*,100)'Enter '
            WRITE(*,100) RESERVK(12)(1:TRUELEN(RESERVK(12)))
            WRITE(*,100)' (DD/MM/YY) '
            PACKDATE(1:20)=READC('01/01/99','@')
          END IF
          IF(TIMSTART(1:2).EQ.'XX')THEN
            WRITE(*,100)'Enter '
            WRITE(*,100) RESERVK(13)(1:TRUELEN(RESERVK(13)))
            WRITE(*,100)' (HH:MM:SS.ddd) '
            TIMSTART(1:20)=READC('00:00:00.000','@')
          END IF
          IF(TIMEND(1:2).EQ.'XX')THEN
            WRITE(*,100)'Enter '
            WRITE(*,100) RESERVK(14)(1:TRUELEN(RESERVK(14)))
            WRITE(*,100)' (HH:MM:SS.ddd) '
            TIMEND(1:20)=READC('00:00:00.000','@')
          END IF
        END IF
C------------------------------------------------------------------------------
C elegimos observatorio
        CALL OBSERVAT(OBS,LAT,LONG,ALTOBS)
        CALL DECIMAL(LAT,LAT)
        CALL DECIMAL(LONG,LONG)
C------------------------------------------------------------------------------
C extraemos la informacion necesaria de las Keywords
        IF(CPACKDATE.EQ.'1')THEN
          READ(PACKDATE,'(I2)') ANO
          READ(PACKDATE,'(3X,I2)') MES
          READ(PACKDATE,'(6X,I2)') DIA
        ELSE
          READ(PACKDATE,'(I2)') DIA
          READ(PACKDATE,'(3X,I2)') MES
          READ(PACKDATE,'(6X,I2)') ANO
        END IF
        WRITE(*,100)'> Year.: '
        ANO=ANO+1900
        WRITE(*,*) ANO
        WRITE(*,100)'> Month: '
        WRITE(*,*) MES
        WRITE(*,100)'> Day..: '
        WRITE(*,*) DIA
C
        WRITE(*,100)'> UT initial: '
        READ(TIMSTART,'(I2,1X,I2,1X,F6.3)')IXX1,IXX2,FXX3
        WRITE(*,'(I2,1X,I2,1X,F6.3)')IXX1,IXX2,FXX3
        TU1=REAL(IXX1)+REAL(IXX2)/60.+FXX3/3600.
        WRITE(*,100)'> UT final..: '
        READ(TIMEND,'(I2,1X,I2,1X,F6.3)')IXX1,IXX2,FXX3
        WRITE(*,'(I2,1X,I2,1X,F6.3)')IXX1,IXX2,FXX3
        TU2=REAL(IXX1)+REAL(IXX2)/60.+FXX3/3600.
        IF(TU2.LT.TU1) TU2=TU2+24.
        TU=(TU1+TU2)/2.
        IF(TU.GT.24.) TU=TU-24.
        WRITE(*,100)'> Mean UT...: '
        WRITE(*,*) TU
C
        READ(EQUINOX,'(1X,F6.1)') TII
C
        WRITE(*,100)'> R.A. ('//EQUINOX(1:TRUELEN(EQUINOX))//'): '
        LCOLON=INDEX(RA,':')
        IF(LCOLON.EQ.3)THEN
          READ(RA,'(I2,1X,I2,1X,F5.2)')IXX1,IXX2,FXX3
        ELSEIF(LCOLON.EQ.4)THEN
          READ(RA,'(1X,I2,1X,I2,1X,F5.2)')IXX1,IXX2,FXX3
        ELSE
          WRITE(*,100)'LCOLON: '
          WRITE(*,*)LCOLON
          STOP'FATAL ERROR: invalid LCOLON'
        END IF
        WRITE(*,'(I2.2,1X,I2.2,1X,F5.2)')IXX1,IXX2,FXX3
        ARINI=REAL(IXX1)+REAL(IXX2)/60.+FXX3/3600.
        WRITE(*,100)'> DEC. ('//EQUINOX(1:TRUELEN(EQUINOX))//'): '
        L1=INDEX(DEC,':')
        L2=L1+INDEX(DEC(L1+1:),':')
        IF((L1.EQ.0).OR.(L2.EQ.0))THEN
          WRITE(*,100)'FATAL ERROR: unknown declination format.'
          STOP
        END IF
        IF((DEC(1:1).EQ.'-').OR.(DEC(1:1).EQ.'+'))THEN
          READ(DEC(2:L1-1),*)IXX1
          READ(DEC(L1+1:L2-1),*)IXX2
          READ(DEC(L2+1:),*)FXX3
          WRITE(*,'(A1,I2.2,1X,I2.2,1X,F5.2)')DEC(1:1),IXX1,IXX2,FXX3
          DECINI=REAL(IXX1)+REAL(IXX2)/60.+FXX3/3600.
          IF(DEC(1:1).EQ.'-') DECINI=-DECINI
        ELSE
          READ(DEC(1:L1-1),*)IXX1
          READ(DEC(L1+1:L2-1),*)IXX2
          READ(DEC(L2+1:),*)FXX3
          WRITE(*,'(A1,I2.2,1X,I2.2,1X,F5.2)')' ',IXX1,IXX2,FXX3
          DECINI=REAL(IXX1)+REAL(IXX2)/60.+FXX3/3600.
        END IF
C
        IF(MES.LE.2)THEN
          MONTH=AINT(REAL(MES-1)*63./2.)
        ELSE
          MONTH=AINT(REAL(MES+1)*30.6)-63.
        END IF
        MONTH=MONTH+REAL(DIA)
        TFF=REAL(ANO)+MONTH/365.25
        TI=(TII-2000.)/100.
        TF=(TFF-2000.-100.*TI)/100.
C elementos precesionales en grados
        GIO=((2306.2181+1.39656*TI-.000139*TI*TI)*TF+(.30188-.000344*TI)
     +    *TF*TF+.017998*TF*TF*TF)/3600.
        ZETA=GIO+((.7928+.00041*TI)*TF*TF+.000205*TF*TF*TF)/3600.
        TETA=((2004.3109-.8533*TI-.000217*TI*TI)*TF-(.42665+.000217*TI)
     +     *TF*TF-.041833*TF*TF*TF)/3600.
C los pasamos a radianes
        GIO=GIO*PI/180.
        ZETA=ZETA*PI/180.
        TETA=TETA*PI/180.
C matriz de rotacion
        M(1,1)=-SIN(GIO)*SIN(ZETA)+COS(GIO)*COS(TETA)*COS(ZETA)
        M(1,2)=-COS(GIO)*SIN(ZETA)-SIN(GIO)*COS(ZETA)*COS(TETA)
        M(1,3)=-SIN(TETA)*COS(ZETA)
        M(2,1)=SIN(GIO)*COS(ZETA)+COS(GIO)*COS(TETA)*SIN(ZETA)
        M(2,2)=COS(GIO)*COS(ZETA)-SIN(GIO)*COS(TETA)*SIN(ZETA)
        M(2,3)=-SIN(TETA)*SIN(ZETA)
        M(3,1)=COS(GIO)*SIN(TETA)
        M(3,2)=-SIN(GIO)*SIN(TETA)
        M(3,3)=COS(TETA)
C coordenadas rectangulares del objeto
        X0(1)=COS(DECINI*PI/180.)*COS(ARINI*15.*PI/180.)
        X0(2)=COS(DECINI*PI/180.)*SIN(ARINI*15.*PI/180.)
        X0(3)=SIN(DECINI*PI/180.)
C cambio a coordenadas de la epoca
        DO K=1,3
          X(K)=0.
          DO KK=1,3
            X(K)=X(K)+X0(KK)*M(K,KK)
          END DO
        END DO
        ARFIN=ATAN2(X(2),X(1))
        DECFIN=ASIN(X(3))
C coordenadas en grados
        ARFIN=ARFIN*180./PI/15.
        IF(ARFIN.LT.0.)ARFIN=ARFIN+24.
        DECFIN=DECFIN*180./PI
C calculamos la fecha juliana y el TS a 0h TU en Greenwich (TS0)
        HORA=0.
        CALL FJ_TS0(REAL(ANO),REAL(MES),REAL(DIA),HORA,FJ,TS0)
        WRITE(*,100)'Julian Date: '
        WRITE(*,*) FJ
        WRITE(*,100)'TS0 (TS in Greenwich at 0h UT): '
        WRITE(*,*) TS0
C TS en Greenwich en el instante medio de observacion TU
        TSG=TS0+TU+TU*2.737909E-3
        IF(TSG.GT.24.) TSG=TSG-24.
        WRITE(*,100)'TSG (TS in Greenwich).........: '
        WRITE(*,*) TSG
C TSL para el instante TU
        TSL=TSG+LONG/15.
        IF(TSL.GT.24.) TSL=TSL-24.
        IF(TSL.LT.0.) TSL=TSL+24.
        WRITE(*,100)'TSL (Local TS)................: '
        WRITE(*,*) TSL
C datos del observatorio elegido
        WRITE(*,101)'> Observatory: '//OBS(1:TRUELEN(OBS))
        WRITE(*,100)'> Longitude..: '
        WRITE(*,*) LONG
        WRITE(*,100)'> Latitude...: '
        WRITE(*,*) LAT
        WRITE(*,100)'> Height.....: '
        WRITE(*,*) ALTOBS
C angulo horario
        HOR=TSL-ARFIN
        IF(HOR.LT.0.) HOR=HOR+24.
C pasamos a coordenadas horizontales
        CALL CAMBCOOR(HOR,DECFIN,ACIMUT,ALTURA,LAT)
        WRITE(*,100)'> Azimut: '
        WRITE(*,*) ACIMUT
        WRITE(*,100)'> Height: '
        WRITE(*,*) ALTURA
C
C------------------------------------------------------------------------------
100        FORMAT(A,$)
101        FORMAT(A)
        END
C
C **********************************************************************
C                                                      SUBROUTINE FJ_TS0
C                                                      *****************
      SUBROUTINE FJ_TS0(ANO,MES,DIA,HORA,FJ,TS0)
C
C Calcula la fecha juliana y el Tiempo Sidereo en Greenwich a 0h UT a
C partir de los datos de una fecha concreta. En la fecha juliana se
C tiene en cuenta la fraccion del dia correspondiente a la hora
C traspasada a traves de la variable global HORA.
C Los calculos se realizan en doble precision para no perder informacion
C por redondeo.
C
      IMPLICIT NONE
C---> variables globales: ENTRADA
      REAL ANO,MES,DIA,HORA
C---> variables globales: SALIDA
      REAL TS0
      DOUBLE PRECISION FJ
C---> variables locales
      INTEGER*4 FECHA
      DOUBLE PRECISION M,A,B,Y,DT,DTS0
      DOUBLE PRECISION DANO,DMES,DDIA
C-----------------------------------------------------------------------
      DANO=DBLE(ANO)
      DMES=DBLE(MES)
      DDIA=DBLE(DIA)
      FECHA=INT(ANO)*10000+INT(MES)*100+INT(DIA)
      IF(FECHA.GE.15821015)THEN
        A=DINT(DANO/1.D2)
        B=2.D0-A+DINT(A/4.D0)
      ELSE
        B=0.D0
      END IF
      IF(INT(ANO).GE.0)THEN
        A=0.D0
      ELSE
        A=-.75D0
      END IF
      IF(INT(MES).GE.2)THEN
        Y=DANO
        M=DMES
      ELSE
        Y=DANO-1.D0
        M=DMES+1.2D1
      END IF
      FJ=DINT(365.25D0*Y+A)+DINT(30.6001D0*(M+1.D0))+DDIA+B+
     +   1720994.5D0
      DT=(FJ-2451545.D0)/36525.D0
      DTS0=6.D0*3.6D3+41.D0*6.D1+50.54851D0+8640184.812866D0*DT
     +     +.093104D0*DT*DT-6.2D-6*DT*DT*DT
      DTS0=DMOD(DTS0,8.64D4)
      DTS0=DTS0/8.64D4
      DTS0=DTS0*2.4D1
      TS0=REAL(DTS0)
      IF(TS0.LT.0.)TS0=TS0+24.
      IF(TS0.GT.24.)TS0=TS0-24.
      FJ=FJ+DBLE(HORA/24.)
      END
C
C **********************************************************************
C                                                    SUBROTUINE OBSERVAT
C                                                    *******************
      SUBROUTINE OBSERVAT(OBS,LAT,LONG,ALTOBS)
C
C Longitud y latitud (grados), y altura (metros) de los observatorios
C con el formato DD.MMSS
C Por convenio longitudes son positivas hacia el Este y negativas hacia
C el Oeste.
C
      IMPLICIT NONE
      INCLUDE 'futils.inc'
      INTEGER READI
      REAL READF
C---> argumentos ficticios: SALIDA
      REAL LAT,LONG
      REAL ALTOBS
      CHARACTER*15 OBS
C---> parametros locales
      INTEGER NOBS
      PARAMETER(NOBS=6)
C---> variables locales
      INTEGER I,NOPC
      CHARACTER*15 OBSER(NOBS)
C
      DATA OBSER/'CALAR ALTO','IZAÑA','MADRID','LA PALMA',
     + 'SAN PEDRO','LICK'/
C
      WRITE(*,*)
      WRITE(*,160)0,'Enter data from keyboard (Longitude,'//
     +              'Latitude,Height)'
      DO I=1,NOBS
        WRITE(*,160)I,OBSER(I)
      END DO
    5 CONTINUE
      WRITE(*,100)'Observatory number '
      NOPC=READI('4')
      IF((NOPC.LT.0).OR.(NOPC.GT.NOBS))THEN
        WRITE(*,101)'ERROR: option out of range. Try again.'
        GOTO 5
      END IF
      WRITE(*,*)
C
      IF(NOPC.EQ.0)THEN
        WRITE(*,*)
        WRITE(*,100)'Observatory name (max. 15 characters)'
        OBS(1:15)=READC('@','@')
        WRITE(*,100)'Latitude (GG.MMSS)'
        LAT=READF('@')
        WRITE(*,100)'Longitude (GG.MMSS, +E, -W)'
        LONG=READF('@')
        WRITE(*,100)'Height (metres)'
        ALTOBS=READF('@')
      ELSE IF(NOPC.EQ.1)THEN
C Calar Alto
        LAT=37.1305
        LONG=-2.3240
        ALTOBS=2165.
      ELSE IF(NOPC.EQ.2)THEN
C Izaña
        LAT=28.1732
        LONG=-16.2945
        ALTOBS=2400.
      ELSE IF(NOPC.EQ.3)THEN
C Madrid
        LAT=40.2430
        LONG=-3.4115
        ALTOBS=656.
      ELSE IF(NOPC.EQ.4)THEN
C La Palma
        LAT=28.4540
        LONG=-17.5247
        ALTOBS=2334.
      ELSE IF(NOPC.EQ.5)THEN
C San Pedro (datos aproximados extraidos de un mapa mundi)
        LAT=32.0
        LONG=-115
        ALTOBS=2000.
      ELSE IF(NOPC.EQ.6)THEN
C Lick (datos facilitados por Jesus Gallego)
        LAT=37.2036
        LONG=-123.3812
        ALTOBS=1283.
      END IF
      IF(NOPC.NE.0)OBS=OBSER(NOPC)
  100 FORMAT(A,$)
  101 FORMAT(A)
  160 FORMAT('(',I1,')',1X,A)
      END
C **********************************************************************
C                                                    SUBROUTINE CAMBCOOR
C                                                    *******************
      SUBROUTINE CAMBCOOR(HOR,DEC,ACI,ALT,LAT)
C
C Transformaciones de coordenadas:
C (angulo horario,declinacion) ----> (acimut,altura)
C
      IMPLICIT NONE
C---> parametros
      REAL PI
      PARAMETER(PI=3.141593)
C---> argumentos ficticios: ENTRADA
      REAL HOR,DEC
C---> argumentos ficticios: SALIDA
      REAL ACI,ALT
C---> latitud
      REAL LAT
C---> variables locales
      REAL HORR,DECR
      REAL LATR
      REAL SENOAC,COSEAC
C
      HORR=HOR*PI/180.*15.
      DECR=DEC*PI/180.
      LATR=LAT*PI/180.
      SENOAC=COS(DECR)*SIN(HORR)
      COSEAC=-SIN(DECR)*COS(LATR)+COS(DECR)*SIN(LATR)*COS(HORR)
      ACI=ATAN2(SENOAC,COSEAC)
      ACI=ACI*180./PI
      IF(ACI.LT.0.)ACI=ACI+360.
      ALT=ASIN(SIN(LATR)*SIN(DECR)+COS(LATR)*COS(DECR)*COS(HORR))
      ALT=ALT*180./PI
      END
C
C **********************************************************************
C                                                     SUBROUTINE DECIMAL
C                                                     ******************
      SUBROUTINE DECIMAL(A,B)
C
C Paso de DD.MMSS a DD.dddd
C Utilizamos variables INTEGER*4 para evitar los errores de redondeo.
C Es posible hacer una llamada a la subrutina en la que los argumentos
C verdaderos correspondientes a las argumentos ficticios A y B
C coincidan.
C
      IMPLICIT NONE
C---> argumentos ficticios: ENTRADA
      REAL A
C---> argumentos ficticios: SALIDA
      REAL B
C---> variables locales
      INTEGER*4 AA,A1,A2,A3
      CHARACTER*1 SIGNO
C
      SIGNO='+'
      IF(A.LT.0)SIGNO='-'
      AA=ABS(NINT(A*10000))
      A1=AA/10000
      A2=AA-A1*10000
      A2=A2/100
      A3=AA-A1*10000-A2*100
      B=REAL(A1)+REAL(A2)/60.+REAL(A3)/3600.
      IF(SIGNO.EQ.'-')B=-B
      END
C
C **********************************************************************
C                                                       SUBROUTINE SEXAG
C                                                       ****************
      SUBROUTINE SEXAG(A,B,C,D)
C
C Paso de DD.dddd a DD.MMSS
C
      IMPLICIT NONE
C---> argumentos ficticios: ENTRADA
      REAL A
C---> argumentos ficticios: SALIDA
      REAL B,C,D
C---> variables locales
      REAL AA,BR
C
      AA=ABS(A)
      B=AINT(AA)
      BR=(AA-B)*60.
      C=AINT(BR)
      D=ANINT((BR-AINT(BR))*60)
      IF(D.GE.60.)THEN
        C=C+1.
        D=0.
        IF(C.GE.60.)THEN
          B=B+1.
          C=0.
        END IF
      END IF
C hay que tener cuidado con la siguiente asignacion de signos
C (necesaria por si B=0 o/y C=0)
      IF(A.LT.0.)THEN
        B=-B
        C=-C
        D=-D
      END IF
      END
C
C******************************************************************************
C Calcula TIMEND a partir de TIMSTART y TIMEXPOS
C
        SUBROUTINE NEWTIMEND(TIMSTART,TIMEXPOS,TIMEND)
        IMPLICIT NONE
C
        CHARACTER*20 TIMSTART,TIMEND
        REAL TIMEXPOS
C
        INTEGER HH,MM
        REAL SS
        REAL HORA
        REAL FHH,FMM,FSS
C------------------------------------------------------------------------------
        READ(TIMSTART,'(I2,1X,I2,1X,F6.3)')HH,MM,SS
        HORA=REAL(HH)+REAL(MM)/60.+REAL(SS)/3600.
        HORA=HORA+TIMEXPOS/3600.
        CALL SEXAG(HORA,FHH,FMM,FSS)
        WRITE(TIMEND,'(I2,A1,I2,A1,F6.3)')INT(FHH),':',INT(FMM),
     +   ':',FSS
        END
