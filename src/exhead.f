C------------------------------------------------------------------------------
C Version 24-March-1998                                          file: exhead.f
C @ ncl & fjg
C------------------------------------------------------------------------------
Comment
C
C Program: exhead
C Classification: examination & statistics
C Description: Shows and changes image header information. Image data are not
C modified.
C
Comment
C
C Lee la cabecera de una imagen para obtener la informacion de su
C tamanho. Este programa no realiza ninguna modificacion sobre el
C fichero.
C
        PROGRAM EXHEAD
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER TRUELEN
        INTEGER READILIM
        REAL READF
C
        INTEGER I,J,L
        REAL A(NCMAX,NSMAX)
        CHARACTER*1 CCHANGE,CCC
        CHARACTER*12 CLAVE
        CHARACTER*50 CDUMMY
        CHARACTER*80 INFILE,OUTFILE
        LOGICAL LOGFILE
C------------------------------------------------------------------------------
        INFILEX=INFILEX
C
        THISPROGRAM='exhead'
        CALL WELCOME('24-March-1998')
C
        DO I=1,LEN(OBJECT)
          OBJECT(I:I)=CHAR(32)
        END DO
        DO I=1,LEN(FITSFILE)
          FITSFILE(I:I)=CHAR(32)
        END DO
        DO I=1,LEN(COMMENT)
          COMMENT(I:I)=CHAR(32)
        END DO
C
        WRITE(*,100)'File name to be examined'
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
          WRITE(CDUMMY,*)STWV
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          WRITE(*,101)'STWV    : '//CDUMMY(1:L)
          WRITE(CDUMMY,*)DISP
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          WRITE(*,101)'DISP    : '//CDUMMY(1:L)
          READ(20)AIRMASS
          WRITE(CDUMMY,*)AIRMASS
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          WRITE(*,101)'Airmass : '//CDUMMY(1:L)
          READ(20)TIMEXPOS
          WRITE(CDUMMY,*)TIMEXPOS
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          WRITE(*,101)'Timexpos: '//CDUMMY(1:L)
          READ(20)L
          IF(L.GT.0)THEN
            READ(20)OBJECT(1:L)
            WRITE(*,101)'Object  : '//OBJECT(1:TRUELEN(OBJECT))
            WRITE(*,110)'No. of characters: ',L
          ELSE
            OBJECT=CHAR(32)
            WRITE(*,101)'Object  : [not found]'
          END IF
          READ(20)L 
          IF(L.GT.0)THEN
            READ(20)FITSFILE(1:L)
            WRITE(*,101)'FITSfile: '//FITSFILE(1:TRUELEN(FITSFILE))
            WRITE(*,110)'No. of characters: ',L
          ELSE
            FITSFILE=CHAR(32)
            WRITE(*,101)'FITSfile: [not found]'
          END IF
          READ(20)L
          IF(L.GT.0)THEN
            READ(20)COMMENT(1:L)
            WRITE(*,101)'Comment : '//COMMENT(1:TRUELEN(COMMENT))
            WRITE(*,110)'No. of characters: ',L
          ELSE
            COMMENT=CHAR(32)
            WRITE(*,101)'Comment : [not found]'
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
          CLOSE(20)
          STOP
        END IF
        IF(NCHAN.GT.NCMAX)THEN
          WRITE(*,101)'ERROR: NCHAN.GT.NCMAX'
          CLOSE(20)
          STOP
        END IF
C Descomentar siguiente linea y saltar IF para el programa exheads.f
C        GOTO 40
        WRITE(*,100)'Change head information (y/n) '
        CCHANGE(1:1)=READC('n','yn')
        IF(CCHANGE.EQ.'y')THEN
          WRITE(CDUMMY,*)NSCAN
          WRITE(*,100)'NSCAN    '
          NSCAN=READILIM(CDUMMY,1,NSMAX)
c
          WRITE(CDUMMY,*)NCHAN
          WRITE(*,100)'NCHAN    '
          NCHAN=READILIM(CDUMMY,1,NCMAX)
c
          WRITE(CDUMMY,*)STWV
          WRITE(*,100)'STWV     '
          STWV=READF(CDUMMY)
c
          WRITE(CDUMMY,*)DISP
          WRITE(*,100)'DISP     '
          DISP=READF(CDUMMY)
c
          WRITE(CDUMMY,*)AIRMASS
          WRITE(*,100)'AIRMASS  '
          AIRMASS=READF(CDUMMY)
c
          WRITE(CDUMMY,*)TIMEXPOS
          WRITE(*,100)'TIMEXPOS '
          TIMEXPOS=READF(CDUMMY)
c
          L=TRUELEN(OBJECT)
          IF(L.EQ.0)THEN
            WRITE(*,101)'OBJECT  : [not found] '
          ELSE
            WRITE(*,101)'OBJECT  : '//OBJECT(1:L)
          END IF
          WRITE(*,100)'Change OBJECT value (y/n/e=add @ERROR@) '
          CCC(1:1)=READC('n','yne')
          IF(CCC.EQ.'y')THEN
            WRITE(*,100)'OBJECT   '
            OBJECT=READC('@','@')
          ELSEIF(CCC.EQ.'e')THEN
            IF(L.EQ.0)THEN
              OBJECT(1:7)='@ERROR@'
            ELSE
              OBJECT(L+1:L+8)=' @ERROR@'
            END IF
          END IF
c
          L=TRUELEN(FITSFILE)
          IF(L.EQ.0)THEN
            WRITE(*,101)'FITSFILE: [not found] '
          ELSE
            WRITE(*,101)'FITSFILE: '//FITSFILE(1:L)
          END IF
          WRITE(*,100)'Change FITSFILE value (y/n) '
          CCC(1:1)=READC('n','yn')
          IF(CCC.EQ.'y')THEN
            WRITE(*,100)'FITSFILE '
            FITSFILE=READC('@','@')
          END IF
c
          L=TRUELEN(COMMENT)
          IF(L.EQ.0)THEN
            WRITE(*,100)'COMMENT : [not found] '
          ELSE
            WRITE(*,101)'COMMENT : '
            WRITE(*,101)COMMENT(1:L)
          END IF
          WRITE(*,100)'Change COMMENT value (y/n) '
          CCC(1:1)=READC('n','yn')
          IF(CCC.EQ.'y')THEN
            WRITE(*,101)'COMMENT  '
            COMMENT=READC('@','@')
          END IF
c
          DO I=1,NSCAN
            READ(20) (A(J,I),J=1,NCHAN)
          END DO
          WRITE(*,100)'Output file name'
          OUTFILE=OUTFILEX(30,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
          DO I=1,NSCAN
            WRITE(30) (A(J,I),J=1,NCHAN)
          END DO
          CLOSE(30)
        END IF
C
ccc40      CLOSE(20)
        CLOSE(20)
        STOP
C
100     FORMAT(A,$)
101     FORMAT(A)
110     FORMAT(A,I6)
        END
