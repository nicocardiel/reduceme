C------------------------------------------------------------------------------
C Version 01-April-2004                                     file: leefits_new.f
C @ ncl & fjg
C------------------------------------------------------------------------------
Comment
C
C Program: leefits_new
C Classification: input/output
C Description: Reads a FITS file and creates a new file with REDUCEME format.
C
Comment
C
C Este programa lee un fichero en formato FITS normal y lo reescribe en otro
C fichero con el formato para la reduccion con nuestros programas.
C
        PROGRAM LEEFITS_NEW
        IMPLICIT NONE
C
        INTEGER MAXNAXIS
        PARAMETER(MAXNAXIS=3)               !valor maximo admisible para NAXIS
C
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER TRUELEN
        INTEGER READILIM
        REAL READF
C
        INTEGER IUNIT
        INTEGER BITPIX
        INTEGER NAXIS(0:MAXNAXIS)                               !NAXIS=NAXIS(0)
        INTEGER I,J,K
        INTEGER NS1,NS2
        INTEGER ISTATUS,IREADWRITE,NULLVAL
        INTEGER BLOCKSIZE
        INTEGER NHDU,HDUTYPE,HDUNUM
        INTEGER NFOUND,FIRSTPIX
        INTEGER JROW(NCMAX)
        REAL CRPIX1,CRVAL1,CDELT1
        REAL FNULL
        REAL IMAGEN(NCMAX,NSMAX)
        REAL DATAMIN,DATAMAX
        REAL FROW(NCMAX)
        DOUBLE PRECISION DROW(NCMAX)
        CHARACTER*1 CWHOLE
        CHARACTER*50 CDUMMY
        CHARACTER*50 MYCOMMENT
        CHARACTER*80 INFILE,OUTFILE
        LOGICAL LNULL(NCMAX,NSMAX)
        LOGICAL ANYNULL,LANYNULL
        LOGICAL EXTEND
        LOGICAL LOBJECT,LEXPTIME,LAIRMASS
        LOGICAL LCRPIX1,LCRVAL1,LCDELT1
        LOGICAL LROW(NCMAX)
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        THISPROGRAM='leefits_new'
        CALL WELCOME('01-April-2004')
C
        WRITE(*,*)
        WRITE(*,100) 'Input FITS file name'
        INFILE=INFILEX(20,'@',0,0,0.,0.,4,.FALSE.)
        CLOSE(20)                  !lo cerramos porque FITSIO lo vuelve a abrir
C------------------------------------------------------------------------------
C inicializamos variables
        ISTATUS=0               !controla posibles errores durante la ejecucion
        IREADWRITE=0                      !la imagen se abrira en modo READONLY
        NULLVAL=-999
        LANYNULL=.FALSE.
        NHDU=1
C..............................................................................
C localizamos un numero de unidad de fichero no utilizada
        IUNIT=80 !ojo, IUNIT=99 entra en conflicto con el Postcript y es
                 !precisamente este numero el que toma esta funcion
C abrimos el fichero
        CALL FTOPEN(IUNIT,INFILE,IREADWRITE,BLOCKSIZE,ISTATUS)
        CALL FTTHDU(IUNIT,HDUNUM,ISTATUS)
        WRITE(*,*)
        WRITE(*,100) '-> number of HDUs: '
        WRITE(*,*) HDUNUM
C averiguamos si el fichero tiene extensiones
        CALL FTGKYL(IUNIT,'EXTEND',EXTEND,MYCOMMENT,ISTATUS)
        IF(EXTEND)THEN
          WRITE(*,101) '-> INFO: this file contains extensions'
        END IF
        WRITE(*,*)
C..............................................................................
C leemos los keywords buscados en la cabecera del primer HDU (si hay algun
C keyword que no se encuentra, y el fichero contiene extensiones, el programa
C tratara de buscar los keywords que faltan en la cabecera de la extension)
c
        CALL FTGKYS(IUNIT,'OBJECT',OBJECT,MYCOMMENT,ISTATUS)
        IF(ISTATUS.EQ.0)THEN
          LOBJECT=.TRUE.
          WRITE(*,100) '> OBJECT.: '
          WRITE(*,101) OBJECT(1:TRUELEN(OBJECT))
        ELSE
          OBJECT='[not found]'
          LOBJECT=.FALSE.
          ISTATUS=0
        END IF
c
        CALL FTGKYE(IUNIT,'EXPTIME',TIMEXPOS,MYCOMMENT,ISTATUS)
        IF(ISTATUS.EQ.0)THEN
          LEXPTIME=.TRUE.
          WRITE(*,100) '> EXPTIME: '
          WRITE(*,*) TIMEXPOS
        ELSE
          TIMEXPOS=-999
          LEXPTIME=.FALSE.
          ISTATUS=0
        END IF
c
        CALL FTGKYE(IUNIT,'AIRMASS',AIRMASS,MYCOMMENT,ISTATUS)
        IF(ISTATUS.EQ.0)THEN
          LAIRMASS=.TRUE.
          WRITE(*,100) '> AIRMASS: '
          WRITE(*,*) AIRMASS
        ELSE
          AIRMASS=-999
          LAIRMASS=.FALSE.
          ISTATUS=0
        END IF
c
        CALL FTGKYE(IUNIT,'CRPIX1',CRPIX1,MYCOMMENT,ISTATUS)
        IF(ISTATUS.EQ.0)THEN
          LCRPIX1=.TRUE.
          WRITE(*,100) '> CRPIX1: '
          WRITE(*,*) CRPIX1
        ELSE
          CRPIX1=0.0
          LCRPIX1=.FALSE.
          ISTATUS=0
        END IF
c
        CALL FTGKYE(IUNIT,'CRVAL1',CRVAL1,MYCOMMENT,ISTATUS)
        IF(ISTATUS.EQ.0)THEN
          LCRVAL1=.TRUE.
          WRITE(*,100) '> CRVAL1: '
          WRITE(*,*) CRVAL1
        ELSE
          CRVAL1=0.0
          LCRVAL1=.FALSE.
          ISTATUS=0
        END IF
c
        CALL FTGKYE(IUNIT,'CDELT1',CDELT1,MYCOMMENT,ISTATUS)
        IF(ISTATUS.EQ.0)THEN
          LCDELT1=.TRUE.
          WRITE(*,100) '> CDELT1: '
          WRITE(*,*) CDELT1
        ELSE
          CDELT1=0.0
          LCDELT1=.FALSE.
          ISTATUS=0
        END IF
c..............................................................................
C si hay extensiones, saltamos a la primera extension de tipo IMAGE
        IF(EXTEND)THEN
          HDUTYPE=-1
          DO WHILE(HDUTYPE.NE.0)
            NHDU=NHDU+1
            CALL FTMAHD(IUNIT,NHDU,HDUTYPE,ISTATUS)
            IF(ISTATUS.NE.0)THEN
              WRITE(*,101) 'FATAL ERROR while skipping HDU with FTMAHD'
              CALL PRINTERROR(ISTATUS)
              STOP
            END IF
          END DO
c intentamos leer los keywords que no se han encontrado en el primer HDU
          IF(.NOT.LOBJECT)THEN
            CALL FTGKYS(IUNIT,'OBJECT',OBJECT,MYCOMMENT,ISTATUS)
            IF(ISTATUS.EQ.0)THEN
              LOBJECT=.TRUE.
              WRITE(*,100) '> OBJECT.: '
              WRITE(*,101) OBJECT(1:TRUELEN(OBJECT))
            ELSE
              OBJECT='[not found]'
              LOBJECT=.FALSE.
              ISTATUS=0
            END IF
          END IF
c
          IF(.NOT.LEXPTIME)THEN
            CALL FTGKYE(IUNIT,'EXPTIME',TIMEXPOS,MYCOMMENT,ISTATUS)
            IF(ISTATUS.EQ.0)THEN
              LEXPTIME=.TRUE.
              WRITE(*,100) '> EXPTIME: '
              WRITE(*,*) TIMEXPOS
            ELSE
              TIMEXPOS=-999
              LEXPTIME=.FALSE.
              ISTATUS=0
            END IF
          END IF
c
          IF(.NOT.LAIRMASS)THEN
            CALL FTGKYE(IUNIT,'AIRMASS',AIRMASS,MYCOMMENT,ISTATUS)
            IF(ISTATUS.EQ.0)THEN
              LAIRMASS=.TRUE.
              WRITE(*,100) '> AIRMASS: '
              WRITE(*,*) AIRMASS
            ELSE
              AIRMASS=-999
              LAIRMASS=.FALSE.
              ISTATUS=0
            END IF
          END IF
c
          IF(.NOT.LCRPIX1)THEN
            CALL FTGKYE(IUNIT,'CRPIX1',CRPIX1,MYCOMMENT,ISTATUS)
            IF(ISTATUS.EQ.0)THEN
              LCRPIX1=.TRUE.
              WRITE(*,100) '> CRPIX1: '
              WRITE(*,*) CRPIX1
            ELSE
              CRPIX1=0.0
              LCRPIX1=.FALSE.
              ISTATUS=0
            END IF
          END IF
c
          IF(.NOT.LCRVAL1)THEN
            CALL FTGKYE(IUNIT,'CRVAL1',CRVAL1,MYCOMMENT,ISTATUS)
            IF(ISTATUS.EQ.0)THEN
              LCRVAL1=.TRUE.
              WRITE(*,100) '> CRVAL1: '
              WRITE(*,*) CRVAL1
            ELSE
              CRVAL1=0.0
              LCRVAL1=.FALSE.
              ISTATUS=0
            END IF
          END IF
c
          IF(.NOT.LCDELT1)THEN
            CALL FTGKYE(IUNIT,'CDELT1',CDELT1,MYCOMMENT,ISTATUS)
            IF(ISTATUS.EQ.0)THEN
              LCDELT1=.TRUE.
              WRITE(*,100) '> CDELT1: '
              WRITE(*,*) CDELT1
            ELSE
              CDELT1=0.0
              LCDELT1=.FALSE.
              ISTATUS=0
            END IF
          END IF
c
        END IF
c..............................................................................
C leemos BITPIX
        CALL FTGKYJ(IUNIT,'BITPIX',BITPIX,MYCOMMENT,ISTATUS)
        WRITE(*,*)
        WRITE(*,100) '> BITPIX: '
        WRITE(*,*) BITPIX
C leemos NAXIS
        CALL FTGKYJ(IUNIT,'NAXIS',NAXIS(0),MYCOMMENT,ISTATUS)
        WRITE(*,100) '> NAXIS: '
        WRITE(*,*) NAXIS(0)
        IF(NAXIS(0).GT.MAXNAXIS)THEN
          WRITE(*,100)'FATAL ERROR: NAXIS >'
          WRITE(*,*)MAXNAXIS
          CALL FTCLOS(IUNIT,ISTATUS)
          STOP
        END IF
        IF(NAXIS(0).GE.3)THEN
          DO K=3,NAXIS(0)
            IF(NAXIS(K).NE.1)THEN
              WRITE(*,'(A,I1,A)') 'ERROR: NAXIS(',K,') .NE. 1'
              STOP
            END IF
          END DO
        END IF
C leemos NAXIS1... [notar que el quinto parametro es NAXIS(1) en lugar
C de NAXIS para asi recuperar NAXIS(1)...NAXIS(3)]
        CALL FTGKNJ(IUNIT,'NAXIS',1,NAXIS(0),NAXIS(1),NFOUND,ISTATUS)
        IF(NAXIS(0).EQ.1) NAXIS(2)=1
        WRITE(*,100) '> NAXIS1: '
        WRITE(*,*) NAXIS(1)
        WRITE(*,100) '> NAXIS2: '
        WRITE(*,*) NAXIS(2)
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
        WRITE(*,*)
        WRITE(*,100) 'Please wait (reading FITS file)...'
        IF(BITPIX.EQ.16)THEN
          DO I=1,NAXIS(2)
            FIRSTPIX=(I-1)*NAXIS(1)+1          !+(NAXIS3-1)*(NAXIS(1)*NAXIS(2))
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
            FIRSTPIX=(I-1)*NAXIS(1)+1          !+(NAXIS3-1)*(NAXIS(1)*NAXIS(2))
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
            FIRSTPIX=(I-1)*NAXIS(1)+1          !+(NAXIS3-1)*(NAXIS(1)*NAXIS(2))
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
            FIRSTPIX=(I-1)*NAXIS(1)+1          !+(NAXIS3-1)*(NAXIS(1)*NAXIS(2))
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
        IF(LCRPIX1.AND.LCRVAL1.AND.LCDELT1)THEN
          IF(CRPIX1.EQ.1.0)THEN
            STWV=CRVAL1
            DISP=CDELT1
          ELSE
            STWV=CRVAL1+CDELT1*(1.-CRPIX1)
            DISP=CDELT1
          END IF
        ELSEIF(LCRVAL1.AND.LCDELT1)THEN
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
C------------------------------------------------------------------------------
        WRITE(*,100) 'Output file name'
        OUTFILE=OUTFILEX(30,'@',NSCAN,NAXIS(1),STWV,DISP,1,.FALSE.)
        WRITE(*,100) '>>> DATAMIN: '
        WRITE(*,*) DATAMIN
        WRITE(*,100) '>>> DATAMAX: '
        WRITE(*,*) DATAMAX
        IF(ANYNULL)THEN
          WRITE(*,*)
          WRITE(*,101) '*********************************************'
          WRITE(*,101) 'WARNING: this image contains undefined pixels'
          WRITE(*,101) '*********************************************'
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
C------------------------------------------------------------------------------
        STOP
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
