!------------------------------------------------------------------------------
! Version 24-November-2026                                file: leefits_new.f90
! @ ncl & fjg
!------------------------------------------------------------------------------
!Comment
!
! Program: leefits_new
! Classification: input/output
! Description: Reads a FITS file and creates a new file with REDUCEME format.
!
!Comment
!
! Este programa lee un fichero en formato FITS normal y lo reescribe en otro
! fichero con el formato para la reduccion con nuestros programas.
!
        PROGRAM LEEFITS_NEW
        USE RED_M
        USE FUTILS_M
        USE REDLIB_M
        IMPLICIT NONE
!
        INTEGER MAXNAXIS
        PARAMETER(MAXNAXIS=3)               !valor maximo admisible para NAXIS
!
        INCLUDE 'redlib.inc'
!
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
        REAL, ALLOCATABLE :: IMAGEN(:,:)
        REAL DATAMIN,DATAMAX
        REAL FROW(NCMAX)
        DOUBLE PRECISION DROW(NCMAX)
        CHARACTER*1 CWHOLE
        CHARACTER*50 CDUMMY
        CHARACTER*50 MYCOMMENT
        CHARACTER*255 INFILE,OUTFILE
        LOGICAL, ALLOCATABLE :: LNULL(:,:)
        LOGICAL ANYNULL,LANYNULL
        LOGICAL EXTEND
        LOGICAL LOBJECT,LEXPTIME,LAIRMASS
        LOGICAL LCRPIX1,LCRVAL1,LCDELT1
        LOGICAL LROW(NCMAX)
!------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
!
        THISPROGRAM='leefits_new'
        CALL WELCOME('24-November-2026')
!
        WRITE(*,*)
        WRITE(*,100,ADVANCE='NO') 'Input FITS file name'
        INFILE=INFILEX(20,'@',0,0,0.,0.,4,.FALSE.)
        CLOSE(20)                  !lo cerramos porque FITSIO lo vuelve a abrir
!------------------------------------------------------------------------------
! inicializamos variables
        ISTATUS=0               !controla posibles errores durante la ejecucion
        IREADWRITE=0                      !la imagen se abrira en modo READONLY
        NULLVAL=-999
        LANYNULL=.FALSE.
        NHDU=1
!..............................................................................
! localizamos un numero de unidad de fichero no utilizada
        IUNIT=80 !ojo, IUNIT=99 entra en conflicto con el Postcript y es
                 !precisamente este numero el que toma esta funcion
! abrimos el fichero
        CALL FTOPEN(IUNIT,INFILE,IREADWRITE,BLOCKSIZE,ISTATUS)
        CALL FTTHDU(IUNIT,HDUNUM,ISTATUS)
        WRITE(*,*)
        WRITE(*,100,ADVANCE='NO') '-> number of HDUs: '
        WRITE(*,*) HDUNUM
! averiguamos si el fichero tiene extensiones
        CALL FTGKYL(IUNIT,'EXTEND',EXTEND,MYCOMMENT,ISTATUS)
! keyword EXTEND not found (KEY_NO_EXIST): no extensions; reset the status,
! otherwise the next FITSIO calls do nothing
        IF(ISTATUS.EQ.202)THEN
          EXTEND=.FALSE.
          ISTATUS=0
        END IF
        IF(EXTEND)THEN
          WRITE(*,101) '-> INFO: this file contains extensions'
        END IF
        WRITE(*,*)
!..............................................................................
! leemos los keywords buscados en la cabecera del primer HDU (si hay algun
! keyword que no se encuentra, y el fichero contiene extensiones, el programa
! tratara de buscar los keywords que faltan en la cabecera de la extension)
!
        CALL FTGKYS(IUNIT,'OBJECT',OBJECT,MYCOMMENT,ISTATUS)
        IF(ISTATUS.EQ.0)THEN
          LOBJECT=.TRUE.
          WRITE(*,100,ADVANCE='NO') '> OBJECT.: '
          WRITE(*,101) OBJECT(1:TRUELEN(OBJECT))
        ELSE
          OBJECT='[not found]'
          LOBJECT=.FALSE.
          ISTATUS=0
        END IF
!
        CALL FTGKYE(IUNIT,'EXPTIME',TIMEXPOS,MYCOMMENT,ISTATUS)
        IF(ISTATUS.EQ.0)THEN
          LEXPTIME=.TRUE.
          WRITE(*,100,ADVANCE='NO') '> EXPTIME: '
          WRITE(*,*) TIMEXPOS
        ELSE
          TIMEXPOS=-999
          LEXPTIME=.FALSE.
          ISTATUS=0
        END IF
!
        CALL FTGKYE(IUNIT,'AIRMASS',AIRMASS,MYCOMMENT,ISTATUS)
        IF(ISTATUS.EQ.0)THEN
          LAIRMASS=.TRUE.
          WRITE(*,100,ADVANCE='NO') '> AIRMASS: '
          WRITE(*,*) AIRMASS
        ELSE
          AIRMASS=-999
          LAIRMASS=.FALSE.
          ISTATUS=0
        END IF
!
        CALL FTGKYE(IUNIT,'CRPIX1',CRPIX1,MYCOMMENT,ISTATUS)
        IF(ISTATUS.EQ.0)THEN
          LCRPIX1=.TRUE.
          WRITE(*,100,ADVANCE='NO') '> CRPIX1: '
          WRITE(*,*) CRPIX1
        ELSE
          CRPIX1=0.0
          LCRPIX1=.FALSE.
          ISTATUS=0
        END IF
!
        CALL FTGKYE(IUNIT,'CRVAL1',CRVAL1,MYCOMMENT,ISTATUS)
        IF(ISTATUS.EQ.0)THEN
          LCRVAL1=.TRUE.
          WRITE(*,100,ADVANCE='NO') '> CRVAL1: '
          WRITE(*,*) CRVAL1
        ELSE
          CRVAL1=0.0
          LCRVAL1=.FALSE.
          ISTATUS=0
        END IF
!
        CALL FTGKYE(IUNIT,'CDELT1',CDELT1,MYCOMMENT,ISTATUS)
        IF(ISTATUS.EQ.0)THEN
          LCDELT1=.TRUE.
          WRITE(*,100,ADVANCE='NO') '> CDELT1: '
          WRITE(*,*) CDELT1
        ELSE
          CDELT1=0.0
          LCDELT1=.FALSE.
          ISTATUS=0
        END IF
!..............................................................................
! si hay extensiones, saltamos a la primera extension de tipo IMAGE
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
! intentamos leer los keywords que no se han encontrado en el primer HDU
          IF(.NOT.LOBJECT)THEN
            CALL FTGKYS(IUNIT,'OBJECT',OBJECT,MYCOMMENT,ISTATUS)
            IF(ISTATUS.EQ.0)THEN
              LOBJECT=.TRUE.
              WRITE(*,100,ADVANCE='NO') '> OBJECT.: '
              WRITE(*,101) OBJECT(1:TRUELEN(OBJECT))
            ELSE
              OBJECT='[not found]'
              LOBJECT=.FALSE.
              ISTATUS=0
            END IF
          END IF
!
          IF(.NOT.LEXPTIME)THEN
            CALL FTGKYE(IUNIT,'EXPTIME',TIMEXPOS,MYCOMMENT,ISTATUS)
            IF(ISTATUS.EQ.0)THEN
              LEXPTIME=.TRUE.
              WRITE(*,100,ADVANCE='NO') '> EXPTIME: '
              WRITE(*,*) TIMEXPOS
            ELSE
              TIMEXPOS=-999
              LEXPTIME=.FALSE.
              ISTATUS=0
            END IF
          END IF
!
          IF(.NOT.LAIRMASS)THEN
            CALL FTGKYE(IUNIT,'AIRMASS',AIRMASS,MYCOMMENT,ISTATUS)
            IF(ISTATUS.EQ.0)THEN
              LAIRMASS=.TRUE.
              WRITE(*,100,ADVANCE='NO') '> AIRMASS: '
              WRITE(*,*) AIRMASS
            ELSE
              AIRMASS=-999
              LAIRMASS=.FALSE.
              ISTATUS=0
            END IF
          END IF
!
          IF(.NOT.LCRPIX1)THEN
            CALL FTGKYE(IUNIT,'CRPIX1',CRPIX1,MYCOMMENT,ISTATUS)
            IF(ISTATUS.EQ.0)THEN
              LCRPIX1=.TRUE.
              WRITE(*,100,ADVANCE='NO') '> CRPIX1: '
              WRITE(*,*) CRPIX1
            ELSE
              CRPIX1=0.0
              LCRPIX1=.FALSE.
              ISTATUS=0
            END IF
          END IF
!
          IF(.NOT.LCRVAL1)THEN
            CALL FTGKYE(IUNIT,'CRVAL1',CRVAL1,MYCOMMENT,ISTATUS)
            IF(ISTATUS.EQ.0)THEN
              LCRVAL1=.TRUE.
              WRITE(*,100,ADVANCE='NO') '> CRVAL1: '
              WRITE(*,*) CRVAL1
            ELSE
              CRVAL1=0.0
              LCRVAL1=.FALSE.
              ISTATUS=0
            END IF
          END IF
!
          IF(.NOT.LCDELT1)THEN
            CALL FTGKYE(IUNIT,'CDELT1',CDELT1,MYCOMMENT,ISTATUS)
            IF(ISTATUS.EQ.0)THEN
              LCDELT1=.TRUE.
              WRITE(*,100,ADVANCE='NO') '> CDELT1: '
              WRITE(*,*) CDELT1
            ELSE
              CDELT1=0.0
              LCDELT1=.FALSE.
              ISTATUS=0
            END IF
          END IF
!
        END IF
!..............................................................................
! leemos BITPIX
        CALL FTGKYJ(IUNIT,'BITPIX',BITPIX,MYCOMMENT,ISTATUS)
        WRITE(*,*)
        WRITE(*,100,ADVANCE='NO') '> BITPIX: '
        WRITE(*,*) BITPIX
! leemos NAXIS
        CALL FTGKYJ(IUNIT,'NAXIS',NAXIS(0),MYCOMMENT,ISTATUS)
        WRITE(*,100,ADVANCE='NO') '> NAXIS: '
        WRITE(*,*) NAXIS(0)
        IF(NAXIS(0).GT.MAXNAXIS)THEN
          WRITE(*,100,ADVANCE='NO')'FATAL ERROR: NAXIS >'
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
! leemos NAXIS1... [notar que el quinto parametro es NAXIS(1) en lugar
! de NAXIS para asi recuperar NAXIS(1)...NAXIS(3)]
        CALL FTGKNJ(IUNIT,'NAXIS',1,NAXIS(0),NAXIS(1),NFOUND,ISTATUS)
        IF(NAXIS(0).EQ.1) NAXIS(2)=1
        WRITE(*,100,ADVANCE='NO') '> NAXIS1: '
        WRITE(*,*) NAXIS(1)
        WRITE(*,100,ADVANCE='NO') '> NAXIS2: '
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
! leemos la imagen
        WRITE(*,*)
        ! allocate with the actual image dimensions (static arrays
        ! were limited by NCMAX,NSMAX); initialised as before
        IF(ALLOCATED(IMAGEN)) DEALLOCATE(IMAGEN)
        IF(ALLOCATED(LNULL)) DEALLOCATE(LNULL)
        ALLOCATE(IMAGEN(NAXIS(1),NAXIS(2)),LNULL(NAXIS(1),NAXIS(2)))
        IMAGEN=0.0
        LNULL=.FALSE.
        WRITE(*,100,ADVANCE='NO') 'Please wait (reading FITS file)...'
        IF(BITPIX.EQ.16)THEN
          DO I=1,NAXIS(2)
            FIRSTPIX=(I-1)*NAXIS(1)+1          !+(NAXIS3-1)*(NAXIS(1)*NAXIS(2))
            CALL FTGPFJ(IUNIT,1,FIRSTPIX,NAXIS(1),JROW(1),LROW(1),&
     &       ANYNULL,ISTATUS)
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
            CALL FTGPFE(IUNIT,1,FIRSTPIX,NAXIS(1),FROW(1),LROW(1),&
     &       ANYNULL,ISTATUS)
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
!cc       CALL FTG2DE(IUNIT,1,NULLVAL,NCMAX,NAXIS(1),NAXIS(2),
!cc     +     IMAGEN,ANYNULL,ISTATUS)
        ELSEIF(BITPIX.EQ.-32)THEN
          DO I=1,NAXIS(2)
            FIRSTPIX=(I-1)*NAXIS(1)+1          !+(NAXIS3-1)*(NAXIS(1)*NAXIS(2))
            CALL FTGPFE(IUNIT,1,FIRSTPIX,NAXIS(1),FROW(1),LROW(1),&
     &       ANYNULL,ISTATUS)
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
!cc       CALL FTG2DE(IUNIT,1,NULLVAL,NCMAX,NAXIS(1),NAXIS(2),
!cc     +     IMAGEN,ANYNULL,ISTATUS)
        ELSEIF(BITPIX.EQ.-64)THEN
          DO I=1,NAXIS(2)
            FIRSTPIX=(I-1)*NAXIS(1)+1          !+(NAXIS3-1)*(NAXIS(1)*NAXIS(2))
            CALL FTGPFD(IUNIT,1,FIRSTPIX,NAXIS(1),DROW(1),LROW(1),&
     &       ANYNULL,ISTATUS)
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
          WRITE(*,100,ADVANCE='NO')&
     &       'FATAL ERROR in subroutine SLEEFITS: BITPIX ='
          WRITE(*,*) BITPIX
          CALL FTCLOS(IUNIT,ISTATUS)
          STOP
        END IF
        WRITE(*,101) '  ..OK!'
! cerramos el fichero
        CALL FTCLOS(IUNIT,ISTATUS)
! liberamos el numero de unidad del fichero utilizado
        CALL FTFIOU(IUNIT,ISTATUS)
! calculamos maximo y minimo
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
! chequeamos si se ha producido algun error
        IF(ISTATUS.GT.0)THEN
          CALL PRINTERROR(ISTATUS)
        END IF
        ANYNULL=LANYNULL                  !basta que haya ocurrido una sola vez
!------------------------------------------------------------------------------
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
        WRITE(*,100,ADVANCE='NO') '> STWV: '
        WRITE(*,*) STWV
        WRITE(*,100,ADVANCE='NO') '> DISP: '
        WRITE(*,*) DISP
        WRITE(*,*)
!
        IF(NAXIS(2).GT.1)THEN
          WRITE(*,100,ADVANCE='NO') 'Save whole frame (y/n) '
          CWHOLE(1:1)=READC('y','yn')
          IF(CWHOLE.EQ.'n')THEN
            WRITE(*,100,ADVANCE='NO') 'First scan: '
            NS1=READILIM('1',1,NAXIS(2))
            WRITE(*,100,ADVANCE='NO') 'Last  scan: '
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
!------------------------------------------------------------------------------
        WRITE(*,100,ADVANCE='NO') 'Output file name'
        OUTFILE=OUTFILEX(30,'@',NSCAN,NAXIS(1),STWV,DISP,1,.FALSE.)
        WRITE(*,100,ADVANCE='NO') '>>> DATAMIN: '
        WRITE(*,*) DATAMIN
        WRITE(*,100,ADVANCE='NO') '>>> DATAMAX: '
        WRITE(*,*) DATAMAX
        IF(ANYNULL)THEN
          WRITE(*,*)
          WRITE(*,101) '*********************************************'
          WRITE(*,101) 'WARNING: this image contains undefined pixels'
          WRITE(*,101) '*********************************************'
          WRITE(*,*)
          WRITE(*,100,ADVANCE='NO') 'Real value for undefined pixels '
          FNULL=READF('@')
          DO I=NS1,NS2
            DO J=1,NAXIS(1)
              IF(LNULL(J,I)) IMAGEN(J,I)=FNULL
            END DO
          END DO
        END IF
        WRITE(*,100,ADVANCE='NO') 'Please wait (saving file)...'
        DO I=NS1,NS2
          WRITE(30) (IMAGEN(J,I),J=1,NAXIS(1))
        END DO
        WRITE(*,101) '  ..OK!'
        CALL CLOSEFILEX(30)
!------------------------------------------------------------------------------
        STOP
100     FORMAT(A)
101     FORMAT(A)
        END
!
!******************************************************************************
!
        SUBROUTINE PRINTERROR(ISTATUS)
! Print out the FITSIO error messages to the user
        INTEGER ISTATUS
        CHARACTER ERRTEXT*30,ERRMESSAGE*80
! Check if status is OK (no error); if so, simply return
        IF(ISTATUS.LE.0) RETURN
! Get the text string which describes the error
        CALL FTGERR(ISTATUS,ERRTEXT)
        WRITE(*,'(A)',ADVANCE='NO')'FITSIO Error Status = '
        WRITE(*,*)ISTATUS
        WRITE(*,'(A)')ERRTEXT
! Read and print out all the error messages on the FITSIO stack
        CALL FTGMSG(ERRMESSAGE)
        DO WHILE(ERRMESSAGE.NE.' ')
          WRITE(*,'(A)') ERRMESSAGE
          CALL FTGMSG(ERRMESSAGE)
        END DO
        END
