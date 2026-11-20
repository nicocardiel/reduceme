!------------------------------------------------------------------------------
! Version 20-November-2026                                 file:replacefits.f90
! @ ncl & fjg
!------------------------------------------------------------------------------
!Comment
!
! Program: replacefits
! Classification: input/output
! Description: Put data from a REDUCEME file into a previously existing FITS 
! file.
!
! Note: this program requires the FITSIO subroutine package
!
!Comment
!
        PROGRAM REPLACEFITS
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER READILIM
!
        INTEGER SYSTEMFUNCTION
        INTEGER TRUELEN
!
        INTEGER I,J
        INTEGER L1,L2
        INTEGER ISYSTEM
        INTEGER IA(NCMAX,NSMAX)
        REAL A(NCMAX,NSMAX)
        CHARACTER COPC
        CHARACTER*70 COMMENTLINE
        CHARACTER*75 FILENAME,FILECOPY
        CHARACTER*255 INFILE
        CHARACTER*255 COMANDO
        LOGICAL LOGFILE
!
        INTEGER ISTATUS,NAXES(2),NAXIS,BITPIX,BLOCKSIZE,IUNIT
        INTEGER HDUNUM,NHDU,HDUTYPE
        INTEGER NFOUND,GROUP
        INTEGER IREADWRITE
!------------------------------------------------------------------------------
        THISPROGRAM='replacefits'
        CALL WELCOME('20-November-2026')
!------------------------------------------------------------------------------
        OUTFILEX=OUTFILEX !evita un WARNING de compilacion
! leemos imagen FITS original
        LOGFILE=.FALSE.
        DO WHILE(.NOT.LOGFILE)
          WRITE(*,100,ADVANCE='NO')'Input FITS file name....? '
          READ(*,'(A)') FILENAME
          INQUIRE(FILE=FILENAME,EXIST=LOGFILE)
          IF(.NOT.LOGFILE)THEN
            WRITE(*,101,ADVANCE='NO')&
     &         'ERROR: this file does not exist. Try again.'
          END IF
        END DO
! pedimos nombre de imagen FITS de salida (no tiene que existir el fichero)
        LOGFILE=.TRUE.
        DO WHILE(LOGFILE)
          WRITE(*,100,ADVANCE='NO')'Output FITS file name...? '
          READ(*,'(A)') FILECOPY
          INQUIRE(FILE=FILECOPY,EXIST=LOGFILE)
          IF(LOGFILE)THEN
            WRITE(*,101,ADVANCE='NO')&
     &         'ERROR: this file already exist. Try again.'
          END IF
        END DO
!
        L1=TRUELEN(FILENAME)
        L2=TRUELEN(FILECOPY)
! copiamos
        COMANDO='cp '//FILENAME(1:L1)//' '//FILECOPY(1:L2)
        WRITE(*,101,ADVANCE='NO')&
     &     'Executing: '//COMANDO(1:TRUELEN(COMANDO))
        ISYSTEM=SYSTEMFUNCTION(COMANDO)
! damos atributo de escritura al nuevo fichero
        COMANDO='chmod u+w '//FILECOPY(1:L2)
        WRITE(*,101,ADVANCE='NO')&
     &     'Executing: '//COMANDO(1:TRUELEN(COMANDO))
        ISYSTEM=SYSTEMFUNCTION(COMANDO)
!------------------------------------------------------------------------------
! leemos imagen REDUCEME
        WRITE(*,100,ADVANCE='NO')'Input REDUCEME file name'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        DO I=1,NSCAN
          READ(20) (A(J,I),J=1,NCHAN)
        END DO
        CLOSE(20)
!------------------------------------------------------------------------------
        ISTATUS=0
        IREADWRITE=1 !READONLY=0, READWRITE=1
        IUNIT=80
! abrimos el fichero
        CALL FTOPEN(IUNIT,FILECOPY,IREADWRITE,BLOCKSIZE,ISTATUS)
! determinamos numero de extensiones
        CALL FTTHDU(IUNIT,HDUNUM,ISTATUS)
        IF(ISTATUS.NE.0)THEN
          CALL PRINTERROR(ISTATUS)
          STOP
        END IF
! si hay más de una extensión, pedimos número de extensión a modificar
        IF(HDUNUM.GT.1)THEN
          WRITE(*,101,ADVANCE='NO')&
     &       'WARNING: This file contains extensions'
          WRITE(*,100,ADVANCE='NO')'HDU number to be modified '
          NHDU=READILIM('@',1,HDUNUM)
        ELSE
          NHDU=1
        END IF
! establecemos la extensión por defecto
        CALL FTMAHD(IUNIT,NHDU,HDUTYPE,ISTATUS)
! leemos BITPIX
        CALL FTGKYJ(IUNIT,'BITPIX',BITPIX,COMMENTLINE,ISTATUS)
        WRITE(*,*)
        WRITE(*,100,ADVANCE='NO') '>>> BITPIX: '
        WRITE(*,*) BITPIX
! leemos NAXIS y comprobamos que es 1 ó 2
        CALL FTGKYJ(IUNIT,'NAXIS',NAXIS,COMMENTLINE,ISTATUS)
        IF((NAXIS.LT.1).OR.(NAXIS.GT.2))THEN
          WRITE(*,100,ADVANCE='NO') 'NAXIS: '
          WRITE(*,*) NAXIS
          WRITE(*,101,ADVANCE='NO')'FATAL ERROR: invalid NAXIS found.'
          CALL FTCLOS(IUNIT,ISTATUS)
          STOP
        ELSEIF(NAXIS.EQ.1)THEN
          CALL FTGKNJ(IUNIT,'NAXIS',1,1,NAXES,NFOUND,ISTATUS)
          NAXES(2)=1
        ELSEIF(NAXIS.EQ.2)THEN
          CALL FTGKNJ(IUNIT,'NAXIS',1,2,NAXES,NFOUND,ISTATUS)
        ELSE
          CALL FTCLOS(IUNIT,ISTATUS)
          STOP 'FATAL ERROR: this line should never be executed.'
        END IF
        WRITE(*,*)
        WRITE(*,100,ADVANCE='NO') '>>> NAXIS2: '
        WRITE(*,*) NAXES(2)
        WRITE(*,100,ADVANCE='NO') '>>> NAXIS1: '
        WRITE(*,*) NAXES(1)
! comprobamos que las dimensiones coinciden
        IF((NAXES(1).NE.NCHAN).AND.(NAXES(2).NE.NSCAN))THEN
          WRITE(*,100,ADVANCE='NO')'NAXIS1,NCHAN: '
          WRITE(*,*) NAXES(1),NCHAN
          WRITE(*,100,ADVANCE='NO')'NAXIS2,NSCAN: '
          WRITE(*,*) NAXES(2),NSCAN
          WRITE(*,101,ADVANCE='NO')&
     &       'FATAL ERROR: dimensions do not match.'
          CALL FTCLOS(IUNIT,ISTATUS)
          STOP
        END IF
! establecemos el nuevo BITPIX
        IF(BITPIX.EQ.16)THEN
          WRITE(*,101,ADVANCE='NO') '* Select BITPIX for output image:'
          WRITE(*,101,ADVANCE='NO') '(1) BITPIX=16'
          WRITE(*,101,ADVANCE='NO') '(2) BITPIX=-32'
          WRITE(*,100,ADVANCE='NO') 'Option (1/2) '
          COPC(1:1)=READC('@','12')
          IF(COPC.EQ.'1')THEN
            !vamos a respetar BITPIX 16
            CALL FTRSIM(IUNIT,16,2,NAXES,ISTATUS)
          ELSE
            !pasamos a -32
            CALL FTRSIM(IUNIT,-32,2,NAXES,ISTATUS)
          END IF
        ELSE
          IF(BITPIX.NE.-32)THEN
            WRITE(*,101,ADVANCE='NO')&
     &         'WARNING: BITPIX will be changed to -32'
            CALL FTRSIM(IUNIT,-32,2,NAXES,ISTATUS)
          END IF
        END IF
! escribimos los datos
        GROUP=1
        IF((BITPIX.EQ.16).AND.(COPC.EQ.'1'))THEN
          WRITE(*,101,ADVANCE='NO')&
     &       'WARNING: output data will be rounded'//&
     &     ' to INTEGER'
          DO I=1,NSCAN
            DO J=1,NCHAN
              IA(J,I)=NINT(A(J,I))
            END DO
          END DO
          CALL FTP2DJ(IUNIT,GROUP,NCMAX,NCHAN,NSCAN,IA,ISTATUS)
        ELSE
          CALL FTP2DE(IUNIT,GROUP,NCMAX,NCHAN,NSCAN,A,ISTATUS)
        END IF
! cerramos el fichero FITS
        CALL FTCLOS(IUNIT,ISTATUS)
        IF(ISTATUS.GT.0) CALL PRINTERROR(ISTATUS)
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
