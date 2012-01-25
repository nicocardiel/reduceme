C------------------------------------------------------------------------------
C Version 05-December-2007                                   file:replacefits.f
C @ ncl & fjg
C------------------------------------------------------------------------------
Comment
C
C Program: replacefits
C Classification: input/output
C Description: Put data from a REDUCEME file into a previously existing FITS 
C file.
C
C Note: this program requires the FITSIO subroutine package
C
Comment
C
        PROGRAM REPLACEFITS
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER READILIM
C
        INTEGER SYSTEMFUNCTION
        INTEGER TRUELEN
C
        INTEGER I,J
        INTEGER L1,L2
        INTEGER ISYSTEM
        INTEGER IA(NCMAX,NSMAX)
        REAL A(NCMAX,NSMAX)
        CHARACTER COPC
        CHARACTER*70 COMMENTLINE
        CHARACTER*75 FILENAME,FILECOPY,INFILE
        CHARACTER*255 COMANDO
        LOGICAL LOGFILE
C
        INTEGER ISTATUS,NAXES(2),NAXIS,BITPIX,BLOCKSIZE,IUNIT
        INTEGER HDUNUM,NHDU,HDUTYPE
        INTEGER NFOUND,GROUP
        INTEGER IREADWRITE
C------------------------------------------------------------------------------
        THISPROGRAM='replacefits'
        CALL WELCOME('05-December-2007')
C------------------------------------------------------------------------------
        OUTFILEX=OUTFILEX !evita un WARNING de compilacion
C leemos imagen FITS original
        LOGFILE=.FALSE.
        DO WHILE(.NOT.LOGFILE)
          WRITE(*,100)'Input FITS file name....? '
          READ(*,'(A)') FILENAME
          INQUIRE(FILE=FILENAME,EXIST=LOGFILE)
          IF(.NOT.LOGFILE)THEN
            WRITE(*,101)'ERROR: this file does not exist. Try again.'
          END IF
        END DO
C pedimos nombre de imagen FITS de salida (no tiene que existir el fichero)
        LOGFILE=.TRUE.
        DO WHILE(LOGFILE)
          WRITE(*,100)'Output FITS file name...? '
          READ(*,'(A)') FILECOPY
          INQUIRE(FILE=FILECOPY,EXIST=LOGFILE)
          IF(LOGFILE)THEN
            WRITE(*,101)'ERROR: this file already exist. Try again.'
          END IF
        END DO
C
        L1=TRUELEN(FILENAME)
        L2=TRUELEN(FILECOPY)
C copiamos
        COMANDO='cp '//FILENAME(1:L1)//' '//FILECOPY(1:L2)
        WRITE(*,101) 'Executing: '//COMANDO(1:TRUELEN(COMANDO))
        ISYSTEM=SYSTEMFUNCTION(COMANDO)
C damos atributo de escritura al nuevo fichero
        COMANDO='chmod u+w '//FILECOPY(1:L2)
        WRITE(*,101) 'Executing: '//COMANDO(1:TRUELEN(COMANDO))
        ISYSTEM=SYSTEMFUNCTION(COMANDO)
C------------------------------------------------------------------------------
C leemos imagen REDUCEME
        WRITE(*,100)'Input REDUCEME file name'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        DO I=1,NSCAN
          READ(20) (A(J,I),J=1,NCHAN)
        END DO
        CLOSE(20)
C------------------------------------------------------------------------------
        ISTATUS=0
        IREADWRITE=1 !READONLY=0, READWRITE=1
        IUNIT=80
C abrimos el fichero
        CALL FTOPEN(IUNIT,FILECOPY,IREADWRITE,BLOCKSIZE,ISTATUS)
C determinamos numero de extensiones
        CALL FTTHDU(IUNIT,HDUNUM,ISTATUS)
        IF(ISTATUS.NE.0)THEN
          CALL PRINTERROR(ISTATUS)
          STOP
        END IF
C si hay más de una extensión, pedimos número de extensión a modificar
        IF(HDUNUM.GT.1)THEN
          WRITE(*,101)'WARNING: This file contains extensions'
          WRITE(*,100)'HDU number to be modified '
          NHDU=READILIM('@',1,HDUNUM)
        ELSE
          NHDU=1
        END IF
C establecemos la extensión por defecto
        CALL FTMAHD(IUNIT,NHDU,HDUTYPE,ISTATUS)
C leemos BITPIX
        CALL FTGKYJ(IUNIT,'BITPIX',BITPIX,COMMENTLINE,ISTATUS)
        WRITE(*,*)
        WRITE(*,100) '>>> BITPIX: '
        WRITE(*,*) BITPIX
C leemos NAXIS y comprobamos que es 1 ó 2
        CALL FTGKYJ(IUNIT,'NAXIS',NAXIS,COMMENTLINE,ISTATUS)
        IF((NAXIS.LT.1).OR.(NAXIS.GT.2))THEN
          WRITE(*,100) 'NAXIS: '
          WRITE(*,*) NAXIS
          WRITE(*,101)'FATAL ERROR: invalid NAXIS found.'
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
        WRITE(*,100) '>>> NAXIS2: '
        WRITE(*,*) NAXES(2)
        WRITE(*,100) '>>> NAXIS1: '
        WRITE(*,*) NAXES(1)
C comprobamos que las dimensiones coinciden
        IF((NAXES(1).NE.NCHAN).AND.(NAXES(2).NE.NSCAN))THEN
          WRITE(*,100)'NAXIS1,NCHAN: '
          WRITE(*,*) NAXES(1),NCHAN
          WRITE(*,100)'NAXIS2,NSCAN: '
          WRITE(*,*) NAXES(2),NSCAN
          WRITE(*,101)'FATAL ERROR: dimensions do not match.'
          CALL FTCLOS(IUNIT,ISTATUS)
          STOP
        END IF
C establecemos el nuevo BITPIX
        IF(BITPIX.EQ.16)THEN
          WRITE(*,101) '* Select BITPIX for output image:'
          WRITE(*,101) '(1) BITPIX=16'
          WRITE(*,101) '(2) BITPIX=-32'
          WRITE(*,100) 'Option (1/2) '
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
            WRITE(*,101)'WARNING: BITPIX will be changed to -32'
            CALL FTRSIM(IUNIT,-32,2,NAXES,ISTATUS)
          END IF
        END IF
C escribimos los datos
        GROUP=1
        IF((BITPIX.EQ.16).AND.(COPC.EQ.'1'))THEN
          WRITE(*,101)'WARNING: output data will be rounded'//
     +     ' to INTEGER'
          DO I=1,NSCAN
            DO J=1,NCHAN
              IA(J,I)=NINT(A(J,I))
            END DO
          END DO
          CALL FTP2DJ(IUNIT,GROUP,NCMAX,NCHAN,NSCAN,IA,ISTATUS)
        ELSE
          CALL FTP2DE(IUNIT,GROUP,NCMAX,NCHAN,NSCAN,A,ISTATUS)
        END IF
C cerramos el fichero FITS
        CALL FTCLOS(IUNIT,ISTATUS)
        IF(ISTATUS.GT.0) CALL PRINTERROR(ISTATUS)
C------------------------------------------------------------------------------
        STOP
100     FORMAT(A,$)
101     FORMAT(A,$)
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
