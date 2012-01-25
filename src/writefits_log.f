C------------------------------------------------------------------------------
C Version 14-June-2003                                     file:writefits_log.f
C @ ncl & fjg
C------------------------------------------------------------------------------
Comment
C
C Program: writefits_log
C Classification: input/output
C Description: Reads a file with REDUCEME format and creates a new file with
C FITS format and CTYPE1="WAVE-LOG".
C
C Note: this program requires the FITSIO subroutine package
C
Comment
C
        PROGRAM WRITEFITS_LOG
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER READI
C
        INTEGER I,J
        REAL A(NCMAX,NSMAX)
        CHARACTER*1 CHEADER
        CHARACTER*80 CFECHA,RECORD
        CHARACTER*70 COMMENTLINE
        CHARACTER*75 FILENAME,OUTFILE,TEMPFILE
        LOGICAL LOGFILE
C
        INTEGER ISTATUS,NAXES(2),NAXIS,BITPIX,BLOCKSIZE,UNIT,IUNIT
        INTEGER READWRITE,NKEYS,NSPACE,NSKIP
        INTEGER GROUP,IDAY,IMONTH,IYEAR
        LOGICAL SIMPLE,EXTEND
C------------------------------------------------------------------------------
        OUTFILEX=OUTFILEX
C
        THISPROGRAM='writefits_log'
        CALL WELCOME('14-June-2003')
C
        WRITE(*,100)'Input file name'
        FILENAME=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        DO I=1,NSCAN
          READ(20) (A(J,I),J=1,NCHAN)
        END DO
        CLOSE(20)
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
        ISTATUS=0
C salvamos fichero
10      WRITE(*,100)'Output file name? '
        READ(*,'(A)')OUTFILE
        INQUIRE(FILE=OUTFILE,EXIST=LOGFILE)
        IF(LOGFILE)THEN
          WRITE(*,101)'ERROR: this file already exist. Try again.'
          GOTO 10
        END IF
C get an unused Logical Unit Number to use to open the FITS file
        CALL FTGIOU(UNIT,ISTATUS)
        CALL FTGIOU(IUNIT,ISTATUS)
C create the new empty FITS file
        BLOCKSIZE=1
        CALL FTINIT(UNIT,OUTFILE,BLOCKSIZE,ISTATUS)
C initialize parameters about the FITS image
        SIMPLE=.TRUE.
        BITPIX=-32
        NAXIS=2
        NAXES(1)=NCHAN
        NAXES(2)=NSCAN
        EXTEND=.FALSE.
C write the required header keywords
        CALL FTPHPR(UNIT,SIMPLE,BITPIX,NAXIS,NAXES,0,1,EXTEND,ISTATUS)
        WRITE(*,100) 'Are you including the FITS header from '
        WRITE(*,100) 'another file (y/n) '
        CHEADER(1:1)=READC('n','yn')
        IF(CHEADER.EQ.'n')THEN
C write REDUCEME required keywords
          DO I=1,70
            COMMENTLINE(I:I)='-'
          END DO
          CALL FTPCOM(UNIT,COMMENTLINE,ISTATUS)
          CALL FTPCOM(UNIT,'************************** '//
     +     'REDUCEME HEADER '//
     +     '***************************',ISTATUS)
          CALL FTPCOM(UNIT,COMMENTLINE,ISTATUS)
          CALL FTGSDT(IDAY,IMONTH,IYEAR,ISTATUS)
          WRITE(CFECHA,'(I2.2,A1,I2.2,A1,I2.2)')IDAY,'/',IMONTH,'/',
     +     IYEAR
          CALL FTPHIS(UNIT,'Date: '//CFECHA,ISTATUS)
          CALL FTPKYS(UNIT,'CTYPE1','WAVE-LOG',' ',ISTATUS)
          CALL FTPKYS(UNIT,'CUNIT1','Angstrom',' ',ISTATUS)
          CALL FTPKYF(UNIT,'CRPIX1',1.0,1,' ',ISTATUS)
          CALL FTPKYF(UNIT,'CRVAL1',STWV,4,
     +     'central wavelength of first pixel',ISTATUS)
          CALL FTPKYF(UNIT,'CDELT1',DISP,6,
     +     'linear dispersion (Angstrom/pixel)',ISTATUS)
          CALL FTPKYS(UNIT,'OBJECT',OBJECT,'Object name',ISTATUS)
          CALL FTPKYS(UNIT,'FITSfi',FITSFILE,'Object name',ISTATUS)
          CALL FTPCOM(UNIT,COMMENT,ISTATUS)
          CALL FTPKYF(UNIT,'AIRMASS',AIRMASS,5,'Airmass',ISTATUS)
          CALL FTPKYF(UNIT,'TIMEXPO',TIMEXPOS,1,'Timexpos',ISTATUS)
        ELSE
C read FITS header from external FITS file
20        WRITE(*,100) 'Header template file name? '
          READ(*,'(A)') TEMPFILE
          INQUIRE(FILE=TEMPFILE,EXIST=LOGFILE)
          IF(.NOT.LOGFILE)THEN
            WRITE(*,101) 'ERROR: this file does not exist. Try again.'
            GOTO 20
          END IF
          READWRITE=0 !solo lectura
          CALL FTOPEN(IUNIT,TEMPFILE,READWRITE,BLOCKSIZE,ISTATUS)
          CALL FTGHSP(IUNIT,NKEYS,NSPACE,ISTATUS)
          WRITE(*,100) 'KEYS='
          WRITE(*,*) NKEYS
          WRITE(*,101) 'NOTE: the first keywords of the header template'
          WRITE(*,101) 'will be skipped since they have already been'
          WRITE(*,101) 'written by this program. In particular:'
          WRITE(*,101) 'SIMPLE  =                    T'
          WRITE(*,101) 'BITPIX  =                  -32'
          WRITE(*,101) 'NAXIS   =                    2'
          WRITE(*,101) 'NAXIS1  =                 ????'
          WRITE(*,101) 'NAXIS2  =                 ????'
          WRITE(*,101) ' '
          WRITE(*,101) 'So, you must know the actual number of'
          WRITE(*,101) 'keywords including the NAXIS2 description.'
          WRITE(*,100) 'No. of initial keywords to skip '
          NSKIP=READI('5')
          DO I=1,NSKIP
            CALL FTGREC(IUNIT,I,RECORD,ISTATUS)
          END DO
          DO I=NSKIP+1,NKEYS
            CALL FTGREC(IUNIT,I,RECORD,ISTATUS)
            CALL FTPREC(UNIT,RECORD,ISTATUS) 
          END DO
          CALL FTCLOS(IUNIT,ISTATUS)
          CALL FTFIOU(IUNIT,ISTATUS)
        END IF
C write the array to the FITS file
        GROUP=1
        CALL FTP2DE(UNIT,GROUP,NCMAX,NCHAN,NSCAN,A,ISTATUS)
C close the file and free the unit number
        CALL FTCLOS(UNIT,ISTATUS)
        CALL FTFIOU(UNIT,ISTATUS)
        IF(ISTATUS.GT.0) CALL PRINTERROR(ISTATUS)
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
