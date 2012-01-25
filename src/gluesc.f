C------------------------------------------------------------------------------
C Version 6-December-1996                                        file: gluesc.f
C @ ncl & fjg
C------------------------------------------------------------------------------
Comment
C
C Program: gluesc
C Classification: arithmetic & manipulations
C Description: Takes different spectra and creates a new frame (in which each 
C individual spectrum correspond to the former single spectra).
C
Comment
C
C "Pega" varios espectros para formar una unica imagen
C
        PROGRAM GLUESC
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER READILIM
C
        INTEGER I,J,K,II
        INTEGER NCHAN2,NSKIP
        REAL S(NCMAX),A(NCMAX,NSMAX),ERR(NCMAX,NSMAX)
        REAL STWV2,DISP2
        CHARACTER*1 COUTHEADER,CERR
        CHARACTER*75 FILENAME,ERRFILE,OUTFILE,HEADFILE
C------------------------------------------------------------------------------
        THISPROGRAM='gluesc'
        CALL WELCOME('6-December-1996')
C
        WRITE(*,100)'Work with error images (y/n) '
        CERR(1:1)=READC('n','yn')
C
        WRITE(*,100)'No. of spectra '
        K=READILIM('@',1,NSMAX)
C
        WRITE(*,100)'Create output file with header description (y/n) '
        COUTHEADER(1:1)=READC('n','yn')
        IF(COUTHEADER.EQ.'y')THEN
          WRITE(*,100)'Header description file name'
          HEADFILE=OUTFILEX(25,'@',0,0,.0,.0,3,.FALSE.)
        END IF
C------------------------------------------------------------------------------
        DO I=1,K
          WRITE(*,117)'Spectrum #',I,'  -> Input file name'
          IF(I.EQ.1)THEN
            FILENAME=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
          ELSE
            FILENAME=INFILEX(20,'@',NSCAN,NCHAN2,STWV2,DISP2,1,.FALSE.)
            IF(NCHAN2.NE.NCHAN)THEN
              WRITE(*,101)'WARNING: NCHAN in last spectrum is '//
     +         'different.'
            END IF
            IF((STWV.NE.STWV2).OR.(DISP.NE.DISP2))THEN
              WRITE(*,101)'WARNING: STWV and/or DISP in last '//
     +         'spectrum are different.'
            END IF
          END IF
          IF(NSCAN.NE.1)THEN
            WRITE(*,101)'WARNING: last file is an image and not '//
     +       'a spectrum.'
            WRITE(*,100)'Which spectrum do you want to read '
            NSKIP=READILIM('1',1,NSCAN)
          ELSE
            NSKIP=1
          END IF
          DO II=1,NSKIP
            READ(20) (S(J),J=1,NCHAN)
          END DO
          DO J=1,NCHAN
            A(J,I)=S(J)
          END DO
          CLOSE(20)
          IF(CERR.EQ.'y')THEN
            WRITE(*,100)'Error file name '
            CALL GUESSEF(FILENAME,ERRFILE)
            IF(I.EQ.1)THEN
              ERRFILE=INFILEX(20,ERRFILE,NSCAN,NCHAN,STWV,DISP,21,
     +         .TRUE.)                                               !match
            ELSE
              ERRFILE=INFILEX(20,ERRFILE,NSCAN,NCHAN2,STWV2,DISP2,21,
     +         .TRUE.)                                               !match
            END IF
            DO II=1,NSKIP
              READ(20) (S(J),J=1,NCHAN)
            END DO
            DO J=1,NCHAN
              ERR(J,I)=S(J)
            END DO
            CLOSE(20)
          END IF
          IF(COUTHEADER.EQ.'y')THEN
            WRITE(25,'(A20,1X,A20,$)')OBJECT,FITSFILE
            WRITE(25,'(2(1X,I6),$)')NSCAN,NCHAN
            WRITE(25,'(1X,F8.2,1X,F8.3,$)')STWV,DISP
            IF(AIRMASS.LE.0.)THEN
              WRITE(25,'(1X,F8.6,1X,F8.1)')0.0,TIMEXPOS
            ELSE
              WRITE(25,'(1X,F8.6,1X,F8.1)')AIRMASS,TIMEXPOS
            END IF
          END IF
        END DO
C------------------------------------------------------------------------------
        IF(COUTHEADER.EQ.'y')THEN
          CLOSE(25)
        END IF
C------------------------------------------------------------------------------
C alteramos descriptores de cabecera
        AIRMASS=0.
        TIMEXPOS=0.
        OBJECT='[from gluesc]'
        FITSFILE='[from gluesc]'
C
        WRITE(*,100)'Outpuf file name'
        OUTFILE=OUTFILEX(30,'@',K,NCHAN,STWV,DISP,1,.FALSE.)
        DO I=1,K
          WRITE(30) (A(J,I),J=1,NCHAN)
        END DO
        CLOSE(30)
        IF(CERR.EQ.'y')THEN
          WRITE(*,100)'Output error file name '
          CALL GUESSEF(OUTFILE,ERRFILE)
          ERRFILE=OUTFILEX(30,ERRFILE,K,NCHAN,STWV,DISP,1,.TRUE.)
          DO I=1,K
            WRITE(30) (ERR(J,I),J=1,NCHAN)
          END DO
          CLOSE(30)
        END IF
C
        WRITE(*,110)'No. of saved scans: ',K
C------------------------------------------------------------------------------
        STOP
100     FORMAT(A,$)
101     FORMAT(A)
110     FORMAT(A,I6)
117     FORMAT(A,I5,3X,A,$)
        END
