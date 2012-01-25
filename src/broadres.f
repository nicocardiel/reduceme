C------------------------------------------------------------------------------
C Version 20-July-1999                                         file: broadres.f
C------------------------------------------------------------------------------
C Copyright N. Cardiel & J. Gorgas, Departamento de Astrofisica
C Universidad Complutense de Madrid, 28040-Madrid, Spain
C E-mail: ncl@astrax.fis.ucm.es or fjg@astrax.fis.ucm.es
C------------------------------------------------------------------------------
C This program is free software; you can redistribute it and/or modify it
C under the terms of the GNU General Public License as published by the Free
C Software Foundation; either version 2 of the License, or (at your option) any
C later version. See the file gnu-public-license.txt for details.
C------------------------------------------------------------------------------
Comment
C
C Program: broadres
C Classification: arithmetic & manipulations
C Description: Broadens selected spectra of an image to produce spectra at a 
C fixed spectral resolution. It needs a table with velocity dispersions for
C each spectrum (scan).
C
Comment
C
C Ensancha espectros de una imagen introduciendo una dispersion de velocidades.
C
        PROGRAM BROADRES
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        REAL READF
C
        INTEGER I,J,L
        INTEGER N1,N2
        INTEGER NADD
        REAL S(NCMAX),SS(NCMAX)
        REAL ES(NCMAX),ESS(NCMAX)
        REAL SIGMA,SIGFINAL,SIGMAVEC(NCMAX)
        REAL SIGI(NSMAX)
        INTEGER IDUM
        CHARACTER*10 CHDEF,CHDEF2
        CHARACTER*1 CERR,CTYPEERR
        CHARACTER*75 INFILE,OUTFILE,ERRFILE,TABLEFILE,CHDUM
        LOGICAL IFSCAN(NSMAX)
C------------------------------------------------------------------------------
        THISPROGRAM='broadres'
        CALL WELCOME('20-July-1999')
C
        WRITE(*,100)'Work with error images (y/n) '
        CERR(1:1)=READC('n','yn')
C------------------------------------------------------------------------------
C Fichero de trabajo
        WRITE(*,100)'Input file name'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        IF(CERR.EQ.'y')THEN
          WRITE(*,100)'Input error file name '
          CALL GUESSEF(INFILE,ERRFILE)
          ERRFILE=INFILEX(22,ERRFILE,NSCAN,NCHAN,STWV,DISP,21,.TRUE.)!....match
        END IF
C------------------------------------------------------------------------------
C Si es una imagen, seleccionamos espectros a ensanchar
        IF(NSCAN.EQ.1)THEN
          IFSCAN(1)=.TRUE.
          GOTO 30
        END IF
C
        WRITE(CHDEF,240) NSCAN
 240    FORMAT('1,',I3)
        CALL RMBLANK(CHDEF,CHDEF2,L)
        DO I=1,NSCAN
          IFSCAN(I)=.FALSE.
        END DO
10      WRITE(*,101)'* Define scans to be employed: '
11      WRITE(*,100)'Scans (0,0=EXIT) '
        CALL READ2I(CHDEF2(1:L),N1,N2)
        IF((N1.EQ.0).AND.(N2.EQ.0)) GOTO 20
        IF((N1.LT.1).OR.(N2.GT.NSCAN).OR.(N1.GT.N2))THEN
          WRITE(*,101)'ERROR: numbers out of range. Try again.'
          GOTO 11
        END IF
        DO I=N1,N2
          IFSCAN(I)=.TRUE.
        END DO
        GOTO 11
20      NADD=0
        DO I=1,NSCAN
          IF(IFSCAN(I)) NADD=NADD+1
        END DO
        IF(NADD.EQ.0)THEN
          WRITE(*,101)'ERROR: number of scans to be employed = 0!'
          GOTO 10
        END IF
C------------------------------------------------------------------------------
        WRITE(*,100)'Input file with individual sigmas'
        TABLEFILE(1:75)=READC('@','@')
        OPEN(UNIT=24,FILE=TABLEFILE,STATUS='OLD',FORM='FORMATTED')
C This file must go from N1 to N2, with free format (two columns: I, sigma)
        DO I=1,NSCAN
           IF(IFSCAN(I)) THEN
             READ(24,'(A)') CHDUM
             READ(CHDUM,*) IDUM,SIGI(I)
             IF(IDUM.NE.I) THEN
               WRITE(*,101) 'ERROR IN SCAN NUMBERS (CHECK INPUT TABLE)' 
               STOP
             END IF
           END IF
        END DO
        CLOSE(24)
C------------------------------------------------------------------------------
        IF(CERR.EQ.'y')THEN
          WRITE(*,101)'1 - compute error spectrum like a normal'//
     +     ' spectrum'
          WRITE(*,101)'2 - compute error of the broadened spectrum'
          WRITE(*,100)'Option (1/2) '
          CTYPEERR(1:1)=READC('1','12')
        END IF
C------------------------------------------------------------------------------
30      WRITE(*,100)'Final sigma (km/sec)'
        SIGFINAL=READF('@')
C------------------------------------------------------------------------------
        WRITE(*,100)'Output file name'
        OUTFILE=OUTFILEX(30,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        IF(CERR.EQ.'y')THEN
          WRITE(*,100)'Output error file name '
          CALL GUESSEF(OUTFILE,ERRFILE)
          ERRFILE=OUTFILEX(32,ERRFILE,NSCAN,NCHAN,STWV,DISP,1,.TRUE.)
        END IF
C------------------------------------------------------------------------------
        DO I=1,NSCAN
          READ(20) (S(J),J=1,NCHAN)
          IF(CERR.EQ.'y')THEN
            READ(22) (ES(J),J=1,NCHAN)
          END IF
          IF(IFSCAN(I))THEN
            IF(SIGFINAL.GT.SIGI(I)) THEN
              SIGMA=SQRT(SIGFINAL*SIGFINAL-SIGI(I)*SIGI(I))
              WRITE(*,210) I,SIGMA
 210          FORMAT('Broadening spectrum no. ',I3,' with sigma =',
     +         F5.1)
              DO J=1,NCHAN
                SIGMAVEC(J)=SIGMA
              END DO
              CALL BROADEN(S,SS,NCHAN,STWV,DISP,SIGMAVEC,.FALSE.)
              WRITE(30) (SS(J),J=1,NCHAN)
              IF(CERR.EQ.'y')THEN
                IF(CTYPEERR.EQ.'1')THEN
                  CALL BROADEN(ES,ESS,NCHAN,STWV,DISP,SIGMAVEC,.FALSE.)
                ELSE
                  CALL BROADEN(ES,ESS,NCHAN,STWV,DISP,SIGMAVEC,.TRUE.)
                END IF
                WRITE(32) (ESS(J),J=1,NCHAN)
              END IF
            ELSE
              WRITE(*,220) I
 220          FORMAT('Cannot broaden spectrum no. ',I3)
              WRITE(30) (S(J),J=1,NCHAN)
              IF(CERR.EQ.'y')THEN
                WRITE(32) (ES(J),J=1,NCHAN)
              END IF
            END IF
          ELSE
            WRITE(30) (S(J),J=1,NCHAN)
            IF(CERR.EQ.'y')THEN
              WRITE(32) (ES(J),J=1,NCHAN)
            END IF
          END IF
        END DO
C------------------------------------------------------------------------------
        CLOSE(20)
        CLOSE(30)
        IF(CERR.EQ.'y')THEN
          CLOSE(22)
          CLOSE(32)
        END IF
C
        STOP
100     FORMAT(A,$)
101     FORMAT(A)
        END
        
