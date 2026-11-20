!------------------------------------------------------------------------------
! Version 20-November-2026                                   file: broadres.f90
!------------------------------------------------------------------------------
! Copyright N. Cardiel & J. Gorgas, Departamento de Astrofisica
! Universidad Complutense de Madrid, 28040-Madrid, Spain
! E-mail: ncl@astrax.fis.ucm.es or fjg@astrax.fis.ucm.es
!------------------------------------------------------------------------------
! This program is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the Free
! Software Foundation; either version 2 of the License, or (at your option) any
! later version. See the file gnu-public-license.txt for details.
!------------------------------------------------------------------------------
!Comment
!
! Program: broadres
! Classification: arithmetic & manipulations
! Description: Broadens selected spectra of an image to produce spectra at a 
! fixed spectral resolution. It needs a table with velocity dispersions for
! each spectrum (scan).
!
!Comment
!
! Ensancha espectros de una imagen introduciendo una dispersion de velocidades.
!
        PROGRAM BROADRES
        USE RED_M
        USE FSPEC_M
        USE FUTILS_M
        USE REDLIB_M
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
!
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
        CHARACTER*75 CHDUM
        CHARACTER*255 TABLEFILE
        CHARACTER*255 INFILE,OUTFILE,ERRFILE
        LOGICAL IFSCAN(NSMAX)
!------------------------------------------------------------------------------
        THISPROGRAM='broadres'
        CALL WELCOME('20-November-2026')
!
        WRITE(*,100,ADVANCE='NO')'Work with error images (y/n) '
        CERR(1:1)=READC('n','yn')
!------------------------------------------------------------------------------
! Fichero de trabajo
        WRITE(*,100,ADVANCE='NO')'Input file name'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        IF(CERR.EQ.'y')THEN
          WRITE(*,100,ADVANCE='NO')'Input error file name '
          CALL GUESSEF(INFILE,ERRFILE)
          ERRFILE=INFILEX(22,ERRFILE,NSCAN,NCHAN,STWV,DISP,21,.TRUE.)!....match
        END IF
!------------------------------------------------------------------------------
! Si es una imagen, seleccionamos espectros a ensanchar
        IF(NSCAN.EQ.1)THEN
          IFSCAN(1)=.TRUE.
          GOTO 30
        END IF
!
        WRITE(CHDEF,240) NSCAN
 240    FORMAT('1,',I3)
        CALL RMBLANK(CHDEF,CHDEF2,L)
        DO I=1,NSCAN
          IFSCAN(I)=.FALSE.
        END DO
10      WRITE(*,101)'* Define scans to be employed: '
11      WRITE(*,100,ADVANCE='NO')'Scans (0,0=EXIT) '
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
!------------------------------------------------------------------------------
        WRITE(*,100,ADVANCE='NO')'Input file with individual sigmas'
        TABLEFILE=READC_FILE('@','@')
        OPEN(UNIT=24,FILE=TABLEFILE,STATUS='OLD',FORM='FORMATTED')
! This file must go from N1 to N2, with free format (two columns: I, sigma)
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
!------------------------------------------------------------------------------
        IF(CERR.EQ.'y')THEN
          WRITE(*,101)'1 - compute error spectrum like a normal'//&
     &     ' spectrum'
          WRITE(*,101)'2 - compute error of the broadened spectrum'
          WRITE(*,100,ADVANCE='NO')'Option (1/2) '
          CTYPEERR(1:1)=READC('1','12')
        END IF
!------------------------------------------------------------------------------
30      WRITE(*,100,ADVANCE='NO')'Final sigma (km/sec)'
        SIGFINAL=READF('@')
!------------------------------------------------------------------------------
        WRITE(*,100,ADVANCE='NO')'Output file name'
        OUTFILE=OUTFILEX(30,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        IF(CERR.EQ.'y')THEN
          WRITE(*,100,ADVANCE='NO')'Output error file name '
          CALL GUESSEF(OUTFILE,ERRFILE)
          ERRFILE=OUTFILEX(32,ERRFILE,NSCAN,NCHAN,STWV,DISP,1,.TRUE.)
        END IF
!------------------------------------------------------------------------------
        DO I=1,NSCAN
          READ(20) (S(J),J=1,NCHAN)
          IF(CERR.EQ.'y')THEN
            READ(22) (ES(J),J=1,NCHAN)
          END IF
          IF(IFSCAN(I))THEN
            IF(SIGFINAL.GT.SIGI(I)) THEN
              SIGMA=SQRT(SIGFINAL*SIGFINAL-SIGI(I)*SIGI(I))
              WRITE(*,210) I,SIGMA
 210          FORMAT('Broadening spectrum no. ',I3,' with sigma =',     &
     &         F5.1)
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
!------------------------------------------------------------------------------
        CLOSE(20)
        CLOSE(30)
        IF(CERR.EQ.'y')THEN
          CLOSE(22)
          CLOSE(32)
        END IF
!
        STOP
100     FORMAT(A)
101     FORMAT(A)
        END
        
