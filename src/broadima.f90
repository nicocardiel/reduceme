!------------------------------------------------------------------------------
! Version 20-November-2026                                   file: broadima.f90
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
! Program: broadima
! Classification: arithmetic & manipulations
! Description: Broaden selected spectra of an image.
!
!Comment
!
! Ensancha espectros de una imagen introduciendo una dispersion de velocidades.
!
        PROGRAM BROADIMA
        USE RED_M
        USE FSPEC_M
        USE FUTILS_M
        USE REDLIB_M
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
!
        REAL C
        PARAMETER (C=299792.458)
!
        INTEGER I,J
        INTEGER N1,N2
        INTEGER NADD,NSCURRENT
        INTEGER NEXTINFO
        INTEGER NSCAN_,NCHAN_
        REAL STWV_,DISP_,FLDO
        REAL S(NCMAX),SS(NCMAX)
        REAL ES(NCMAX),ESS(NCMAX)
        REAL SIGMA,SIGMAVEC(NCMAX)
        CHARACTER*1 CERR,CTYPEERR
        CHARACTER*255 INFILE,OUTFILE,ERRFILE,POLFILE
        LOGICAL IFSCAN(NSMAX)
!------------------------------------------------------------------------------
        THISPROGRAM='broadima'
        CALL WELCOME('20-November-2026')
        CALL SHOWHLP('explanation')
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
        DO I=1,NSCAN
          IFSCAN(I)=.FALSE.
        END DO
10      WRITE(*,101)'* Define scans to be employed: '
11      WRITE(*,100,ADVANCE='NO')'Scans (0,0=EXIT) '
        CALL READ2I('0,0',N1,N2)
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
        IF(CERR.EQ.'y')THEN
          CALL SHOWHLP('error types')
          WRITE(*,101)'1- compute error spectrum like a normal spectrum'
          WRITE(*,101)'2- compute error of the broadened spectrum'
          WRITE(*,100,ADVANCE='NO')'Option (1/2) '
          CTYPEERR(1:1)=READC('1','12')
        END IF
!------------------------------------------------------------------------------
30      WRITE(*,100,ADVANCE='NO')&
     & 'Sigma (km/sec) (-1=from external file)'
        SIGMA=READF('@')
        IF(SIGMA.LT.0)THEN
          WRITE(*,100,ADVANCE='NO')&
     &       'File name with spectrum with FWHM(Angs)'
          WRITE(*,100,ADVANCE='NO') ' for each channel'
          POLFILE=INFILEX(30,'@',NSCAN_,NCHAN_,STWV_,DISP_,1,.FALSE.)
          IF(NSCAN_.NE.1)THEN
            WRITE(*,101) 'ERROR: NSCAN.NE.1'
            CLOSE(30)
            STOP
          ELSEIF(NCHAN_.NE.NCHAN)THEN
            WRITE(*,101) 'ERROR: NCHAN_.NE.NCHAN'
            CLOSE(30)
            STOP
          END IF
          READ(30) (SIGMAVEC(J),J=1,NCHAN)
          CLOSE(30)
          DO J=1,NCHAN
            SIGMAVEC(J)=SIGMAVEC(J)/2.35 !FWHM (Angs.) -> sigma (Angs.)
            FLDO=REAL(J-1)*DISP+STWV
            SIGMAVEC(J)=SIGMAVEC(J)/FLDO*C !sigma (Angs.) -> sigma (km/s)
          END DO
        ELSE
          DO J=1,NCHAN
            SIGMAVEC(J)=SIGMA
          END DO
        END IF
!------------------------------------------------------------------------------
        WRITE(*,100,ADVANCE='NO')'Output file name'
        OUTFILE=OUTFILEX(30,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        IF(CERR.EQ.'y')THEN
          WRITE(*,100,ADVANCE='NO')'Output error file name '
          CALL GUESSEF(OUTFILE,ERRFILE)
          ERRFILE=OUTFILEX(32,ERRFILE,NSCAN,NCHAN,STWV,DISP,1,.TRUE.)
        END IF
!------------------------------------------------------------------------------
        NSCURRENT=0
        DO I=1,NSCAN
          READ(20) (S(J),J=1,NCHAN)
          IF(CERR.EQ.'y')THEN
            READ(22) (ES(J),J=1,NCHAN)
          END IF
          IF(IFSCAN(I))THEN
            NSCURRENT=NSCURRENT+1
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
            WRITE(30) (S(J),J=1,NCHAN)
            IF(CERR.EQ.'y')THEN
              WRITE(32) (ES(J),J=1,NCHAN)
            END IF
          END IF
          CALL SHOWPERC(1,NADD,1,NSCURRENT,NEXTINFO)
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
