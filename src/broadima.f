C------------------------------------------------------------------------------
C Version 30-January-1997                                      file: broadima.f
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
C Program: broadima
C Classification: arithmetic & manipulations
C Description: Broaden selected spectra of an image.
C
Comment
C
C Ensancha espectros de una imagen introduciendo una dispersion de velocidades.
C
        PROGRAM BROADIMA
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        REAL READF
C
        REAL C
        PARAMETER (C=299792.458)
C
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
        CHARACTER*75 INFILE,OUTFILE,ERRFILE,POLFILE
        LOGICAL IFSCAN(NSMAX)
C------------------------------------------------------------------------------
        THISPROGRAM='broadima'
        CALL WELCOME('30-January-1997')
        CALL SHOWHLP('explanation')
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
        DO I=1,NSCAN
          IFSCAN(I)=.FALSE.
        END DO
10      WRITE(*,101)'* Define scans to be employed: '
11      WRITE(*,100)'Scans (0,0=EXIT) '
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
C------------------------------------------------------------------------------
        IF(CERR.EQ.'y')THEN
          CALL SHOWHLP('error types')
          WRITE(*,101)'1- compute error spectrum like a normal spectrum'
          WRITE(*,101)'2- compute error of the broadened spectrum'
          WRITE(*,100)'Option (1/2) '
          CTYPEERR(1:1)=READC('1','12')
        END IF
C------------------------------------------------------------------------------
30      WRITE(*,100)'Sigma (km/sec) (-1=from external file)'
        SIGMA=READF('@')
        IF(SIGMA.LT.0)THEN
          WRITE(*,100) 'File name with spectrum with FWHM(Angs)'
          WRITE(*,100) ' for each channel'
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
C------------------------------------------------------------------------------
        WRITE(*,100)'Output file name'
        OUTFILE=OUTFILEX(30,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        IF(CERR.EQ.'y')THEN
          WRITE(*,100)'Output error file name '
          CALL GUESSEF(OUTFILE,ERRFILE)
          ERRFILE=OUTFILEX(32,ERRFILE,NSCAN,NCHAN,STWV,DISP,1,.TRUE.)
        END IF
C------------------------------------------------------------------------------
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
