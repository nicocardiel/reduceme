!------------------------------------------------------------------------------
! Version 20-November-2026                                  file:testrandom.f90
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
! Program: testrandom
! Classification: miscellany
! Description: Program to test the random number generator.
!
!Comment
!
        PROGRAM TESTRANDOM
        USE RED_M
        USE FSPEC_M
        USE FUTILS_M
        USE REDLIB_M
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
!
!
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        INTEGER I,I0,NSEED,NSIMUL
        REAL X(1001),Y(1001)
        REAL X0
        REAL MEAN0
        REAL MEAN,SIGMA
        REAL XMIN,XMAX,YMIN,YMAX
        CHARACTER*50 CDUMMY
        LOGICAL LCOLOR(MAX_ID_RED)
        LOGICAL LUNO
!------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
!
        THISPROGRAM='testrandom'
        CALL WELCOME('20-November-2026')
!
        WRITE(*,100,ADVANCE='NO')'Test #1...'
        NSEED=-1
        LUNO=.FALSE.
        X0=RANRED(NSEED)
        IF(X0.GT.1)THEN
          WRITE(*,101)'FATAL ERROR: random number .gt. 1.0'
          STOP
        ELSEIF(X0.EQ.1.0)THEN
          LUNO=.TRUE.
        END IF
        DO I=1,1000
          X0=RANRED(NSEED)
          IF(X0.GT.1)THEN
            WRITE(*,101)'FATAL ERROR: random number .gt. 1.0'
            STOP
          END IF
        END DO
        WRITE(*,101)'   ...OK!'
        IF(LUNO) WRITE(*,101)'WARNING: 1.0 value reached.'
!------------------------------------------------------------------------------
        WRITE(*,101)'Test #2...'
        WRITE(*,100,ADVANCE='NO')'No. of simulations '
        NSIMUL=READILIM('10000',100,100000)
        MEAN0=REAL(NSIMUL)/1000.
        CALL PIDEGTER(NTERM,IDN,LCOLOR)
        WRITE(*,100,ADVANCE='NO')'Xmin '
        XMIN=READF('-0.1')
        WRITE(*,100,ADVANCE='NO')'Xmax '
        XMAX=READF('+1.1')
        WRITE(*,100,ADVANCE='NO')'Ymin '
        YMIN=READF('0.0')
        WRITE(*,100,ADVANCE='NO')'Ymax '
        WRITE(CDUMMY,*) 2.*MEAN0
        YMAX=READF(CDUMMY)
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,0)
          IF(LCOLOR(ITERM)) CALL PGSCI(2)
        END DO
!
        DO I=1,1000
          X(I)=REAL(I-1)/1000.
          Y(I)=0.
        END DO
!
        DO I=1,NSIMUL
          X0=RANRED(NSEED)
          I0=INT(X0*1000)+1
          Y(I0)=Y(I0)+1
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGMOVE(X(I0),0.)
            CALL PGDRAW(X(I0),Y(I0))
          END DO
        END DO
!
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM)) CALL PGSCI(3)
          CALL PGMOVE(XMIN,MEAN0)
          CALL PGDRAW(XMAX,MEAN0)
        END DO
!
        MEAN=0.
        DO I=1,1000
          MEAN=MEAN+Y(I)
        END DO
        MEAN=MEAN/1000.
        SIGMA=0.
        DO I=1,1000
          SIGMA=SIGMA+(Y(I)-MEAN)*(Y(I)-MEAN)
        END DO
        SIGMA=SQRT(SIGMA/999.)
        WRITE(*,100,ADVANCE='NO')'>> Mean value........: '
        WRITE(*,*) MEAN
        WRITE(*,100,ADVANCE='NO')'>> Standard deviation: '
        WRITE(*,*) SIGMA
!
        CALL PGEND
!------------------------------------------------------------------------------
        STOP
100     FORMAT(A)
101     FORMAT(A)
        END
