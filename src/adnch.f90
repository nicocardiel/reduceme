!------------------------------------------------------------------------------
! Version 20-November-2026                                       file:adnch.f90
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
! Program: adnch
! Classification: arithmetic & manipulations
! Description: Add channels of an image, generating a new image with the 
! same NSCAN and variable NCHAN.
!
!Comment
!
        PROGRAM ADNCH
        USE RED_M
        USE FUTILS_M
        USE REDLIB_M
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
!
        INTEGER I,J
        INTEGER K
        INTEGER N,N1,N2
        REAL S(NSMAX)
        REAL, ALLOCATABLE :: A(:,:),ERR(:,:)
        REAL ES(NSMAX)
        REAL, ALLOCATABLE :: B(:,:),ERRB(:,:)
        REAL, ALLOCATABLE :: AUX(:,:)
        REAL FACTOR
        CHARACTER*1 CERR
        CHARACTER*255 FILENAME,ERRFILE,OUTFILE
!------------------------------------------------------------------------------
        THISPROGRAM='adnch'
        CALL WELCOME('20-November-2026')
        CALL SHOWHLP('explanation')
!
        WRITE(*,100,ADVANCE='NO')'Work with error images (y/n) '
        CERR(1:1)=READC('n','yn')
!
        WRITE(*,100,ADVANCE='NO')'Input file name'
        FILENAME=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        ! allocate with the actual image dimensions (static arrays
        ! were limited by NCMAX,NSMAX); initialised as before
        IF(ALLOCATED(A)) DEALLOCATE(A)
        IF(ALLOCATED(ERR)) DEALLOCATE(ERR)
        IF(ALLOCATED(B)) DEALLOCATE(B)
        IF(ALLOCATED(ERRB)) DEALLOCATE(ERRB)
        ALLOCATE(A(NCHAN,NSCAN),ERR(NCHAN,NSCAN),B(NCHAN,NSCAN),ERRB(NCHAN,NSCAN))
        A=0.0
        ERR=0.0
        B=0.0
        ERRB=0.0
        DO I=1,NSCAN
          READ(20) (A(J,I),J=1,NCHAN)
        END DO
        CLOSE(20)
        IF(CERR.EQ.'y')THEN
          WRITE(*,100,ADVANCE='NO')'Input error file name '
          CALL GUESSEF(FILENAME,ERRFILE)
          ERRFILE=INFILEX(21,ERRFILE,NSCAN,NCHAN,STWV,DISP,21,.TRUE.)!....match
          DO I=1,NSCAN
            READ(21) (ERR(J,I),J=1,NCHAN)
          END DO
          CLOSE(21)
        END IF
!------------------------------------------------------------------------------
        K=0
!
50      CALL SHOWHLP('selection')
        WRITE(*,100,ADVANCE='NO')'1st Channel (0=Exit) '
        N1=READILIM('@',0,NCHAN)
        IF(N1.EQ.0) GO TO 60
        WRITE(*,100,ADVANCE='NO')'No. of Channels '
        N=READILIM('@',1,NCHAN-N1+1)
        WRITE(*,100,ADVANCE='NO')'Factor'
        FACTOR=READF('@')
        N2=N1+N-1
        DO I=1,NSCAN
          S(I)=0.
        END DO
        DO I=1,NSCAN
          DO J=N1,N2
            S(I)=S(I)+A(J,I)*FACTOR
          END DO
        END DO
        IF(CERR.EQ.'y')THEN
          DO I=1,NSCAN
            ES(I)=0.
          END DO
          DO I=1,NSCAN
            DO J=N1,N2
              ES(I)=ES(I)+ERR(J,I)*ERR(J,I)
            END DO
          END DO
        END IF
        K=K+1
        IF(K.GT.NCMAX)THEN
          WRITE(*,101)'ERROR: maximum number of channels excedeed'
          STOP
        END IF
! B and ERRB grow as channels are added (doubling their size)
        IF(K.GT.SIZE(B,1))THEN
          ALLOCATE(AUX(2*SIZE(B,1),NSCAN))
          AUX=0.0
          AUX(1:SIZE(B,1),:)=B
          CALL MOVE_ALLOC(AUX,B)
          ALLOCATE(AUX(2*SIZE(ERRB,1),NSCAN))
          AUX=0.0
          AUX(1:SIZE(ERRB,1),:)=ERRB
          CALL MOVE_ALLOC(AUX,ERRB)
        END IF
        DO I=1,NSCAN
          B(K,I)=S(I)
        END DO
        IF(CERR.EQ.'y')THEN
          DO I=1,NSCAN
            ERRB(K,I)=FACTOR*SQRT(ES(I))
          END DO
        END IF
        GO TO 50
  60    CONTINUE
!------------------------------------------------------------------------------
        WRITE(*,100,ADVANCE='NO')'Output file name'
        OUTFILE=OUTFILEX(30,'@',NSCAN,K,0.,0.,1,.FALSE.)
        DO I=1,NSCAN
          WRITE(30) (B(J,I),J=1,K)
        END DO
        CLOSE(30)
        IF(CERR.EQ.'y')THEN
          WRITE(*,100,ADVANCE='NO')'Output error file name '
          CALL GUESSEF(OUTFILE,ERRFILE)
          OUTFILE=OUTFILEX(31,ERRFILE,NSCAN,K,0.,0.,1,.TRUE.)
          DO I=1,NSCAN
            WRITE(31) (ERRB(J,I),J=1,K)
          END DO
          CLOSE(31)
        END IF
!
        WRITE(*,110)'No. of saved channels: ',K
        STOP
100     FORMAT(A)
101     FORMAT(A)
110     FORMAT(A,I6)
        END
