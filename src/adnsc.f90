!------------------------------------------------------------------------------
! Version 24-November-2026                                      file: adnsc.f90
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
! Program: adnsc
! Classification: arithmetic & manipulations
! Description: Add scans of an image, generating a new image with the same 
! NCHAN and variable NSCAN.
!
!Comment
!
        PROGRAM ADNSC
        USE RED_M
        USE FUTILS_M
        USE REDLIB_M
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
!
        INTEGER I,J,K,L
        INTEGER N,N1,N2
        REAL S(NCMAX)
        REAL, ALLOCATABLE :: A(:,:),B(:,:)
        REAL, ALLOCATABLE :: AUX(:,:)
        REAL ES(NCMAX)
        REAL, ALLOCATABLE :: ERR(:,:),ERRB(:,:)
        REAL FACTOR
        CHARACTER*1 CERR,CSINGLE,CNOR,CFHEAD2
        CHARACTER*15 CSCANBINNING
        CHARACTER*255 FILENAME,ERRFILE,HEADFILE,OUTFILE
!------------------------------------------------------------------------------
        THISPROGRAM='adnsc'
        CALL WELCOME('24-November-2026')
        CALL SHOWHLP('explanation')
!
        WRITE(*,100,ADVANCE='NO')'Work with error images (y/n) '
        CERR(1:1)=READC('y','yn')
!
        CALL SHOWHLP('single spectrum')
        WRITE(*,100,ADVANCE='NO')'Output into a single spectrum (y/n) '
        CSINGLE(1:1)=READC('y','yn')
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
        IF(CSINGLE.EQ.'n')THEN
          CALL SHOWHLP('file to index')
          WRITE(*,100,ADVANCE='NO')&
     &       'Create log file (with binning regions) '//&
     &     'for index.f (y/n) '
          CFHEAD2(1:1)=READC('n','yn')
          IF(CFHEAD2.EQ.'y')THEN
            WRITE(*,100,ADVANCE='NO')'Log file name'
            HEADFILE=OUTFILEX(27,'@',0,0,0.,0.,3,.FALSE.)
          END IF
          GOTO 40
        END IF
!------------------------------------------------------------------------------
! Salida en un unico espectro
        DO J=1,NCHAN
          S(J)=0.
        END DO
        IF(CERR.EQ.'y')THEN
          DO J=1,NCHAN
            ES(J)=0.
          END DO
        END IF
!
        WRITE(*,*)
        CALL SHOWHLP('factors versus normalization')
        WRITE(*,101)'* If you are using different factors for '//&
     &   'each spectrum, DO NOT normalize'
        WRITE(*,100,ADVANCE='NO')'Normalize final spectrum (y/n) '
        CNOR(1:1)=READC('y','yn')
        IF(CNOR.EQ.'y')THEN
          FACTOR=1.
        END IF
!
        N=0
10      CALL SHOWHLP('enter scan region')
        WRITE(*,100,ADVANCE='NO')'Scan region (0,0=EXIT) '
        CALL READ2I('0,0',N1,N2)
        IF((N1.EQ.0).AND.(N2.EQ.0)) GOTO 20
        IF((N1.LT.0).OR.(N2.GT.NSCAN).OR.(N1.GT.N2))THEN
          WRITE(*,101)'ERROR: invalid entry. Try again.'
          GOTO 10
        END IF
        IF(CNOR.EQ.'n')THEN
          WRITE(*,100,ADVANCE='NO')'Factor '
          FACTOR=READF('1.0')
        END IF
        DO I=N1,N2
          DO J=1,NCHAN
            S(J)=S(J)+A(J,I)*FACTOR
          END DO
          IF(CERR.EQ.'y')THEN
            DO J=1,NCHAN
              ES(J)=ES(J)+ERR(J,I)*ERR(J,I)*FACTOR*FACTOR
            END DO
          END IF
        END DO
        N=N+(N2-N1+1)
        WRITE(*,110)'>>> Number of scans added......: ',N2-N1+1
        WRITE(*,110)'>>> Total number of scans added: ',N
        GOTO 10
20      IF(N.EQ.0)THEN
          WRITE(*,101)'ERROR: number of scans added = 0. Try again.'
          GOTO 10
        END IF
        IF(CNOR.EQ.'y')THEN
          DO J=1,NCHAN
            S(J)=S(J)/REAL(N)
          END DO
          IF(CERR.EQ.'y')THEN
            DO J=1,NCHAN
              ES(J)=SQRT(ES(J))/REAL(N)
            END DO
          END IF
        ELSE
          IF(CERR.EQ.'y')THEN
            DO J=1,NCHAN
              ES(J)=SQRT(ES(J))
            END DO
          END IF
        END IF
        K=1
        GOTO 70
!------------------------------------------------------------------------------
40      K=0
50      CALL SHOWHLP('scan region with increment')
        WRITE(*,100,ADVANCE='NO')'1st Scan (0=Exit) '
        N1=READILIM('@',0,NSCAN)
        IF(N1.EQ.0) GO TO 60
        WRITE(*,100,ADVANCE='NO')'No. of Scans '
        N=READILIM('@',1,NSCAN-N1+1)
        WRITE(*,100,ADVANCE='NO')'Factor '
        FACTOR=READF('1.0')
        N2=N1+N-1
!
        IF(CFHEAD2.EQ.'y')THEN
          WRITE(CSCANBINNING,'(I7,A1,I7)')N1,',',N2
          CALL RMBLANK(CSCANBINNING,CSCANBINNING,L)
          WRITE(27,'(A15)') CSCANBINNING(1:L)
        END IF
!
        DO J=1,NCHAN
          S(J)=0.0
        END DO
        DO I=N1,N2
          DO J=1,NCHAN
            S(J)=S(J)+A(J,I)*FACTOR
          END DO
        END DO
        IF(CERR.EQ.'y')THEN
          DO J=1,NCHAN
            ES(J)=0.0
          END DO
          DO I=N1,N2
            DO J=1,NCHAN
              ES(J)=ES(J)+ERR(J,I)*ERR(J,I)
            END DO
          END DO
        END IF
!
        K=K+1
! B and ERRB grow as scans are added (doubling their size); the limit
! NSMAX below is kept
        IF(K.GT.SIZE(B,2))THEN
          ALLOCATE(AUX(NCHAN,2*SIZE(B,2)))
          AUX=0.0
          AUX(:,1:SIZE(B,2))=B
          CALL MOVE_ALLOC(AUX,B)
          ALLOCATE(AUX(NCHAN,2*SIZE(ERRB,2)))
          AUX=0.0
          AUX(:,1:SIZE(ERRB,2))=ERRB
          CALL MOVE_ALLOC(AUX,ERRB)
        END IF
        DO J=1,NCHAN
          B(J,K)=S(J)
        END DO
        IF(CERR.EQ.'y')THEN
          DO J=1,NCHAN
            ERRB(J,K)=FACTOR*SQRT(ES(J))
          END DO
        END IF
        IF(K.EQ.NSMAX)THEN
          WRITE(*,101)'WARNING: No. of scans = NSMAX'
          WRITE(*,101)'         Output file must be saved.'
          GOTO 60
        END IF
        GO TO 50
60      IF(CFHEAD2.EQ.'y') CALL CLOSEFILEX(27)
!------------------------------------------------------------------------------
! salvamos fichero
70      WRITE(*,100,ADVANCE='NO')'Output file name'
        OUTFILE=OUTFILEX(30,'@',K,NCHAN,STWV,DISP,1,.FALSE.)
        IF(CSINGLE.EQ.'y')THEN
          WRITE(30) (S(J),J=1,NCHAN)
        ELSE
          DO I=1,K
            WRITE(30) (B(J,I),J=1,NCHAN)
          END DO
        END IF
        CALL CLOSEFILEX(30)
!
        IF(CERR.EQ.'y')THEN
          WRITE(*,100,ADVANCE='NO')'Output error file name '
          CALL GUESSEF(OUTFILE,ERRFILE)
          OUTFILE=OUTFILEX(31,ERRFILE,K,NCHAN,STWV,DISP,1,.TRUE.)
          IF(CSINGLE.EQ.'y')THEN
            WRITE(31) (ES(J),J=1,NCHAN)
          ELSE
            DO I=1,K
              WRITE(31) (ERRB(J,I),J=1,NCHAN)
            END DO
          END IF
          CALL CLOSEFILEX(31)
        END IF
!
        STOP
100     FORMAT(A)
101     FORMAT(A)
110     FORMAT(A,I6)
        END
