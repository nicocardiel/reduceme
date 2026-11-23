!------------------------------------------------------------------------------
! Version 24-November-2026                                    file: findmax.f90
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
! Program: findmax
! Classification: wavelengths
! Description: Automatic detection of line peaks in a spectrum.
!
!Comment
!
! busca lineas en un espectro
!
        PROGRAM FINDMAX
        USE RED_M
        USE FUTILS_M
        USE REDLIB_M
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
!
        INTEGER I,J,K
        INTEGER N,NMED
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        REAL S(NCMAX)
        REAL X(NCMAX)
        REAL XMIN,XMAX,YMIN,YMAX,DX,DY
        REAL XX,YY
        REAL DATAMIN
        CHARACTER*1 COPC,CSAVE,CREPEAT
        CHARACTER*50 CDUMMY
        CHARACTER*255 INFILE,OUTFILE
        LOGICAL LCOLOR(MAX_ID_RED),LGOOD
!------------------------------------------------------------------------------
        THISPROGRAM='findmax'
        CALL WELCOME('24-November-2026')
!
        CALL PIDEGTER(NTERM,IDN,LCOLOR)
!
        WRITE(*,100,ADVANCE='NO')'Spectrum file name'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        IF(TRUELEN(OBJECT).GT.0)THEN
          INFILE=INFILE(1:TRUELEN(INFILE))//' ['//&
     &     OBJECT(1:TRUELEN(OBJECT))//']'
        END IF
        READ(20) (S(I),I=1,NCHAN)
        CLOSE(20)
!
        DO I=1,NCHAN
          X(I)=REAL(I)
        END DO
!
        YMIN=S(1)
        YMAX=YMIN
        DO I=2,NCHAN
          IF(S(I).LT.YMIN) YMIN=S(I)
          IF(S(I).GT.YMAX) YMAX=S(I)
        END DO
        XMIN=1.
        XMAX=REAL(NCHAN)
        DX=XMAX-XMIN
        XMIN=XMIN-DX/50.
        XMAX=XMAX+DX/50.
        DY=YMAX-YMIN
        YMIN=YMIN-DY/50.
        YMAX=YMAX+DY/50.
!
10      DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,0)
          CALL PGIDEN_RED
          CALL PGBIN(NCHAN,X,S,.TRUE.)
          CALL PGLABEL('channel','No. counts',INFILE)
        END DO
!
20      WRITE(*,100,ADVANCE='NO')'Change Y-limits (y/n) '
        COPC(1:1)=READC('n','yn')
        IF(COPC.EQ.'y')THEN
          WRITE(*,100,ADVANCE='NO')'YMIN '
          WRITE(CDUMMY,*) YMIN
          YMIN=READF(CDUMMY)
          WRITE(*,100,ADVANCE='NO')'YMAX '
          WRITE(CDUMMY,*) YMAX
          YMAX=READF(CDUMMY)
          DY=YMAX-YMIN
          YMIN=YMIN-DY/50.
          YMAX=YMAX+DY/50.
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,0)
            CALL PGIDEN_RED
            CALL PGBIN(NCHAN,X,S,.TRUE.)
            CALL PGLABEL('channel','No. counts',INFILE)
          END DO
          GOTO 20
        END IF
!
        WRITE(*,101)'Maximum will be searched in a region of N '//&
     &   'channels.'
        WRITE(*,100,ADVANCE='NO')'N (must be odd) '
        N=READI('5')
        IF(MOD(N,2).EQ.0)THEN
          N=N+1
          WRITE(*,'(A,I3)')'Effective N to be employed: ',N
        END IF
        NMED=N/2
        WRITE(*,'(A,I3,A)')'First and last ',NMED,' channels will '//&
     &   'be ignored.'
!
        WRITE(*,100,ADVANCE='NO')'Minimum threshold'
        DATAMIN=READF('@')
!
        K=0
        DO I=NMED+1,NCHAN-NMED-1
          LGOOD=.TRUE.
          J=0
          DO WHILE((J.LT.NMED).AND.(LGOOD))
            J=J+1
            IF(S(I-NMED+J-1).GT.S(I-NMED+J)) LGOOD=.FALSE.
          END DO
          J=NMED
          DO WHILE((J.LT.N-1).AND.(LGOOD))
            J=J+1
            IF(S(I-NMED+J-1).LT.S(I-NMED+J)) LGOOD=.FALSE.
          END DO
          IF((S(I).GE.DATAMIN).AND.(LGOOD))THEN
            K=K+1
            XX=REAL(I)
            YY=S(I)
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              IF(LCOLOR(ITERM)) CALL PGSCI(2)
              CALL PGPOINT(1,XX,YY,4)
              IF(LCOLOR(ITERM)) CALL PGSCI(1)
            END DO
          END IF
        END DO
        WRITE(*,110)'No. of lines found: ',K
!
        WRITE(*,100,ADVANCE='NO')'Repeat search (y/n) '
        CREPEAT(1:1)=READC('n','yn')
        IF(CREPEAT.EQ.'y') GOTO 10
!
        WRITE(*,100,ADVANCE='NO')'Save output into a file (y/n) '
        CSAVE(1:1)=READC('n','yn')
        IF(CSAVE.EQ.'y')THEN
          WRITE(*,100,ADVANCE='NO')'Output file name'
          OUTFILE=OUTFILEX(30,'@',0,0,.0,.0,3,.FALSE.)
          K=0
          DO I=NMED+1,NCHAN-NMED-1
            LGOOD=.TRUE.
            J=0
            DO WHILE((J.LT.NMED).AND.(LGOOD))
              J=J+1
              IF(S(I-NMED+J-1).GT.S(I-NMED+J)) LGOOD=.FALSE.
            END DO
            J=NMED
            DO WHILE((J.LT.N-1).AND.(LGOOD))
              J=J+1
              IF(S(I-NMED+J-1).LT.S(I-NMED+J)) LGOOD=.FALSE.
            END DO
            IF((S(I).GE.DATAMIN).AND.(LGOOD)) WRITE(30,*)I
          END DO
          CALL CLOSEFILEX(30)
        END IF
!
        CALL PGEND
        STOP
100     FORMAT(A)
101     FORMAT(A)
110     FORMAT(A,I6)
        END
