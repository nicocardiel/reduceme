!------------------------------------------------------------------------------
! Version 20-November-2026                                        file: cdisc.f
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
! Program: cdisc
! Classification: distortion
! Description: Correct C-distortion using the polynomial fits of fitcdis.
!
!Comment
!
! Corrige de distorsion C empleando los ajustes polinomicos realizados
! por el programa fitcdis.f
!
        PROGRAM CDISC
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INTEGER READILIM
!
        INTEGER I,J,K,MJ
        INTEGER I0,K0
        INTEGER IBEG,IEND
        INTEGER NSPECTRA
        INTEGER NORDER,NTERMS
        REAL S(NCMAX),S2(NCMAX),COE(20,NSMAX),X(NCMAX),D(NCMAX)
        REAL XJ,XI
        REAL X1,X2,XT
        REAL FI,FI1,FI2
        REAL XMIN,XMAX
        REAL POL,FAC
        REAL SUM,SUM2,DIF,SUMX
        CHARACTER*80 CDUMMY
        CHARACTER*255 POLFILE,INFILE,OUTFILE
!------------------------------------------------------------------------------
        THISPROGRAM='cdisc'
        CALL WELCOME('20-November-2026')
        CALL SHOWHLP('explanation')
!
        WRITE(*,100,ADVANCE='NO')'Arc file name'
        INFILE=INFILEX(14,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
!
        WRITE(*,100,ADVANCE='NO')'Spectra from... to... '
        WRITE(CDUMMY,'(I4,A,I4)')1,',',NSCAN
        CALL READ2I(CDUMMY,IBEG,IEND)
        NSPECTRA=IEND-IBEG+1
!
        WRITE(*,101)'(From fitcdis) FIT of line peak deviations vs. '//&
     &   'channel number:'
        WRITE(*,100,ADVANCE='NO')'Degree of fitted polynomial '
        NORDER=READILIM('@',0,19)
        NTERMS=NORDER+1
!
        WRITE(*,100,ADVANCE='NO')'Polynomial file name (from fitcdis)'
        POLFILE=INFILEX(12,'@',0,0,0.,0.,3,.FALSE.)
!
        DO I=IBEG,IEND
          READ(12,112) XJ,(COE(K,I),K=1,NTERMS)
          IF(XJ.NE.FLOAT(I)) STOP 'ERROR IN SCAN NO.'
        END DO
112     FORMAT(1X,F5.0,20E20.9)
        CLOSE(12)
!
        DO I=1,IBEG-1
          DO K=1,NTERMS
            COE(K,I)=COE(K,IBEG)
          END DO
        END DO
        DO I=IEND+1,NSCAN
          DO K=1,NTERMS
            COE(K,I)=COE(K,IEND)
          END DO
        END DO
!
        WRITE(*,100,ADVANCE='NO')'Output file name'
        OUTFILE=OUTFILEX(15,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
!
        DO 80 MJ=1,NSCAN
!
          READ(14) (S(J),J=1,NCHAN)
!
          DO I=1,NCHAN
            XI=FLOAT(I)
            POL=COE(NTERMS,MJ)
            DO J=NTERMS-1,1,-1
              POL=POL*XI+COE(J,MJ)
            END DO
            X(I)=XI-POL
          END DO
          DO I=2,NCHAN
            D(I)=(X(I)-X(I-1))/2.
            IF(D(I).LE.0.) STOP 'ERROR IN DIFFERENCES'
          END DO
          D(1)=D(2)
!
          DO I=1,NCHAN
            IF((X(I)+D(I+1)).GT.0.5) GO TO 10
          END DO
10        I0=I
!
          DO K=1,NCHAN
            IF((X(I0)-D(I0)).LT.(FLOAT(K)+0.5)) GO TO 20
            S2(K)=0.
          END DO
20        K0=K
          DO I=K0,NCHAN
            S2(I)=0.
          END DO
!
          DO 30 I=I0,NCHAN
!
            X1=X(I)-D(I)
            X2=X(I)+D(I+1)
            XT=X2-X1
!
40          FI=FLOAT(K0)
            FI1=FI-0.5
            FI2=FI+0.5
            XMIN=AMAX1(X1,FI1)
            XMAX=AMIN1(X2,FI2)
            FAC=(XMAX-XMIN)/XT
!
            S2(K0)=S2(K0)+FAC*S(I)
!
            IF(FI2.LE.X2) THEN
              K0=K0+1
              IF(FI2.EQ.X2) GO TO 30
              IF(K0.GT.NCHAN) GO TO 50
              GO TO 40
            END IF
30      CONTINUE
50      CONTINUE
!
          SUM=0.
          SUM2=0.
          DO I=1,NCHAN
            SUM=SUM+S(I)
            SUM2=SUM2+S2(I)
          END DO
          IF(SUM.LE.0.) GO TO 60
          DIF=ABS(SUM-SUM2)/SUM
          SUMX=SUM/FLOAT(NCHAN)
!          TYPE*,MJ,SUM,SUM2,DIF,ABS(SUM-SUM2)/SUMX
!
60      WRITE(15) (S2(J),J=1,NCHAN)
!
80      CONTINUE
!
        CLOSE(14)
        CLOSE(15)
!
        STOP
100     FORMAT(A)
101     FORMAT(A)
        END
