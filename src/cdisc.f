C------------------------------------------------------------------------------
C Version 28-November-1996                                        file: cdisc.f
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
C Program: cdisc
C Classification: distortion
C Description: Correct C-distortion using the polynomial fits of fitcdis.
C
Comment
C
C Corrige de distorsion C empleando los ajustes polinomicos realizados
C por el programa fitcdis.f
C
        PROGRAM CDISC
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INTEGER READILIM
C
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
        CHARACTER*75 POLFILE,INFILE,OUTFILE
C------------------------------------------------------------------------------
        THISPROGRAM='cdisc'
        CALL WELCOME('28-November-1996')
        CALL SHOWHLP('explanation')
C
        WRITE(*,100)'Arc file name'
        INFILE=INFILEX(14,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
C
        WRITE(*,100)'Spectra from... to... '
        WRITE(CDUMMY,'(I4,A,I4)')1,',',NSCAN
        CALL READ2I(CDUMMY,IBEG,IEND)
        NSPECTRA=IEND-IBEG+1
C
        WRITE(*,101)'(From fitcdis) FIT of line peak deviations vs. '//
     +   'channel number:'
        WRITE(*,100)'Degree of fitted polynomial '
        NORDER=READILIM('@',0,19)
        NTERMS=NORDER+1
C
        WRITE(*,100)'Polynomial file name (from fitcdis)'
        POLFILE=INFILEX(12,'@',0,0,0.,0.,3,.FALSE.)
C
        DO I=IBEG,IEND
          READ(12,112) XJ,(COE(K,I),K=1,NTERMS)
          IF(XJ.NE.FLOAT(I)) STOP 'ERROR IN SCAN NO.'
        END DO
112     FORMAT(1X,F5.0,20E20.9)
        CLOSE(12)
C
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
C
        WRITE(*,100)'Output file name'
        OUTFILE=OUTFILEX(15,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
C
        DO 80 MJ=1,NSCAN
C
          READ(14) (S(J),J=1,NCHAN)
C
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
C
          DO I=1,NCHAN
            IF((X(I)+D(I+1)).GT.0.5) GO TO 10
          END DO
10        I0=I
C
          DO K=1,NCHAN
            IF((X(I0)-D(I0)).LT.(FLOAT(K)+0.5)) GO TO 20
            S2(K)=0.
          END DO
20        K0=K
          DO I=K0,NCHAN
            S2(I)=0.
          END DO
C
          DO 30 I=I0,NCHAN
C
            X1=X(I)-D(I)
            X2=X(I)+D(I+1)
            XT=X2-X1
C
40          FI=FLOAT(K0)
            FI1=FI-0.5
            FI2=FI+0.5
            XMIN=AMAX1(X1,FI1)
            XMAX=AMIN1(X2,FI2)
            FAC=(XMAX-XMIN)/XT
C
            S2(K0)=S2(K0)+FAC*S(I)
C
            IF(FI2.LE.X2) THEN
              K0=K0+1
              IF(FI2.EQ.X2) GO TO 30
              IF(K0.GT.NCHAN) GO TO 50
              GO TO 40
            END IF
30      CONTINUE
50      CONTINUE
C
          SUM=0.
          SUM2=0.
          DO I=1,NCHAN
            SUM=SUM+S(I)
            SUM2=SUM2+S2(I)
          END DO
          IF(SUM.LE.0.) GO TO 60
          DIF=ABS(SUM-SUM2)/SUM
          SUMX=SUM/FLOAT(NCHAN)
C          TYPE*,MJ,SUM,SUM2,DIF,ABS(SUM-SUM2)/SUMX
C
60      WRITE(15) (S2(J),J=1,NCHAN)
C
80      CONTINUE
C
        CLOSE(14)
        CLOSE(15)
C
        STOP
100     FORMAT(A,$)
101     FORMAT(A)
        END
