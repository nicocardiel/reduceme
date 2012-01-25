C------------------------------------------------------------------------------
C Version 10-November-1997                                      file: rebincw.f
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
C Program: rebincw
C Classification: wavelengths
C Description: Performs the wavelength calibration of an image by using the 
C polynomial fits obtained with fitlin.
C
Comment
C
C Programa para calibrar en longitud de onda utilizando los resultados
C de fitlin.
C NOTA: para evitar ampliar el numero de matrices, se ha introducido la
C correccion simultanea de las imagenes de errores haciendo que el programa
C recorra dos veces el mismo codigo. 
C
         PROGRAM REBINCW
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER READI
        INTEGER READILIM
        REAL READF
C
        INTEGER I,J,K,II,KK,JMAX,L
        INTEGER IB0,IB1,J1,J2,MJ,MJMAX
        INTEGER IBEG,IEND,NSPECTRA,NTERMS,NTERMS0,INOPOL
        INTEGER NAR,NARC
        REAL S(NCMAX),SA(NCMAX,NSMAX),SB(NCMAX,NSMAX),SS(NCMAX)
        REAL ERR(NCMAX,NSMAX)
        REAL COEI(2,20,NSMAX),COE1(20),AI(2,20),COE(20,NSMAX)
        REAL XJ,A(20),ARCN
        REAL STWV0,DISP0,WAV,SUM,XA,XB,X1,POLA,POLB,X1M,X2M
        REAL F1,F2,A1,A2,A0,PIXELM,ECH
        REAL ST,ST2
        REAL MAXERROR
        REAL RADVEL
        DOUBLE PRECISION PIXEL,PIXEL0,PIXLAM,Q,XACC,P,ERW
        CHARACTER*1 CSHIFT,CCDIS,CERR
        CHARACTER*50 CDUMMY
        CHARACTER*75 INFILE,POLFILE0,POLFILE,OUTFILE,ERRFILE
        CHARACTER*79 SEPARADOR
        LOGICAL LOGFILE,LOGFILE2
        LOGICAL LFIRST_PASS             !determina si es la primera pasada o no
C
        COMMON/BLKI/NTERMS,JMAX
        COMMON/BLKF/A
        COMMON/BLKD/XACC
C
        COMMON/BLK1/S,SS
        COMMON/BLK2/NCHAN,STWV,DISP
        COMMON/BLK3/RADVEL
C------------------------------------------------------------------------------
        THISPROGRAM='rebincw'
        CALL WELCOME('8-December-1996')
C------------------------------------------------------------------------------
        LFIRST_PASS=.TRUE.
C
        DO I=1,79
          SEPARADOR(I:I)='-'
        END DO
C
        WRITE(*,*)
        WRITE(*,101)SEPARADOR
        WRITE(*,101)'                               Program REBINCW'
        WRITE(*,101)'    This program needs data (files) generated '//
     +   'through FITCDIS and FITLIN.'
        WRITE(*,101)'    Since these files are suposed to be in a d'//
     +   'irectory containing all the'
        WRITE(*,101)'    files from arc reduction, REBINCW will sea'//
     +   'rch for the output files of'
        WRITE(*,101)'    FITCIDS and FITLIN  in the current directo'//
     +   'ry first, and,  if they are'
        WRITE(*,101)'    not present there, the search will continu'//
     +   'e in the directory ../cuar/.'
        WRITE(*,101)SEPARADOR
        WRITE(*,*)
C------------------------------------------------------------------------------
        WRITE(*,100)'Work with error images (y/n) '
        CERR(1:1)=READC('n','yn')
C------------------------------------------------------------------------------
C Leemos la imagen a calibrar
        WRITE(*,100)'Input file name'
        INFILE=INFILEX(14,'@',NSCAN,NCHAN,STWV0,DISP0,1,.FALSE.)
        DO I=1,NSCAN
          READ(14) (S(J),J=1,NCHAN)
          DO J=1,NCHAN
            SA(J,I)=S(J)
            SB(J,I)=0.0
          END DO
        END DO
        CLOSE(14)
        IF(CERR.EQ.'y')THEN
          WRITE(*,100)'Input error file name '
          CALL GUESSEF(INFILE,ERRFILE)
          ERRFILE(1:75)=
     +     INFILEX(20,ERRFILE,NSCAN,NCHAN,STWV0,DISP0,21,.TRUE.)  !match
          DO I=1,NSCAN
            READ(20) (ERR(J,I),J=1,NCHAN)
          END DO
          CLOSE(20)
        END IF
C------------------------------------------------------------------------------
C Definimos si vamos a corregir de distorsion C
        WRITE(*,*)
        WRITE(*,100)'Correct from C-distortion (y/n) '
        CCDIS(1:1)=READC('y','yn')
        IF(CCDIS.EQ.'n')THEN
          IBEG=1
          IEND=NSCAN
          NSPECTRA=IEND-IBEG+1
          NTERMS0=1
          INOPOL=0
        END IF
C------------------------------------------------------------------------------
C Introducimos numero de arcos que van a utilizarse
        WRITE(*,100)'Number of arcs to be used '
        NARC=READILIM('1',1,2)
C------------------------------------------------------------------------------
C Bucle en numero de arcos
        DO NAR=1,NARC
          WRITE(*,*)
          WRITE(*,110)'ARC #',NAR
          IF(CCDIS.EQ.'n') GOTO 29
C..............................................................................
          WRITE(*,*)
          WRITE(*,101)'>>> C-DISTORTION CORRECTION:'
          IF(NAR.EQ.1) THEN
            WRITE(*,100)'Spectra for which the corr. was derived: '//
     +       'from... to... '
            WRITE(CDUMMY,'(I1,A1,I8)')1,',',NSCAN
            CALL RMBLANK(CDUMMY,CDUMMY,L)
            CALL READ2I(CDUMMY(1:L),IBEG,IEND)
            NSPECTRA=IEND-IBEG+1
            WRITE(*,100)'Order of polynomial '
            NTERMS0=READILIM('@',0,19)
            NTERMS0=NTERMS0+1
          END IF
          WRITE(*,100)'How many polynomial to be applied '//
     +     '(1=no iteration) '
          INOPOL=READI('1')
          DO I=1,NTERMS0
            DO J=IBEG,IEND
              COEI(NAR,I,J)=0.
            END DO
          END DO
          DO II=1,INOPOL
            WRITE(*,'(A,I2,A,$)')'Iteration #',II,': '
20          WRITE(*,100)'Polynomial file name (from fitcdis)'
            POLFILE0(1:75)=READC('@','@')
            INQUIRE(FILE=POLFILE0,EXIST=LOGFILE)
            IF(LOGFILE)THEN
              OPEN(12,FILE=POLFILE0,STATUS='OLD',FORM='FORMATTED')
            ELSE 
              INQUIRE(FILE='../cuar/'//POLFILE0,EXIST=LOGFILE2)
              IF(LOGFILE2)THEN
                OPEN(12,FILE='../cuar/'//POLFILE0,STATUS='OLD',
     +           FORM='FORMATTED')
              ELSE
                WRITE(*,101)'ERROR: this file does not exist.'//
     +           ' Try again.'
                GOTO 20
              END IF
            END IF
            DO I=IBEG,IEND
              READ(12,112) XJ,(COE1(K),K=1,NTERMS0)
              IF(INT(XJ).NE.I) STOP 'FATAL ERROR: in scan no.'
              DO K=1,NTERMS0
                COEI(NAR,K,I)=COEI(NAR,K,I)+COE1(K)
              END DO
            END DO
112         FORMAT(1X,F5.0,20E20.9)
            CLOSE(12)
          END DO
C..............................................................................
29        WRITE(*,*)
          WRITE(*,101)'>>> WAVELENGTH CALIBRATION:'
!         WRITE(*,100)'Order of polynomial '
!         NTERMS=READILIM('@',0,19)
!         NTERMS=NTERMS+1
30        WRITE(*,100)'Polynomial file name (from fitlin)'
          POLFILE(1:75)=READC('@','@')
          INQUIRE(FILE=POLFILE,EXIST=LOGFILE)
          IF(LOGFILE)THEN
            OPEN(15,FILE=POLFILE,STATUS='OLD',FORM='FORMATTED')
          ELSE
            INQUIRE(FILE='../cuar/'//POLFILE,EXIST=LOGFILE2)
            IF(LOGFILE2)THEN
              OPEN(15,FILE='../cuar/'//POLFILE,STATUS='OLD',
     +         FORM='FORMATTED')
            ELSE
              WRITE(*,101)'ERROR: this file does not exist. Try again.'
              GOTO 30
            END IF
          END IF
          K=1
31        READ(15,*,END=32) KK,AI(NAR,K)
          K=K+1
          GOTO 31
32        CLOSE(15)
          NTERMS=K-1
          WRITE(*,100) '>>> Polynomial degree: '
          WRITE(*,*) NTERMS-1
        END DO
        WRITE(*,*)
C------------------------------------------------------------------------------
C Promediamos los coeficientes de distorsion C y calibracion en l.d.o.
        ARCN=FLOAT(NARC)
C..............................................................................
        DO I=1,NTERMS0
          DO J=IBEG,IEND
            COE(I,J)=0.                      !coeficientes para la distorsion C
          END DO
        END DO
        IF(CCDIS.EQ.'y')THEN
          DO NAR=1,NARC
            DO I=1,NTERMS0
              DO J=IBEG,IEND
                COE(I,J)=COE(I,J)+COEI(NAR,I,J)
              END DO
            END DO
          END DO
          DO I=1,NTERMS0
            DO J=IBEG,IEND
              COE(I,J)=COE(I,J)/ARCN
            END DO
          END DO
        END IF
C..............................................................................
        DO K=1,NTERMS
          A(K)=0.              !coeficientes polinomio de calibracion en l.d.o.
        END DO
        DO NAR=1,NARC
          DO K=1,NTERMS
            A(K)=A(K)+AI(NAR,K)
          END DO
        END DO
        DO K=1,NTERMS
          A(K)=A(K)/ARCN
        END DO
C------------------------------------------------------------------------------
C Si no hemos usado todos los scans para determinar la distorsion C,
C extrapolamos los coeficientes por debajo de IBEG y por encima de IEND
        IF(CCDIS.EQ.'y')THEN
          DO I=1,IBEG-1
            DO K=1,NTERMS0
              COE(K,I)=COE(K,IBEG)
            END DO
          END DO
          DO I=IEND+1,NSCAN
            DO K=1,NTERMS0
              COE(K,I)=COE(K,IEND)
            END DO
          END DO
        END IF
C------------------------------------------------------------------------------
        WRITE(*,*)
        WRITE(*,100)'Enter STWV (central wavelenght of pixel 1)'
        STWV=READF('@')
        WRITE(*,100)'Enter desired dispersion (angstroms/pixel)'
        DISP=READF('@')
        WRITE(*,*)
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
111        XACC=1.D-10
        JMAX=200
C
C   FROM VISTA:
C        To rebin the old spectrum, we calculate the extreme lambdas of each
C        bin in the new spectrum. The location of these lambdas in pixel space
C        of the original spectrum is calculated using the inverse dispersion.
C        Whole pixels lying inside the new bin are directly added. Anyfractional
C        pixel contribution to the bin is estimated by integration of the
C        parabola whose area matches the intensity of the pixel and its two
C        neighbors (Gonzalez-integration).
C        NOTE: this is different from Sympson-integration which assumes the
C        function is only sampled (not binned and sampled). The Sympson-parabola
C        C0+C1*X+C2*X*X passing through (-1,Ym),(0,Y0),(1,Yp), has coefficients
C        C0=Y0, C1=(Yp-Ym)/2, and C2=(Yp+Ym)/2-Y0, while the Gonzalez-parabola
C        A0+A1*X+A2*X*X, has A0=C0-(C2/12), A1=C1, and A2=C2. This scheme takes
C        care of the fact that the original spectrum is also binned data.
C
        DO I=1,NCHAN
          IB0=0
          IB1=0
C-- Find the limits of the new-pixel in the old pixel scale:
          IF (I.GT.1) THEN
            PIXEL=PIXEL0
          ELSE
            WAV = STWV + (FLOAT(I)-1.5)*DISP
            PIXEL = PIXLAM(WAV,0.D0,I)
          END IF
          X1 = SNGL(PIXEL)
          PIXEL0=PIXEL
          WAV = STWV + (FLOAT(I)-0.5)*DISP
          PIXEL = PIXLAM(WAV,PIXEL0,I)
          Q=PIXEL
          P=DBLE(A(NTERMS))
          DO J=NTERMS-1,1,-1
            P=P*Q+DBLE(A(J))
          END DO
          ERW=DABS(DBLE(WAV)-P)
          IF(ERW.GT.5.D-6) THEN
            WRITE(*,*)'Warning: Error in wavelength assignations'
            WRITE(*,'(A,I5,A,$)')'Pixel ',I,'    Diff='
            WRITE(*,*)ERW
            STOP
          END IF
          PIXEL0=PIXEL
          DO MJ=1,NSCAN
            SUM=0.0
            XA=X1
            XB=SNGL(PIXEL)  
            POLA=COE(NTERMS0,MJ)
            POLB=POLA
            DO J=NTERMS0-1,1,-1
              POLA=POLA*XA+COE(J,MJ)
              POLB=POLB*XB+COE(J,MJ)
            END DO
C            TYPE*,I,MJ,POLA,POLB
            X1M=X1+POLA
            PIXELM=SNGL(PIXEL)+POLB
            X2M = MAX(0.5,MIN(NCHAN+0.5,MAX(PIXELM,X1M)))
            X1M = MAX(0.5,MIN(NCHAN+0.5,MIN(PIXELM,X1M)))
            IF (X1M.EQ.X2M) THEN                ! The new pixel maps outside the
              SB(I,MJ) = 0.0E0                  ! range of the old spectrum.
              IF(IB0.EQ.0) IB0=MJ
              IB1=MJ
              GO TO 200
            END IF
            J1 = MIN0(NCHAN-1,MAX0(NINT(X1M),2))    ! This takes care of the
            J2 = MAX0(2,MIN0(NINT(X2M),NCHAN-1))    ! interpolation near edges.
            F1 = X1M - FLOAT(J1)
            F2 = X2M - FLOAT(J2)
            A1 = (SA(J1+1,MJ)-SA(J1-1,MJ))/4.
            A2 = (SA(J1+1,MJ)+SA(J1-1,MJ))/6. - SA(J1,MJ)/3.
            IF (J2 .EQ. J1) THEN
C              A0 = SA(J1,MJ) + A2/4.0                  ! Binned parabola.
              A0 = SA(J1,MJ) - A2/4.
              SUM = (F2-F1)*(A0+(F1+F2)*A1+A2*(F1*F1+F2*F1+F2*F2))
            ELSE
C-- Left-pixel fractional contribution to flux:
              SUM = (0.5-F1)*(SA(J1,MJ)+(F1+0.5)*(A1+A2*F1))
C-- Whole pixel contribution to flux:
              IF (J1.LE.J2-2) THEN                      ! This check was put
                DO J=J1+1,J2-1,1                        ! here for portability.
                  SUM = SUM + SA(J,MJ)
                END DO
              END IF
C-- Right-pixel fractional contribution to flux:
              A1 = (SA(J2+1,MJ)-SA(J2-1,MJ))/4.
              A2 = (SA(J2+1,MJ)+SA(J2-1,MJ))/6. - SA(J2,MJ)/3.
              SUM = SUM + (0.5+F2)*(SA(J2,MJ)+(F2-0.5)*(A1+A2*F2))
            END IF
          
            SB(I,MJ)=SUM
           END DO
200        CONTINUE
C          SPA(I) = SUM/(X2-X1)         ! Mean flux inside the new bin.
          IF(IB0.NE.0) WRITE(*,'(3(A,I5))')'Pixel ',I,'    -------> '//
     +     'Blank from scan ',IB0,' to ',IB1
        END DO
C------------------------------------------------------------------------------
C Conservacion del numero de cuentas
        WRITE(*,*)
        WRITE(*,100)'>>> Checking conservation of number of counts...'
        MAXERROR=0.
        DO MJ=1,NSCAN
          ST=0.
          ST2=0.
          DO I=1,NCHAN
            ST=ST+SA(I,MJ)
            ST2=ST2+SB(I,MJ)
          END DO
          IF(ST.EQ.0) THEN
            IF(ST2.EQ.0) GOTO 456
            WRITE(*,*)
            WRITE(*,110)'Error in scan #',MJ
            WRITE(*,*)ST,ST2
            GOTO 456
          END IF
          ECH=ABS(ST-ST2)/ST*100
          IF(ECH.GT.0.5) THEN
ccc            WRITE(*,*)
ccc            WRITE(*,*)MJ,ST,ST2,ABS(ST-ST2)
ccc            IF(ST.GT.0.) WRITE(*,*)'ERROR = ',ABS(ST-ST2)/ST*100.,'  %'
            IF(ABS(ST-ST2)/ST*100. .GT. MAXERROR)THEN 
              MAXERROR=ABS(ST-ST2)/ST*100.
              MJMAX=MJ
            END IF
          END IF
456       CONTINUE
        END DO
        WRITE(*,101)'   ...OK!'
        WRITE(*,100)'>>> Maximum error (%): '
        WRITE(*,*)MAXERROR
        WRITE(*,110)'    in scan no.: ',MJMAX
        WRITE(*,*)
C------------------------------------------------------------------------------
C Como opcion, podemos trasladar los espectros una cierta velocidad,
C manteniendolos en escala lineal y con la misma STWV y DISP.
        IF(LFIRST_PASS)THEN
          WRITE(*,100)'Shift spectra (radial velocitiy) (y/n) '
          CSHIFT(1:1)=READC('n','yn')
        END IF
        IF(CSHIFT.EQ.'y')THEN
          IF(LFIRST_PASS)THEN
            WRITE(*,101)'(Note: - Vr in file from program rvel)'
            WRITE(*,101)'       ^ (minus)'
            WRITE(*,100)'Radial velocity correction (km/sec)'
            RADVEL=READF('@')
          END IF
          DO MJ=1,NSCAN
            DO I=1,NCHAN
              S(I)=SB(I,MJ)
            END DO
            CALL RVREBIN(RADVEL,NCHAN,S,SS,STWV,DISP)
            DO I=1,NCHAN
              SB(I,MJ)=SS(I)
            END DO
          END DO
        END IF
C------------------------------------------------------------------------------
        IF(CERR.EQ.'y')THEN
          IF(LFIRST_PASS)THEN
            WRITE(*,*)
            WRITE(*,101)'Working with errors...'
            DO MJ=1,NSCAN
              DO I=1,NCHAN
                SA(I,MJ)=SB(I,MJ)
                SB(I,MJ)=ERR(I,MJ)
                ERR(I,MJ)=SA(I,MJ)
                SA(I,MJ)=SB(I,MJ)
                SB(I,MJ)=0.
              END DO
            END DO
            LFIRST_PASS=.FALSE.
            GOTO 111
          ELSE
            DO MJ=1,NSCAN
              DO I=1,NCHAN
                SA(I,MJ)=ERR(I,MJ)
                ERR(I,MJ)=SB(I,MJ)
              END DO
            END DO
          END IF
        ELSE
          DO MJ=1,NSCAN
            DO I=1,NCHAN
              SA(I,MJ)=SB(I,MJ)
            END DO
          END DO
        END IF
C------------------------------------------------------------------------------
        WRITE(*,*)
        WRITE(*,100)'Output file name'
        OUTFILE=OUTFILEX(12,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        DO MJ=1,NSCAN
          DO I=1,NCHAN
            S(I)=SA(I,MJ)
          END DO
          WRITE(12) (S(I),I=1,NCHAN)
        END DO
        CLOSE(12)
        IF(CERR.EQ.'y')THEN
          WRITE(*,100)'Output error file name '
          CALL GUESSEF(OUTFILE,ERRFILE)
          ERRFILE=OUTFILEX(13,ERRFILE,NSCAN,NCHAN,STWV,DISP,1,.TRUE.)
          DO I=1,NSCAN
            WRITE(13) (ERR(J,I),J=1,NCHAN)
          END DO
          CLOSE(13)
        END IF
C
        STOP
100     FORMAT(A,$)
101     FORMAT(A)
110     FORMAT(A,I6)
        END
C
C******************************************************************************
C
        DOUBLE PRECISION FUNCTION PIXLAM(WAV,X0,II)
        IMPLICIT NONE
        INTEGER II,J,JMAX,NTERMS
        REAL WAV,A(20)
        DOUBLE PRECISION XP,X1,X2,PRO,PRO2,DPRO,X0,XACC,RTNEWT,F,DF
        DOUBLE PRECISION DX,DX0
        COMMON/BLKI/NTERMS,JMAX
        COMMON/BLKF/A
        COMMON/BLKD/XACC
C------------------------------------------------------------------------------
        pixlam=0.d0    !ncl: evita un Warning que no tiene sentido
        XP=X0
10      CALL FUNCD(XP,PRO,DPRO,WAV)
        IF(PRO.GT.0.D00) THEN
          XP=XP-1.D00
          GO TO 10
        ELSE
          X1=XP
        END IF
        X2=XP
20      X2=X2+2.D00
        CALL FUNCD(X2,PRO2,DPRO,WAV)
        IF(PRO2.LE.0.D00) GO TO 20
C
C        TYPE*,X1,X2,PRO,PRO2
        RTNEWT=.5D00*(X1+X2)
        DO J=1,JMAX
          CALL FUNCD(RTNEWT,F,DF,WAV)
          DX=F/DF
          IF(J.GT.100) WRITE(*,*)'DX',DX,DX0,DX+DX0
          IF(DABS(DX+DX0).LT.(XACC*1.D-10)) THEN
            WRITE(*,*)II
            WRITE(*,*)'DX',DX,DX0,DX+DX0
            DX=DX/2.D00
            WRITE(*,*)'OOOHO'
          END IF
          DX0=DX
C          TYPE*,RTNEWT,F,DX
          RTNEWT=RTNEWT-DX
          IF((X1-RTNEWT)*(RTNEWT-X2).LT.0.)THEN
            WRITE(*,*)PRO,PRO2
            WRITE(*,*)X1,X2,RTNEWT,WAV
ccc            PAUSE 'WARNING: Jumped out of brackets'
            WRITE(*,'(A)') 'WARNING: Jumped out of brackets'
          END IF
          IF(ABS(DX).LT.XACC) THEN
            PIXLAM=RTNEWT
            RETURN
          END IF
        END DO
        WRITE(*,*)X1,X2,RTNEWT,WAV
        STOP 'PIXLAM exceeding maximum iterations'
        END
C
C******************************************************************************
C
        SUBROUTINE FUNCD(X,F,DF,WAV)
        IMPLICIT NONE
        INTEGER J,NTERMS,JMAX
        REAL WAV,A(20)
        DOUBLE PRECISION X,F,DF,P,DP,XACC
        COMMON/BLKI/NTERMS,JMAX
        COMMON/BLKF/A
        COMMON/BLKD/XACC
C------------------------------------------------------------------------------
        P=DBLE(A(NTERMS))
        DP=0.D00
        DO J=NTERMS-1,1,-1
          DP=DP*X+P
          P=P*X+DBLE(A(J))
        END DO
        F=P-DBLE(WAV)
        DF=DP
        RETURN
        END
