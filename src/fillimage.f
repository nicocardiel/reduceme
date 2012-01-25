C------------------------------------------------------------------------------
C Version 17-December-1999                                    file: fillimage.f
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
C Program: fillimage
C Classification: arithmetic & manipulations
C Description: Fill an image region by extrapolating the signal in a contigous
C region. The program detects the transition between the extrapolated region
C and the fitted region as a sudden signal change.
C
Comment
C
C
        PROGRAM FILLIMAGE
        IMPLICIT NONE
C
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER READILIM
        REAL READF
C
        INTEGER I,J
        INTEGER J1,J2,JJ
        INTEGER I1,I2,ISTEP
        INTEGER L
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        INTEGER NC1,NC2,NS1,NS2
        INTEGER NC0,NS0
        INTEGER NNC1,NNC2,NNS1,NNS2
        INTEGER NB
        INTEGER NCOLOR
        INTEGER NCHAN_AVE,NS_MIN,NS_MAX,NC_MIN,NC_MAX
        INTEGER NAVERAGE
        REAL A(NCMAX,NSMAX),ERR(NCMAX,NSMAX)
        REAL YCUT(NSMAX)
        REAL XBORDER(NCMAX),YBORDER(NCMAX)
        REAL XC,YC
        REAL ZMIN,ZMAX
        REAL XMIN,XMAX,YMIN,YMAX
        REAL FJUMP,SUM,SUM_ERR,SIGMA_CLIPPING,COEF(2)
        REAL TR(6)
        CHARACTER*1 CERR,CH,CSAVE,CSURE,CDIRECTION,CREPEAT,CFIT
        CHARACTER*50 CDUMMY
        CHARACTER*75 INFILE,ERRFILE,OUTFILE
        LOGICAL LCOLOR(MAX_ID_RED),LOOP
C
        COMMON/BLK1/A
        COMMON/BLK2/NC1,NC2,NS1,NS2
        COMMON/BLK3/ZMIN,ZMAX
        COMMON/BLK4/INFILE
        COMMON/BLK5/TR
        COMMON/BLKLIMITS/XMIN,XMAX,YMIN,YMAX
        COMMON/BLKDEVICE1/NTERM,IDN
C------------------------------------------------------------------------------
        THISPROGRAM='fillimage'
        CALL WELCOME('17-December-1999')
C------------------------------------------------------------------------------
        TR(1)=0.
        TR(2)=1.
        TR(3)=0.
        TR(4)=0.
        TR(5)=0.
        TR(6)=1.
C
        ISTEP=1 !avoid compilation warning
        SUM_ERR=0.0 !avoid compilation warning
C------------------------------------------------------------------------------
        WRITE(*,100)'Work with error images (y/n) '
        CERR(1:1)=READC('n','yn')
C------------------------------------------------------------------------------
        CALL RPGBEGIN(NTERM,IDN,LCOLOR)
C------------------------------------------------------------------------------
        WRITE(*,100)'Input file name'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        DO I=1,NSCAN
          READ(20) (A(J,I),J=1,NCHAN)
        END DO
        CLOSE(20)
C
        IF(CERR.EQ.'y')THEN
          CALL GUESSEF(INFILE,ERRFILE)
          WRITE(*,100)'Error file name '
          ERRFILE=INFILEX(30,ERRFILE,NSCAN,NCHAN,STWV,DISP,21,.TRUE.) !...match
          DO I=1,NSCAN
            READ(30) (ERR(J,I),J=1,NCHAN)
          END DO
          CLOSE(30)
        END IF
C------------------------------------------------------------------------------
        ZMIN=A(1,1)
        ZMAX=ZMIN
        DO I=1,NSCAN
          DO J=1,NCHAN
            IF(A(J,I).LT.ZMIN) ZMIN=A(J,I)
            IF(A(J,I).GT.ZMAX) ZMAX=A(J,I)
          END DO
        END DO
C
        NC1=1
        NC2=NCHAN
        NS1=1
        NS2=NSCAN
C
        CALL PLOT_A
C------------------------------------------------------------------------------
        CALL BUTTON(1,'[z]oom (m)',0)
        CALL BUTTON(2,'[w]hole',0)
        CALL BUTTON(3,'[f]ill',0)
        CALL BUTTON(7,'[s]et BG/FG',0)
        CALL BUTTON(8,'Min[,]Max',0)
        CALL BUTTON(12,'[q]uit&save',0)
C------------------------------------------------------------------------------
10      CONTINUE
        CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
        IF(CH.EQ.'z')THEN
          NB=1
        ELSEIF(CH.EQ.'w')THEN
          NB=2
        ELSEIF(CH.EQ.'f')THEN
          NB=3
        ELSEIF(CH.EQ.'s')THEN
          NB=7
        ELSEIF(CH.EQ.',')THEN
          NB=8
        ELSEIF(CH.EQ.'q')THEN
          NB=12
        ELSE
          CALL IFBUTTON(XC,YC,NB)
        END IF
C------------------------------------------------------------------------------
        IF(NB.EQ.0)THEN
          WRITE(*,100) 'Cursor at '
          WRITE(*,*) NINT(XC),NINT(YC),A(NINT(XC),NINT(YC))
C..............................................................................
        ELSEIF(NB.EQ.1)THEN
          CALL BUTTON(1,'[z]oom (m)',5)
          WRITE(*,100)'Press mouse button...'
          CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
          NNC1=NINT(XC)
          NNS1=NINT(YC)
          IF(NNC1.LT.NC1) NNC1=NC1
          IF(NNC1.GT.NC2) NNC1=NC2
          IF(NNS1.LT.NS1) NNS1=NS1
          IF(NNS1.GT.NS2) NNS1=NS2
          WRITE(*,100)'again...'
          IF(LCOLOR(1)) CALL PGSCI(5)
          CALL RPGBAND(2,0,REAL(NNC1),REAL(NNS1),XC,YC,CH)
          IF(LCOLOR(1)) CALL PGSCI(1)
          NNC2=NINT(XC)
          NNS2=NINT(YC)
          IF(NNC2.LT.NC1) NNC2=NC1
          IF(NNC2.GT.NC2) NNC2=NC2
          IF(NNS2.LT.NS1) NNS2=NS1
          IF(NNS2.GT.NS2) NNS2=NS2
          WRITE(*,101)'OK!'
          IF(NNC1.GT.NNC2)THEN
            NC0=NNC1
            NNC1=NNC2
            NNC2=NC0
          END IF
          IF(NNS1.GT.NNS2)THEN
            NS0=NNS1
            NNS1=NNS2
            NNS2=NS0
          END IF
          NC1=NNC1
          NC2=NNC2
          NS1=NNS1
          NS2=NNS2
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(ITERM.EQ.1)THEN
              CALL RPGERASW(0.,1.,0.,0.65)
            ELSE
              CALL PGPAGE
            END IF
          END DO
          CALL PLOT_A
          CALL BUTTON(1,'[z]oom (m)',0)
C..............................................................................
        ELSEIF(NB.EQ.2)THEN
          CALL BUTTON(2,'[w]hole',5)
          NC1=1
          NC2=NCHAN
          NS1=1
          NS2=NSCAN
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(ITERM.EQ.1)THEN
              CALL RPGERASW(0.,1.,0.,0.65)
            ELSE
              CALL PGPAGE
            END IF
          END DO
          CALL PLOT_A
          CALL BUTTON(2,'[w]hole',0)
C..............................................................................
        ELSEIF(NB.EQ.3)THEN
          CALL BUTTON(3,'[f]ill',5)
          WRITE(*,100) 'No. of channels to average spatial profile'
          NCHAN_AVE=READILIM('@',1,NCHAN)
          WRITE(*,100) 'Channel region to be examined '
          WRITE(CDUMMY,'(I10,A1,I10)') 1,',',NCHAN
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          CALL READ2I(CDUMMY(1:L),NC_MIN,NC_MAX)
          WRITE(*,100) 'Scan region to be examined '
          WRITE(CDUMMY,'(I10,A1,I10)') 1,',',NSCAN
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          CALL READ2I(CDUMMY(1:L),NS_MIN,NS_MAX)
          WRITE(*,100) '[u]pwards or [d]ownwards (u/d) '
          CDIRECTION(1:1)=READC('@','ud')
          NCOLOR=1
          FJUMP=15.
20        WRITE(*,100) 'Percentage diference to detect discontinuity '
          WRITE(CDUMMY,*) FJUMP
          FJUMP=READF(CDUMMY)
          NCOLOR=NCOLOR+1
          IF(NCOLOR.GT.6) NCOLOR=2
          CALL PGSCI(NCOLOR)
          DO J=NC_MIN,NC_MAX
            J1=J-NCHAN_AVE/2
            IF(J1.LT.NC_MIN) J1=NC_MIN
            J2=J+NCHAN_AVE/2
            IF(J2.GT.NC_MAX) J2=NC_MAX
            DO I=NS_MIN,NS_MAX
              YCUT(I)=0.
              DO JJ=J1,J2
                YCUT(I)=YCUT(I)+A(JJ,I)
              END DO
              YCUT(I)=YCUT(I)/REAL(J2-J1+1)
            END DO
            IF(CDIRECTION.EQ.'u')THEN
              I1=NS_MIN
              I2=NS_MAX
              ISTEP=1
            ELSE
              I1=NS_MAX
              I2=NS_MIN
              ISTEP=-1
            END IF
            I=I1-ISTEP
            LOOP=.TRUE.
            DO WHILE(LOOP)
              I=I+ISTEP
              LOOP=((YCUT(I)-YCUT(I+ISTEP))/YCUT(I).LT.FJUMP/100.)
              IF(I.EQ.I2-ISTEP) LOOP=.FALSE.
            END DO
            XBORDER(J-NC_MIN+1)=REAL(J)
            YBORDER(J-NC_MIN+1)=REAL(I)
          END DO
          CALL PGLINE(NC_MAX-NC_MIN+1,XBORDER,YBORDER)
          WRITE(*,100) 'Do you want to repeat with a different '
          WRITE(*,100) 'percentage (y/n) '
          CREPEAT(1:1)=READC('n','yn')
          IF(CREPEAT.EQ.'y') GOTO 20
          WRITE(*,100) 'Are you fitting the computed discontinuity '
          WRITE(*,100) 'with a straight line (y/n) '
          CFIT(1:1)=READC('y','yn')
          IF(CFIT.EQ.'y')THEN
            WRITE(*,100) 'Times sigma to reject points '
            SIGMA_CLIPPING=READF('3.0')
            CALL FITLINE(NC_MAX-NC_MIN+1,XBORDER,YBORDER,
     +       SIGMA_CLIPPING,COEF)
            DO J=NC_MIN,NC_MAX
              YBORDER(J-NC_MIN+1)=COEF(1)+COEF(2)*XBORDER(J-NC_MIN+1)
            END DO
            NCOLOR=NCOLOR+1
            IF(NCOLOR.GT.6) NCOLOR=2
            CALL PGSCI(NCOLOR)
            CALL PGLINE(NC_MAX-NC_MIN+1,XBORDER,YBORDER)
          END IF
          CALL PGSCI(1)
          WRITE(*,*)
          WRITE(*,100) 'No. of pixels in front of the discontinuity '
          WRITE(*,101) 'to be averaged'
          WRITE(*,100) 'in the extrapolation (0=EXIT)'
          NAVERAGE=READILIM('@',0,NSCAN)
          IF(NAVERAGE.GT.0)THEN
            DO J=NC_MIN,NC_MAX
              I1=YBORDER(J-NC_MIN+1)
              IF(CDIRECTION.EQ.'u')THEN
                I2=I1-NAVERAGE+1
                IF(I2.LT.1) I2=1
              ELSE
                I2=I1+NAVERAGE-1
                IF(I2.GT.NSCAN) I2=NSCAN
              END IF
              SUM=0.
              DO I=I1,I2,-ISTEP
                SUM=SUM+A(J,I)
              END DO
              SUM=SUM/REAL(IABS(I2-I1)+1)
              IF(CERR.EQ.'y')THEN
                SUM_ERR=0.
                DO I=I1,I2,-ISTEP
                  SUM_ERR=SUM_ERR+ERR(J,I)
                END DO
                SUM_ERR=SUM_ERR/REAL(IABS(I2-I1)+1)
              END IF
              IF(CDIRECTION.EQ.'u')THEN
                I2=NS_MAX
              ELSE
                I2=NS_MIN
              END IF
              DO I=I1,I2,ISTEP
                A(J,I)=SUM
              END DO
              IF(CERR.EQ.'y')THEN
                DO I=I1,I2,ISTEP
                  ERR(J,I)=SUM_ERR
                END DO
              END IF
            END DO
            CALL PGGRAY(A,NCMAX,NSMAX,NC1,NC2,NS1,NS2,ZMAX,ZMIN,TR)
          END IF
          CALL BUTTON(3,'[f]ill',0)
C..............................................................................
        ELSEIF(NB.EQ.7)THEN
          CALL BUTTON(7,'[s]et BG/FG',5)
          WRITE(CDUMMY,*)ZMIN
          WRITE(*,100)'Background '
          ZMIN=READF(CDUMMY)
          WRITE(CDUMMY,*)ZMAX
          WRITE(*,100)'Foreground '
          ZMAX=READF(CDUMMY)
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(ITERM.EQ.1)THEN
              CALL RPGERASW(0.,1.,0.,0.65)
            ELSE
              CALL PGPAGE
            END IF
          END DO
          CALL PLOT_A
          CALL BUTTON(7,'[s]et BG/FG',0)
C..............................................................................
        ELSEIF(NB.EQ.8)THEN
          CALL BUTTON(8,'Min[,]Max',5)
          ZMIN=A(NC1,NS1)
          ZMAX=ZMIN
          DO I=NS1,NS2
            DO J=NC1,NC2
              IF(A(J,I).LT.ZMIN) ZMIN=A(J,I)
              IF(A(J,I).GT.ZMAX) ZMAX=A(J,I)
            END DO
          END DO
          CALL PLOT_A
          CALL BUTTON(8,'Min[,]Max',0)
C..............................................................................
        ELSEIF(NB.EQ.12)THEN
          CALL BUTTON(12,'[q]uit&save',5)
          WRITE(*,100)'Save output file (y/n)'
          CSAVE(1:1)=READC('@','yn')
          IF(CSAVE.EQ.'y')THEN
            WRITE(*,100)'Output file name'
            OUTFILE=OUTFILEX(30,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
            DO I=1,NSCAN
              WRITE(30) (A(J,I),J=1,NCHAN)
            END DO
            CLOSE(30)
            IF(CERR.EQ.'y')THEN
              WRITE(*,100)'Error file name '
              CALL GUESSEF(OUTFILE,ERRFILE)
              OUTFILE=OUTFILEX(30,ERRFILE,NSCAN,NCHAN,STWV,DISP,
     +         1,.TRUE.)
              DO I=1,NSCAN
                WRITE(30) (ERR(J,I),J=1,NCHAN)
              END DO
              CLOSE(30)
            END IF
          END IF
          WRITE(*,100)'Do you really want to quit (y/n) '
          CSURE(1:1)=READC('n','yn')
          IF(CSURE.EQ.'y')THEN
            CALL PGEND
            STOP
          END IF
          CALL BUTTON(12,'[q]uit&save',0)
C..............................................................................
        END IF
C
        GOTO 10
C------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
        END
C
C******************************************************************************
C
C Dibuja la imagen en escala de grises
        SUBROUTINE PLOT_A
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INTEGER TRUELEN
C
        INTEGER NC1,NC2,NS1,NS2
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        REAL A(NCMAX,NSMAX)
        REAL TR(6),ZMIN,ZMAX
        REAL XMIN,XMAX,YMIN,YMAX
        CHARACTER*75 INFILE
C
        COMMON/BLK1/A
        COMMON/BLK2/NC1,NC2,NS1,NS2
        COMMON/BLK3/ZMIN,ZMAX
        COMMON/BLK5/TR
        COMMON/BLK4/INFILE
        COMMON/BLKLIMITS/XMIN,XMAX,YMIN,YMAX
        COMMON/BLKDEVICE1/NTERM,IDN
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        XMIN=REAL(NC1)-0.5
        XMAX=REAL(NC2)+0.5
        YMIN=REAL(NS1)-0.5
        YMAX=REAL(NS2)+0.5
C
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(ITERM.EQ.1)THEN
            CALL RPGERASW(0.,1.,0.,0.80)
            CALL RPGENV(XMIN,XMAX,YMIN,YMAX,0,-2)
          ELSE
            CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,-2)
          END IF
          CALL PGIDEN_RED
          CALL PGSLW(3)
          CALL PGBOX('BCNTSI',0.0,0,'BCNTSI',0.0,0)
          CALL PGLABEL('channel','scan',CHAR(32))
          CALL PGMTEXT('T',1.0,0.5,0.5,
     +     'File: '//INFILE(1:TRUELEN(INFILE)))
          CALL PGSLW(1)
          CALL PGGRAY(A,NCMAX,NSMAX,NC1,NC2,NS1,NS2,ZMAX,ZMIN,TR)
        END DO
        END
C
C******************************************************************************
C Ajusta una recta de N puntos X(N),Y(N), eliminandos puntos que se alejen
C mas de TSGIMA veces
        SUBROUTINE FITLINE(N,X,Y,TSIGMA,A)
        IMPLICIT NONE
C
        INCLUDE 'redlib.inc'
C
        INTEGER N
        REAL X(N),Y(N)
        REAL TSIGMA,A(2)
C
        INTEGER I
        INTEGER NFIT,K
        REAL XFIT(NCMAX),YFIT(NCMAX)
        REAL CHISQR,RMS
        LOGICAL LFIT1(NCMAX),LFIT2(NCMAX)
        LOGICAL LEXIT
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
C primero ajustamos con todos los puntos
        DO I=1,N
          XFIT(I)=X(I)
          YFIT(I)=Y(I)
          LFIT1(I)=.TRUE.
        END DO
        NFIT=N
C ajustamos la recta
10      CALL POLFIT(XFIT,YFIT,YFIT,NFIT,2,0,A,CHISQR)
C calculamos varianza residual (con los puntos ajustados) y los puntos que se 
C desvian (de toda la muestra)
        RMS=0.
        DO I=1,NFIT
          RMS=RMS+(YFIT(I)-A(1)-A(2)*XFIT(I))*
     +            (YFIT(I)-A(1)-A(2)*XFIT(I))
        END DO
        RMS=SQRT(RMS/REAL(NFIT-1))
        DO I=1,N
          LFIT2(I)=(ABS(Y(I)-A(1)-A(2)*X(I)).LT.TSIGMA*RMS) 
        END DO
C comprobamos si ha cambiado algun punto
        LEXIT=.TRUE.
        DO I=1,N
          IF(LFIT1(I).NEQV.LFIT2(I)) LEXIT=.FALSE.
        END DO
C
        IF(.NOT.LEXIT)THEN
          DO I=1,N
            LFIT1(I)=LFIT2(I)
          END DO
          K=0
          DO I=1,N
            IF(LFIT1(I))THEN
              K=K+1
              XFIT(K)=X(I)
              YFIT(K)=Y(I)
            END IF
          END DO
          NFIT=K
          GOTO 10
        END IF
C
        END
