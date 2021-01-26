C------------------------------------------------------------------------------
C Version 6-December-1996                                        file: fitlin.f
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
C Program: fitlin
C Classification: wavelengths
C Description: Calculates the wavelength calibration polynomials.
C
Comment
C
        PROGRAM FITLIN
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER TRUELEN
        INTEGER READI
        INTEGER READILIM
        REAL READF
C
        INTEGER MAXLIN
        PARAMETER(MAXLIN=1000)
C
        INTEGER I,J,K,J1,J2,IJ,IN,II,IJK
        INTEGER NLINES,MLINES,M,M0,NLFIX,MM
        INTEGER LINE(MAXLIN),LINEFIX(MAXLIN)
        INTEGER IEX(MAXLIN),LINEM(4),LINEN(MAXLIN)
        INTEGER ICEN,ISEARCH
        INTEGER IDM,IPARA,IMIN,IMAX,IPARAL,IPARAR,KEYL,KEYR
        INTEGER NPT,NDIM
        INTEGER IPLOT,ILINEQ
        INTEGER NORDER,NTERMS,NPTS,MODE
        INTEGER LLINEM,NOCH
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        INTEGER NEVAL
        REAL S(NCMAX),YMAX,YMIN,RANGE
        REAL FLAM(MAXLIN),FACTOR_FLAM
        REAL DEVMEAN,FML
        REAL YRMSTOL
        REAL XX0(3),DXX0(3),XXF(3),DXXF(3)
        REAL PAR(3,MAXLIN),XPOS(MAXLIN)
        REAL XPOL(MAXLIN),YPOL(MAXLIN),SIGMAY(MAXLIN)
        REAL A(20),B(2),CHISQR
        REAL DEVM(4),DEVI(MAXLIN),DEVI2(MAXLIN)
        REAL POL,SIGMA,XF,ADEVI,TDEVM
        REAL YTE,XXX(NCMAX),YYY(NCMAX)
        REAL SX(9),SY(9)
        EXTERNAL FUNKFITLIN
        REAL FUNKFITLIN
        CHARACTER*1 CREVER,CITER,CRET
        CHARACTER*5 STRIN
        CHARACTER*20 TTER
        CHARACTER*50 CDUMMY
        CHARACTER*75 ARCFILE,TABFILE,POLFILE
        LOGICAL LCOLOR(MAX_ID_RED)
C
        COMMON/BLKS/S
        COMMON/BLKNCHAN/NCHAN
        COMMON/BLKI/IMIN,IMAX
        COMMON/BLKARC/ARCFILE
        COMMON/BLKF1/SX,SY
        COMMON/BLKF2/NPT
        COMMON/BLKP1/LINE,NLINES
        COMMON/BLKP2/PAR,XPOS
        COMMON/BLKTER/TTER
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
C------------------------------------------------------------------------------
        THISPROGRAM='fitlin'
        CALL WELCOME('6-December-1996')
C
        NORDER=0
        YRMSTOL=1.E-6
C
        CALL PIDEGTER(NTERM,IDN,LCOLOR)
C
        WRITE(*,100)'Arc file name'
        ARCFILE=INFILEX(14,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        IF(TRUELEN(OBJECT).GT.0)THEN
          ARCFILE=ARCFILE(1:TRUELEN(ARCFILE))//' ['//
     +     OBJECT(1:TRUELEN(OBJECT))//']'
        END IF
        READ(14) (S(I),I=1,NCHAN)
        CLOSE(14)
C
        CONTINUE
        WRITE(*,100)'File name with line positions and wavelengths '
        TABFILE=INFILEX(12,'fitlin.dat',0,0,.0,.0,3,.FALSE.)
        WRITE(*,100) 'Factor to be applied to wavelength'
        FACTOR_FLAM=READF('1.0')
        I=0
202     CONTINUE
        I=I+1
        READ(12,*,END=203)LINE(I),FLAM(I)
        FLAM(I)=FLAM(I)*FACTOR_FLAM
        GOTO 202
203     CONTINUE
        CLOSE(12)
        NLINES=I-1
        WRITE(*,110)'No. of lines read: ',NLINES
C
        WRITE(*,100)'Do you want to reverse the line positions (y/n) '
        CREVER(1:1)=READC('n','yn')
        IF(CREVER.EQ.'y')THEN
          DO I=1,NLINES
            LINE(I)=NCHAN-LINE(I)+1
          END DO
        END IF
C
        WRITE(*,100)'Centroid (Line=Line+Centroid) '
        ICEN=READI('0')
        IF (ICEN.NE.0) THEN
          DO I=1,NLINES
            LINE(I)=LINE(I)+ICEN
          END DO
        END IF
C
        WRITE(*,100)'Increment to search for peaks (0=no search) '
        ISEARCH=READILIM('2',0,9999)
        IF(ISEARCH.EQ.0) GO TO 204
C
        IDM=0
        DEVMEAN=0.
        DO I=1,NLINES
          FML=0.
          J1=LINE(I)-ISEARCH
          IF(J1.LT.1)J1=1
          J2=LINE(I)+ISEARCH
          IF(J2.GT.NCHAN)J2=NCHAN
          DO J=J1,J2
            IF(S(J).GT.FML) THEN
              IJ=J
              FML=S(J)
            END IF
          END DO
          IF(ABS(LINE(I)-IJ).GT.IDM) IDM=ABS(LINE(I)-IJ)
          WRITE(*,*)I,FLAM(I),LINE(I),IJ,IJ-LINE(I)
          DEVMEAN=DEVMEAN+REAL(IJ-LINE(I))
          LINE(I)=IJ
        END DO
        DEVMEAN=DEVMEAN/REAL(NLINES)
        WRITE(*,*)
        WRITE(*,110)'>>> Max. deviation: ',IDM
        WRITE(*,100)'>>> Mean: '
        WRITE(*,*)DEVMEAN
        WRITE(*,*)
C
204     WRITE(CDUMMY,*)YRMSTOL
        WRITE(*,100)'YRMSTOL for DOWNHILL '
        YRMSTOL=READF(CDUMMY)
C
        WRITE(*,100)'How many points in each side of gaussian '//
     +   '(maximum = 4) '
        IPARA=READI('3')
        DO I=1,NLINES
          IMAX=LINE(I)
          YMAX=S(IMAX)
C Determine the number of points available in each side of gaussian
          IPARAL=0
          IPARAR=0
          KEYL=0
          KEYR=0
          DO IN=1,IPARA
            IF(KEYL.EQ.0)THEN
              IF(S(IMAX-IN).LT.S(IMAX-IN+1))THEN
                IPARAL=IPARAL+1
              ELSE
                KEYL=1
              END IF
            END IF
            IF(KEYR.EQ.0)THEN
              IF(S(IMAX+IN).LT.S(IMAX+IN-1))THEN
                IPARAR=IPARAR+1
              ELSE
                KEYR=1
              END IF
            END IF
          END DO
          NPT=IPARAL+IPARAR+1
          WRITE(*,*)I,IMAX,IPARAL,IPARAR,NPT
          IF(NPT.LT.4)THEN
            WRITE(*,110)'No. points is not enough to fit the line #',I
            WRITE(*,'(A,I3,A)')'Fit forced with ',IPARA,
     +       ' points at both sides of gaussian.'
            IPARAL=IPARA
            IPARAR=IPARA
            NPT=IPARAL+IPARAR+1
          END IF  
C
          DO IN=1,NPT
            II=IMAX-IPARAL+IN-1
            SX(IN)=REAL(II-IMAX)
            SY(IN)=S(II)/YMAX
          END DO
C
          NDIM=3
          XX0(1)=1.0
          XX0(2)=0.
          XX0(3)=2.0
          DXX0(1)=0.1
          DXX0(2)=0.1
          DXX0(3)=0.1
          CALL DOWNHILL(NDIM,XX0,DXX0,FUNKFITLIN,1.0,0.5,2.0,YRMSTOL,
     +     XXF,DXXF,NEVAL)
          IF(NEVAL.GE.2000)THEN
            WRITE(*,110)'DOWNHILL exceeds maximum number of '//
     +       'iterations in line #',I
          END IF
          DO J=1,3
            PAR(J,I)=XXF(J)
          END DO
          IF(ABS(PAR(2,I)).GT.3.) WRITE(*,110)'Bad fit fot line #',I
          XPOS(I)=PAR(2,I)+REAL(LINE(I))
C
        END DO
C
        MLINES=NLINES
        M0=1
        NLFIX=0
205     WRITE(*,100)'Plot fits (-1=no, 0=range, I=line I) '
        IPLOT=READI('-1')
        IF(IPLOT.GE.0) THEN
          IF(IPLOT.EQ.0) THEN
            WRITE(*,100)'Range (channels)'
            CALL READ2I('@',IMIN,IMAX)
            IF(IMIN.LT.1) IMIN=1-10
            IF(IMAX.GT.NCHAN) IMAX=NCHAN+10
          ELSE
            IF((IPLOT.LT.0).OR.(IPLOT.GT.NLINES))THEN
              WRITE(*,101)'ERROR: number out of range. Try again.'
              GOTO 205
            END IF
            IMIN=LINE(IPLOT)-10
            IMAX=LINE(IPLOT)+10
          END IF
          CALL PLOTCH
          GO TO 205
        END IF
C
206     CONTINUE
        WRITE(*,100)'Fixed lines (0=exit) '
        ILINEQ=READILIM('0',0,NLINES)
        IF(ILINEQ.NE.0)THEN
          WRITE(*,'(A,F10.2,5X,A,2X,I5)')'Old position: ',
     +         XPOS(ILINEQ),'New position: ',LINE(ILINEQ)
          XPOS(ILINEQ)=REAL(LINE(ILINEQ))
          NLFIX=NLFIX+1
          LINEFIX(NLFIX)=ILINEQ
          GOTO 206
        END IF
        IF(NLFIX.GT.0)THEN
          WRITE(*,110)'===> Lines fixed: ',NLFIX
          DO IJK=1,NLFIX
            WRITE(*,*)LINEFIX(IJK)
          END DO
        END IF
C
        M=M0
207     WRITE(*,100)'Excluded lines (0=exit) '
        IEX(M)=READI('0')
        IF(IEX(M).NE.0) THEN
          IF(M-1.GT.0)THEN
            DO K=1,M-1
              IF(IEX(K).EQ.IEX(M))THEN
                WRITE(*,101)'ERROR: this line has been already '//
     +           'excluded'
                GOTO 207
              END IF
            END DO
          END IF
          M=M+1
          GO TO 207
        END IF
        M0=M
        M=M-1
        MLINES=NLINES-M
C
        IF(NLINES-MLINES.GT.0)THEN
          WRITE(*,110)'===> Excluded lines: ',NLINES-MLINES
          DO I=1,NLINES-MLINES
            WRITE(*,*)IEX(I)
          END DO
        END IF
C
        K=1
        MM=0
        DO 209 I=1,NLINES
          IF(M.EQ.0) GO TO 208
          DO K=1,M
            IF(I.EQ.IEX(K)) GO TO 209
          END DO
208       MM=MM+1
          LINEN(MM)=I
          XPOL(MM)=XPOS(I)
          YPOL(MM)=FLAM(I)
          SIGMAY(MM)=1.
209     CONTINUE
        IF(MM.NE.MLINES) STOP 'FATAL ERROR: in No. of excluded lines'
C
        WRITE(*,100)'Polynomial degree '
        IF(NORDER.EQ.0)THEN
          NORDER=READILIM('@',0,19)
        ELSE
          WRITE(CDUMMY,*)NORDER
          NORDER=READILIM(CDUMMY,0,19)
        END IF
        NTERMS=NORDER+1
        NPTS=MLINES
        MODE=0
        CALL POLFIT(XPOL,YPOL,SIGMAY,NPTS,NTERMS,MODE,A,CHISQR)
        CALL POLFIT(XPOL,YPOL,SIGMAY,NPTS,2,MODE,B,CHISQR)
C
        WRITE(*,101)'>>> Coefficients of linear approximation '//
     +   '(STWV,DISP):'
        WRITE(*,*)B(1),B(2)
        WRITE(*,*)
        WRITE(*,101)'>>> Coefficients of fitted polynomial:'
        DO K=1,NTERMS
          WRITE(*,*)K,A(K)
        END DO
C
        DO IJ=1,4
         DEVM(IJ)=0.
         LINEM(IJ)=0
        END DO

        SIGMA=0.
        DO J=1,MLINES
          POL=A(NTERMS)
          DO K=NTERMS-1,1,-1
            POL=POL*XPOL(J)+A(K)
          END DO
          DEVI(J)=YPOL(J)-POL
          XF=XPOL(J)
          DEVI2(J)=YPOL(J)-B(1)-XF*B(2)
C looking for the first four lines with the largest deviations
          ADEVI=ABS(DEVI(J))
          IF(ADEVI.GT.DEVM(4)) THEN
            DEVM(4)=ADEVI
            LINEM(4)=LINEN(J)
            IF(DEVM(4).GT.DEVM(3)) THEN
              TDEVM=DEVM(3)
              LLINEM=LINEM(3)
              DEVM(3)=DEVM(4)
              LINEM(3)=LINEM(4)
              DEVM(4)=TDEVM
              LINEM(4)=LLINEM
              IF(DEVM(3).GT.DEVM(2)) THEN
                TDEVM=DEVM(2)
                LLINEM=LINEM(2)
                DEVM(2)=DEVM(3)
                LINEM(2)=LINEM(3)
                DEVM(3)=TDEVM
                LINEM(3)=LLINEM
                IF(DEVM(2).GT.DEVM(1)) THEN
                  TDEVM=DEVM(1)
                  LLINEM=LINEM(1)
                  DEVM(1)=DEVM(2)
                  LINEM(1)=LINEM(2)
                  DEVM(2)=TDEVM
                  LINEM(2)=LLINEM
                END IF
              END IF
            END IF
          END IF
          SIGMA=SIGMA+DEVI(J)**2
        END DO
        SIGMA=SQRT(SIGMA/REAL(MM-1))
        WRITE(*,*)
        WRITE(*,'(A,F6.3)')'>>> Sigma: ',SIGMA
        WRITE(*,*)
        WRITE(*,101)'Largest deviations:'
        DO J=1,4
          WRITE(*,'(A,I3,5X,A,F12.4,5X,A,F12.4)')'Line #',LINEM(J),
     +     'Wavelength: ',FLAM(LINEM(J)),'Deviation: ',DEVM(J)
        END DO
C
        WRITE(*,*)
        WRITE(*,100)'Plots (y/n) '
        CRET(1:1)=READC('y','yn')
        IF(CRET.EQ.'n') GOTO 210
C
        YMAX=1.
        DO J=1,MLINES
          IF(ABS(DEVI2(J)).GT.YMAX) YMAX=ABS(DEVI2(J))
        END DO
        YMIN=-YMAX
        RANGE=YMAX-YMIN
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGSCH(1.5)
          CALL PGSUBP(1,2)
          CALL PGENV(0.,REAL(NCHAN),YMIN,YMAX,0,1)
          CALL PGLABEL('Channel','\\gD\\gl','Plot #1: '//ARCFILE)
          CALL PGPOINT(MLINES,XPOL,DEVI2,17)
          DO I=1,MLINES
            IF(ABS(DEVI(I)).GT.SIGMA) THEN
              J=LINEN(I)
              CALL PGNUMB(J,0,1,STRIN,NOCH)
              YTE=DEVI2(I)+RANGE/30.
              CALL PGPTEXT(XPOL(I),YTE,0.,0.5,STRIN(1:NOCH))
            END IF
          END DO
        END DO
C
        DO J=1,NCHAN
          POL=A(NTERMS)
          DO K=NTERMS-1,1,-1
            POL=POL*REAL(J)+A(K)
          END DO
            XXX(J)=REAL(J)
          YYY(J)=POL-B(1)-XXX(J)*B(2)
        END DO
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM))CALL PGSCI(2)
          CALL PGLINE(NCHAN,XXX,YYY)
          IF(LCOLOR(ITERM))CALL PGSCI(1)
        END DO
C
        YMAX=1.
        DO J=1,MLINES
          IF(ABS(DEVI(J)).GT.YMAX) YMAX=ABS(DEVI(J))
        END DO
        YMIN=-YMAX
        RANGE=YMAX-YMIN
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGENV(0.,REAL(NCHAN),YMIN,YMAX,0,1)
          CALL PGIDEN_RED
          CALL PGLABEL('Channel','\\gD\\gl','Plot #2: '//ARCFILE)
          CALL PGPOINT(MLINES,XPOL,DEVI,17)
          DO I=1,MLINES
            IF(ABS(DEVI(I)).GT.SIGMA) THEN
              J=LINEN(I)
              CALL PGNUMB(J,0,1,STRIN,NOCH)
              YTE=DEVI(I)+RANGE/30.
              CALL PGPTEXT(XPOL(I),YTE,0.,0.5,STRIN(1:NOCH))
            END IF
          END DO
        END DO
C
        WRITE(*,*)
        WRITE(*,101)'>>> PLOT 1: POINTS: differences between '//
     +   'data and straigh line fit'
        WRITE(*,101)'            LINE  : differences between '//
     +   'polynomial and straigh line'
        WRITE(*,*)
        WRITE(*,101)'>>> PLOT 2: POINTS: differences between '//
     +   'data and fit to a polynomial'
        WRITE(*,*)
210     CONTINUE
        WRITE(*,100)'Iterate (y/n) '
        CITER(1:1)=READC('y','yn')
        IF(CITER.EQ.'y')THEN
          GOTO 205
        END IF
C
        WRITE(*,100)'Polynomial (output) file name'
        POLFILE=OUTFILEX(15,'@',0,0,.0,.0,3,.FALSE.)
        DO K=1,NTERMS
          WRITE(15,*) K,A(K)
        END DO
        CLOSE(15)
C
        CALL PGEND
        STOP
100     FORMAT(A,$)
101     FORMAT(A)
110     FORMAT(A,I6)
        END
C
C******************************************************************************
C
        SUBROUTINE FLSORT(CMAX,CMIN)
        IMPLICIT NONE
        REAL CMAX,CMIN
        INCLUDE 'redlib.inc'
C
        INTEGER I
        INTEGER IMIN,IMAX
        REAL S(NCMAX)
        COMMON/BLKS/S
        COMMON/BLKI/IMIN,IMAX
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        CMAX=-1.E7
        DO 10 I=IMIN,IMAX
          IF(S(I).LT.CMAX) GO TO 10
          CMAX=S(I)
  10   CONTINUE
        CMIN=+1.E7
        DO 11 I=IMIN,IMAX
          IF(S(I).GT.CMIN) GO TO 11
          CMIN=S(I)
  11   CONTINUE
        RETURN
        END
C
C******************************************************************************
C
        SUBROUTINE PLOTCH
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
C
        INTEGER MAXLIN
        PARAMETER(MAXLIN=1000)
C
        INTEGER I,J
        INTEGER IMIN,IMAX
        INTEGER LINE(MAXLIN),NLINES
        INTEGER NUMB,NOCH
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        REAL S(NCMAX)
        REAL PAR(3,MAXLIN),XPOS(MAXLIN)
        REAL CMAX,CMIN
        REAL RANGE,CRAN
        REAL XMIN,XMAX,YMIN,YMAX
        REAL X(NCMAX),Y(NCMAX),YTE,XJ
        REAL XP(100),YP(100)
        REAL AMP,XBAR,SIGMA,Z1
        REAL XPO,YPO
        REAL XXPO(1),YYPO(1)
        CHARACTER*5 STRIN
        CHARACTER*20 TTER
        CHARACTER*75 ARCFILE
        LOGICAL LCOLOR(MAX_ID_RED)
C
        COMMON/BLKS/S
        COMMON/BLKNCHAN/NCHAN
        COMMON/BLKI/IMIN,IMAX
        COMMON/BLKP1/LINE,NLINES
        COMMON/BLKP2/PAR,XPOS
        COMMON/BLKARC/ARCFILE
        COMMON/BLKTER/TTER
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        CALL FLSORT(CMAX,CMIN)
C
C
        RANGE=CMAX-CMIN
        YMAX=CMAX+RANGE/50.
        YMIN=CMIN-RANGE/50.
C
        CRAN=REAL(IMAX)-REAL(IMIN)
        XMAX=REAL(IMAX)+CRAN/50.
        XMIN=REAL(IMIN)-CRAN/50.
        IF(XMIN.LT.0.) XMIN=0.
        IF(XMAX.GT.REAL(NCHAN+1)) XMAX=REAL(NCHAN+1)
C
        DO I=1,NCHAN
          X(I)=REAL(I)
          Y(I)=S(I)
        END DO
C
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGSUBP(1,1)
          CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,0)
          CALL PGIDEN_RED
          CALL PGLABEL('Channel','No. of counts',ARCFILE)
          CALL PGBIN(NCHAN,X,Y,.TRUE.)
        END DO
C
        DO I=1,NLINES
          IF((REAL(LINE(I)).GE.XMIN).AND.(REAL(LINE(I)).LE.XMAX))THEN
            NUMB=LINE(I)
            CALL PGNUMB(I,0,1,STRIN,NOCH)
            YTE=S(NUMB)+RANGE/30.
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              CALL PGPTEXT(REAL(NUMB),YTE,0.,0.5,STRIN(1:NOCH))
            END DO
          END IF
        END DO
C
        DO I=1,NLINES
          IF(LINE(I).GT.IMIN.AND.LINE(I).LT.IMAX) THEN
            AMP=PAR(1,I)
            XBAR=PAR(2,I)
            SIGMA=PAR(3,I)
            DO J=1,100
              XJ=(REAL(J)-50.)/50.*4.
              XP(J)=REAL(LINE(I))+XJ
              Z1=AMP*EXP(-(XJ-XBAR)**2/SIGMA/SIGMA)
              YP(J)=S(LINE(I))*Z1
            END DO
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              IF(LCOLOR(ITERM))CALL PGSCI(2)
              CALL PGLINE(100,XP,YP)
              IF(LCOLOR(ITERM))CALL PGSCI(1)
              XPO=XPOS(I)
              YPO=AMP*S(LINE(I))
              XXPO(1)=XPO
              YYPO(1)=YPO
              CALL PGPOINT(1,XXPO,YYPO,124)
            END DO
          END IF
        END DO
C
        RETURN
        END
C
C******************************************************************************
C
        REAL FUNCTION FUNKFITLIN(X)
        IMPLICIT NONE
        REAL X(3)
C
        INTEGER NPT,I
        REAL SX(9),SY(9)
        REAL SM,AMP,XBAR,SIGMA
        REAL X1,Z1
C
        COMMON/BLKF1/SX,SY
        COMMON/BLKF2/NPT
C------------------------------------------------------------------------------
        SM=0.
        AMP=X(1)
        XBAR=X(2)
        SIGMA=X(3)
C
        IF(X(1).LE.0.0) GO TO 10
        IF(X(3).LE.0.0) GO TO 10
        DO I=1,NPT
          X1=SX(I)
          Z1=AMP*EXP(-(X1-XBAR)**2/SIGMA/SIGMA)
          SM=SM+(SY(I)-Z1)**2
        END DO
        FUNKFITLIN=SM
        RETURN
10      FUNKFITLIN=1.E10
        RETURN
        END
