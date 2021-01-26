C------------------------------------------------------------------------------
C Version 30-May-1998                                        file: interpmove.f
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
C Program: interpmove
C Classification: arithmetic & manipulations
C Description: Interpolation/extrapolation of data in an image by using 
C polynomials around a previously fitted polynomial.
C
Comment
C
C Interpola en una imagen alrededor de una cierta estructura que es
C ajustada con la ayuda de un polinomio
C
        PROGRAM INTERPMOVE
        IMPLICIT NONE
C
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER READILIM
        REAL READF
C
        INTEGER NFITMAX
        PARAMETER (NFITMAX=1024)         !maximum number of points to be fitted
C
        INTEGER I,J,K,KK
        INTEGER NC1,NC2,NS1,NS2,NC0,NS0
        INTEGER NNC1,NNC2,NNS1,NNS2
        INTEGER NB
        INTEGER NPT
        INTEGER NDEG,NTERMS
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        REAL A(NCMAX,NSMAX),ERR(NCMAX,NSMAX)
        REAL ZMIN,ZMAX
        REAL XC,YC,XCDUM,YCDUM
        REAL XXC(1),YYC(1)
        REAL XMIN,XMAX,YMIN,YMAX
        REAL XPOL(NCMAX),YPOL(NCMAX)
        REAL XFIT(NFITMAX),YFIT(NFITMAX)
        REAL CHISQR,B(20)
        REAL DIST,DMIN
        REAL FIMIN,FIMAX
        REAL XV1,XV2,YV1,YV2
        CHARACTER*1 CERR,CH,CSURE,CFITOK,CCONT
        CHARACTER*1 CREMOVE,CTYPEPLOT,CSAVE
        CHARACTER*50 CDUMMY
        CHARACTER*75 INFILE,ERRFILE,OUTFILE,MARKFILE
        LOGICAL LEXIT,LBEXIST
        LOGICAL LCOLOR(MAX_ID_RED)
C
        COMMON/BLK0/NSCAN,NCHAN
        COMMON/BLK1/A
        COMMON/BLK2/NC1,NC2,NS1,NS2
        COMMON/BLK3/ZMIN,ZMAX
        COMMON/BLK4/INFILE
        COMMON/BLK5/CERR
        COMMON/BLK6/ERR
        COMMON/BLKLIMITS/XMIN,XMAX,YMIN,YMAX
        COMMON/BLKFIT1/NPT
        COMMON/BLKFIT2/XFIT,YFIT
        COMMON/BLKPOL/XPOL,YPOL
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        THISPROGRAM='interpmove'
        CALL WELCOME('18-Septiembre-1997')
C------------------------------------------------------------------------------
        NPT=0
        CTYPEPLOT='@'
C------------------------------------------------------------------------------
        CALL RPGBEGIN(NTERM,IDN,LCOLOR)
        CALL BUTTSYB(4)
        CALL BUTTQPR(XV1,XV2,YV1,YV2)
        CALL BUTTSPR(XV1,XV2,YV1,0.60)
        CALL BUTTQBR(XV1,XV2,YV1,YV2)
        CALL BUTTSBR(XV1,XV2,0.70,YV2)
C------------------------------------------------------------------------------
        WRITE(*,100)'Work with error images (y/n) '
        CERR(1:1)=READC('n','yn')
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
          ERRFILE=INFILEX(30,ERRFILE,NSCAN,NCHAN,STWV,DISP,21,.TRUE.) !match
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
        CALL BUTTON(7,'[s]et BG/FG',0)
        CALL BUTTON(8,'Min[,]Max',0)
        CALL BUTTON(13,'[p]lot',0)
        CALL BUTTON(14,'[i]nterp.',0)
        CALL BUTTON(19,'[q]uit&save',0)
        CALL BUTTON(20,'[r]emove',0)
C------------------------------------------------------------------------------
10      IF(NPT.GT.0)THEN
          CALL BUTTON(20,'[r]emove',0)
          WRITE(*,101)'* Remeber: key [2] output marks file'
        ELSE
          CALL BUTTON(20,'[r]emove',3)
          WRITE(*,101)'* Remeber: key [1] input marks file'
        END IF
        IF(NPT.GT.1)THEN
          CALL BUTTON(14,'[i]nterp.',0)
        ELSE
          CALL BUTTON(14,'[i]nterp.',3)
        END IF
        CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
        IF(CH.EQ.'z')THEN
          NB=1
        ELSEIF(CH.EQ.'w')THEN
          NB=2
        ELSEIF(CH.EQ.'1')THEN
          NB=3
        ELSEIF(CH.EQ.'s')THEN
          NB=7
        ELSEIF(CH.EQ.',')THEN
          NB=8
        ELSEIF(CH.EQ.'2')THEN
          NB=9
        ELSEIF(CH.EQ.'p')THEN
          NB=13
        ELSEIF(CH.EQ.'i')THEN
          CALL BUTTQEX(14,LBEXIST)
          IF(LBEXIST) NB=14
        ELSEIF(CH.EQ.'q')THEN
          NB=19
        ELSEIF(CH.EQ.'r')THEN
          CALL BUTTQEX(20,LBEXIST)
          IF(LBEXIST) NB=20
        ELSE
          CALL IFBUTTON(XC,YC,NB)
        END IF
C------------------------------------------------------------------------------
        IF(NB.EQ.0)THEN
          IF(CH.EQ.'A')THEN
            IF((XC.GE.XMIN).AND.(XC.LE.XMAX).AND.
     +         (YC.GE.YMIN).AND.(YC.LE.YMAX)) THEN
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                IF(LCOLOR(ITERM)) CALL PGSCI(2)
                CALL PGSCH(1.5)
                XXC(1)=XC
                YYC(1)=YC
                CALL PGPOINT(1,XXC,YYC,5)
                CALL PGSCH(1.0)
              END DO
              WRITE(*,100)'Press mouse button again to confirm '//
     +         'or abort...'
              XCDUM=XC
              YCDUM=YC
              CALL RPGBAND(0,0,0.,0.,XCDUM,YCDUM,CH)
              WRITE(*,101)' thanks!'
              IF(CH.EQ.'A')THEN
                NPT=NPT+1
                XFIT(NPT)=XC
                YFIT(NPT)=YC
                DO ITERM=NTERM,1,-1
                  CALL PGSLCT(IDN(ITERM))
                  IF(LCOLOR(ITERM)) CALL PGSCI(3)
                  CALL PGSCH(1.5)
                  XXC(1)=XC
                  YYC(1)=YC
                  CALL PGPOINT(1,XXC,YYC,5)
                  CALL PGSCH(1.0)
                  IF(LCOLOR(ITERM)) CALL PGSCI(1)
               END DO
              ELSE
                DO ITERM=NTERM,1,-1
                  CALL PGSLCT(IDN(ITERM))
                  CALL PGSCI(0)
                  CALL PGSCH(1.5)
                  XXC(1)=XC
                  YYC(1)=YC
                  CALL PGPOINT(1,XXC,YYC,5)
                  CALL PGSCH(1.0)
                  CALL PGSCI(1)
                END DO
              END IF
            END IF
          END IF
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
          WRITE(*,100)'Input marks file name'
          MARKFILE=INFILEX(10,'@',0,0,0.,0.,3,.FALSE.)
          NPT=0
20        READ(10,*,END=22)XFIT(NPT+1),YFIT(NPT+1)
          NPT=NPT+1
          GOTO 20
22        CLOSE(10)
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(3)
            CALL PGSCH(1.5)
            CALL PGPOINT(NPT,XFIT,YFIT,5)
            CALL PGSCH(1.0)
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END DO
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
        ELSEIF(NB.EQ.9)THEN
          WRITE(*,100)'Output marks file name'
          MARKFILE=OUTFILEX(10,'@',0,0,0.,0.,3,.FALSE.)
          DO K=1,NPT
            WRITE(10,*)XFIT(K),YFIT(K)
          END DO
          CLOSE(10)
C..............................................................................
        ELSEIF(NB.EQ.13)THEN
          CALL BUTTON(13,'[p]lot',5)
          WRITE(*,100)'[x]-plot or [y]-plot (x/y) '
          CTYPEPLOT(1:1)=READC(CTYPEPLOT,'xy')
          IF(CTYPEPLOT.EQ.'x')THEN
            CALL PLOT_X
          ELSE
            CALL PLOT_Y
          END IF
          XMIN=REAL(NC1)-0.5
          XMAX=REAL(NC2)+0.5
          YMIN=REAL(NS1)-0.5
          YMAX=REAL(NS2)+0.5
          CALL BUTTQPR(XV1,XV2,YV1,YV2)
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGVPORT(XV1,XV2,YV1,YV2)
            CALL PGWINDOW(XMIN,XMAX,YMIN,YMAX)
          END DO
          CALL BUTTON(13,'[p]lot',0)
C..............................................................................
        ELSEIF(NB.EQ.14)THEN
          CALL BUTTON(14,'[i]nterp.',5)
          LEXIT=.FALSE.
          DO WHILE(.NOT.LEXIT)
            WRITE(*,100)'Polynomial degree '
            IF(NPT.LT.20)THEN
              NDEG=READILIM('@',0,NPT-1)
            ELSE
              NDEG=READILIM('@',0,19)
            END IF
            NTERMS=NDEG+1
            CALL POLFIT(XFIT,YFIT,YFIT,NPT,NTERMS,0,B,CHISQR)
            DO K=1,NCHAN
              XPOL(K)=REAL(K)
              YPOL(K)=B(NTERMS)
              DO KK=NTERMS-1,1,-1
                YPOL(K)=YPOL(K)*XPOL(K)+B(KK)
              END DO
            END DO
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              IF(LCOLOR(ITERM)) CALL PGSCI(5)
              CALL PGSLW(3)
              CALL PGLINE(NCHAN,XPOL,YPOL)
              CALL PGSLW(1)
              IF(LCOLOR(ITERM)) CALL PGSCI(1)
            END DO
            WRITE(*,100)'Is the fit OK (y/n) '
            CFITOK(1:1)=READC('y','yn')
            IF(CFITOK.EQ.'y')THEN
              LEXIT=.TRUE.
            ELSE
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                CALL PGSCI(0)
                CALL PGSLW(3)
                CALL PGLINE(NCHAN,XPOL,YPOL)
                CALL PGSLW(1)
                CALL PGSCI(1)
              END DO
            END IF
          END DO
          WRITE(*,100)'Continue with y-interpolation (y/n) '
          CCONT(1:1)=READC('y','yn')
          IF(CCONT.EQ.'y')THEN
            FIMIN=YPOL(1)
            FIMAX=FIMIN
            DO J=2,NCHAN
              IF(YPOL(J).LT.FIMIN) FIMIN=YPOL(J)
              IF(YPOL(J).GT.FIMAX) FIMAX=YPOL(J)
            END DO
            WRITE(*,100)'>>> Minimum scan: '
            WRITE(*,*) FIMIN
            WRITE(*,100)'>>> Maximum scan: '
            WRITE(*,*) FIMAX
            CALL INTERP_Y
          END IF
          CALL BUTTON(14,'[i]nterp.',0)
C..............................................................................
        ELSEIF(NB.EQ.19)THEN
          CALL BUTTON(19,'[q]uit&save',5)
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
          CALL BUTTON(19,'[q]uit&save',0)
C..............................................................................
        ELSEIF(NB.EQ.20)THEN
          CALL BUTTON(20,'[r]emove',5)
          WRITE(*,100)'Remove [a]ll or [s]elected marks (a/s) '
          CREMOVE(1:1)=READC('s','as')
          IF(CREMOVE.EQ.'a')THEN
            WRITE(*,100)'Are you sure (y/n) '
            CSURE(1:1)=READC('n','yn')
            IF(CSURE.EQ.'y')THEN
              NPT=0
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                IF(ITERM.EQ.1)THEN
                  CALL RPGERASW(0.,1.,0.,0.65)
                ELSE
                  CALL PGPAGE
                END IF
              END DO
              CALL PLOT_A
            END IF
          ELSE
            LEXIT=.FALSE.
            DO WHILE(.NOT.LEXIT)
              WRITE(*,100)'Press mouse button to remove mark or '//
     +         'abort...'
              CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
              WRITE(*,101)' thanks!'
              IF(CH.EQ.'A')THEN
                IF(NPT.EQ.1)THEN
                  KK=1
                ELSE
                  DMIN=(XC-XFIT(1))*(XC-XFIT(1))+
     +             (YC-YFIT(1))*(YC-YFIT(1))
                  KK=1
                  DO K=2,NPT
                    DIST=(XC-XFIT(K))*(XC-XFIT(K))+
     +               (YC-YFIT(K))*(YC-YFIT(K))
                    IF(DIST.LT.DMIN)THEN
                      DMIN=DIST
                      KK=K
                    END IF
                  END DO
                END IF
                DO ITERM=NTERM,1,-1
                  CALL PGSLCT(IDN(ITERM))
                  CALL PGSCI(0)
                  CALL PGSCH(1.5)
                  CALL PGPOINT(1,XFIT(KK),YFIT(KK),5)
                  CALL PGSCH(1.0)
                  CALL PGSCI(1)
                END DO
                IF(KK.LT.NPT)THEN
                  DO K=KK,NPT-1
                    XFIT(K)=XFIT(K+1)
                    YFIT(K)=YFIT(K+1)
                  END DO
                END IF
                NPT=NPT-1
                IF(NPT.EQ.0)THEN
                  WRITE(*,101)'WARNING: all the marks have been '//
     +             'removed.'
                  LEXIT=.TRUE.
                END IF
              ELSE
                LEXIT=.TRUE.
              END IF
            END DO
          END IF
          CALL BUTTON(20,'[r]emove',0)
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
        INTEGER NFITMAX
        PARAMETER (NFITMAX=1024)         !maximum number of points to be fitted
C
        INTEGER NPT
        INTEGER NC1,NC2,NS1,NS2
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        REAL A(NCMAX,NSMAX)
        REAL TR(6),ZMIN,ZMAX
        REAL XMIN,XMAX,YMIN,YMAX
        REAL XFIT(NFITMAX),YFIT(NFITMAX)
        CHARACTER*75 INFILE
        LOGICAL LCOLOR(MAX_ID_RED)
C
        COMMON/BLK1/A
        COMMON/BLK2/NC1,NC2,NS1,NS2
        COMMON/BLK3/ZMIN,ZMAX
        COMMON/BLK4/INFILE
        COMMON/BLKFIT1/NPT
        COMMON/BLKFIT2/XFIT,YFIT
        COMMON/BLKLIMITS/XMIN,XMAX,YMIN,YMAX
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        TR(1)=0.
        TR(2)=1.
        TR(3)=0.
        TR(4)=0.
        TR(5)=0.
        TR(6)=1.
C
        XMIN=REAL(NC1)-0.5
        XMAX=REAL(NC2)+0.5
        YMIN=REAL(NS1)-0.5
        YMAX=REAL(NS2)+0.5
C
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(ITERM.EQ.1)THEN
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
          IF(LCOLOR(ITERM)) CALL PGSCI(3)
          CALL PGSCH(1.5)
          CALL PGPOINT(NPT,XFIT,YFIT,5)
          CALL PGSCH(1.0)
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
        END DO
C
        END
C
C******************************************************************************
C
C Dibuja un corte promedio en X
        SUBROUTINE PLOT_X
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'futils.inc'
C
        INTEGER NFITMAX
        PARAMETER (NFITMAX=1024)         !maximum number of points to be fitted
C
        INTEGER NPT
        INTEGER I,J
        INTEGER J0,J1,J2
        INTEGER I1,I2,NI
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        REAL MEANI
        REAL A(NCMAX,NSMAX)
        REAL X(NCMAX),Y(NCMAX)
        REAL XMINS,XMAXS,YMINS,YMAXS,DX,DY
        REAL XMIN,XMAX,YMIN,YMAX
        REAL XC,YC
        REAL XFIT(NFITMAX),YFIT(NFITMAX)
        REAL XV1,XV2,YV1,YV2
        REAL XX0(1),YY0(1)
        CHARACTER*1 COPC,CH,CMARKOK
        LOGICAL IFSCAN(NSMAX)
        LOGICAL LEXIT
        LOGICAL LCOLOR(MAX_ID_RED)
C
        COMMON/BLK0/NSCAN,NCHAN
        COMMON/BLK1/A
        COMMON/BLKLIMITS/XMIN,XMAX,YMIN,YMAX
        COMMON/BLKFIT1/NPT
        COMMON/BLKFIT2/XFIT,YFIT
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        WRITE(*,101)'* Plotting X-cut: '
        DO I=1,NSCAN
          IFSCAN(I)=.FALSE.
        END DO
        LEXIT=.FALSE.
        DO WHILE(.NOT.LEXIT)
          WRITE(*,100)'Scan region (0,0=EXIT) '
          CALL READ2I('0,0',I1,I2)
          IF((I1.EQ.0).AND.(I2.EQ.0))THEN
            NI=0
            DO I=1,NSCAN
              IF(IFSCAN(I)) NI=NI+1
            END DO
            IF(NI.EQ.0)THEN
              WRITE(*,101)'ERROR: no. of scans = 0. Try again.'
            ELSE
              LEXIT=.TRUE.
            END IF
          ELSEIF((I1.LT.1).OR.(I2.GT.NSCAN).OR.(I1.GT.I2))THEN
            WRITE(*,101)'ERROR: out of range. Try again.'
          ELSE
            DO I=I1,I2
              IFSCAN(I)=.TRUE.
            END DO
          END IF
        END DO
C------------------------------------------------------------------------------
        DO J=1,NCHAN
          X(J)=REAL(J)
          Y(J)=0.
          DO I=1,NSCAN
            IF(IFSCAN(I)) Y(J)=Y(J)+A(J,I)
          END DO
          Y(J)=Y(J)/REAL(NI)
        END DO
C
        J1=1
        J2=NCHAN
C
10      XMINS=REAL(J1)
        XMAXS=REAL(J2)
        DX=XMAXS-XMINS
        XMINS=XMINS-DX/20.
        XMAXS=XMAXS+DX/20.
C
        CALL FINDMML(NCHAN,J1,J2,Y,YMINS,YMAXS)
        DY=YMAXS-YMINS
        YMINS=YMINS-DY/20.
        YMAXS=YMAXS+DY/20.
C
        CALL BUTTQBR(XV1,XV2,YV1,YV2)
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGVPORT(0.40,XV2,YV1,YV2)
ccc          CALL PGVPORT(0.40,X2VPORT,Y3VPORT,Y4VPORT)
          CALL PGWINDOW(XMINS,XMAXS,YMINS,YMAXS)
          CALL PGSCI(0)
          CALL PGRECT(XMINS,XMAXS,YMINS,YMAXS)
          CALL PGVPORT(0.43,XV2,YV1+0.05,YV2)
ccc          CALL PGVPORT(0.43,X2VPORT,Y3VPORT+0.05,Y4VPORT)
          CALL PGSCI(1)
          CALL PGSCH(0.7)
          CALL PGSLW(3)
          CALL PGBOX('BCNTS',0.0,0,'BCNTS',0.0,0)
          CALL PGSCH(1.0)
          CALL PGMTEXT('B',2.5,0.5,0.5,'channel')
          CALL PGMTEXT('L',2.0,0.5,0.5,'no. counts')
          IF(LCOLOR(ITERM)) CALL PGSCI(2)
          CALL PGBIN(NCHAN,X,Y,.TRUE.)
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
          CALL PGSLW(1)
        END DO
C------------------------------------------------------------------------------
ccc20      WRITE(*,100)'[z]oom, [w]hole, e[x]it, [m]ark (z/w/x/m) '
        WRITE(*,100)'[z]oom, [w]hole, e[x]it, [m]ark (z/w/x/m) '
        COPC(1:1)=READC('x','zwxm')
C
        IF(COPC.EQ.'x')THEN
          RETURN
        ELSEIF(COPC.EQ.'z')THEN
          WRITE(*,100)'Press mouse button...'
          IF(LCOLOR(1)) CALL PGSCI(5)
          CALL RPGBAND(6,0,0.,0.,XC,YC,CH)
          IF(LCOLOR(1)) CALL PGSCI(1)
          J1=NINT(XC)
          IF(J1.LT.1) J1=1
          IF(J1.GT.NCHAN) J1=NCHAN
          WRITE(*,100)'again...'
          IF(LCOLOR(1)) CALL PGSCI(5)
          CALL RPGBAND(4,0,REAL(J1),0.,XC,YC,CH)
          IF(LCOLOR(1)) CALL PGSCI(1)
          J2=NINT(XC)
          IF(J2.LT.1) J2=1
          IF(J2.GT.NCHAN) J2=NCHAN
          WRITE(*,101)'OK!'
          IF(J1.GT.J2)THEN
            J0=J1
            J1=J2
            J2=J0
          END IF
          GOTO 10
        ELSEIF(COPC.EQ.'w')THEN
          J1=1
          J2=NCHAN
          GOTO 10
        ELSEIF(COPC.EQ.'m')THEN
          WRITE(*,100)'Press mouse button...'
          CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
          WRITE(*,101)' thanks!'
          MEANI=0
          DO I=1,NSCAN
            IF(IFSCAN(I)) MEANI=MEANI+REAL(I)
          END DO
          MEANI=MEANI/REAL(NI)
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL BUTTQPR(XV1,XV2,YV1,YV2)
            CALL PGVPORT(XV1,XV2,YV1,YV2)
            CALL PGWINDOW(XMIN,XMAX,YMIN,YMAX)
            IF(LCOLOR(ITERM)) CALL PGSCI(2)
            CALL PGSCH(1.5)
            XX0(1)=XC
            YY0(1)=MEANI
            CALL PGPOINT(1,XX0,YY0,5)
            CALL PGSCH(1.0)
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END DO
          WRITE(*,100)'Is the mark location ok (y/n) '
          CMARKOK(1:1)=READC('y','yn')
          IF(CMARKOK.EQ.'n')THEN
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              CALL PGSCI(0)
              CALL PGSCH(1.5)
              XX0(1)=XC
              YY0(1)=MEANI
              CALL PGPOINT(1,XX0,YY0,5)
              CALL PGSCH(1.0)
              CALL PGSCI(1)
            END DO
          ELSE
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              IF(LCOLOR(ITERM)) CALL PGSCI(3)
              CALL PGSCH(1.5)
              XX0(1)=XC
              YY0(1)=MEANI
              CALL PGPOINT(1,XX0,YY0,5)
              CALL PGSCH(1.0)
              IF(LCOLOR(ITERM)) CALL PGSCI(1)
            END DO
            NPT=NPT+1
            XFIT(NPT)=XC
            YFIT(NPT)=MEANI
          END IF
          GOTO 10
        END IF
C------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
        END
C
C******************************************************************************
C
C Dibuja un corte promedio en Y
        SUBROUTINE PLOT_Y
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'futils.inc'
C
        INTEGER NFITMAX
        PARAMETER (NFITMAX=1024)         !maximum number of points to be fitted
C
        INTEGER NPT
        INTEGER I,J
        INTEGER I0,I1,I2
        INTEGER J1,J2,NJ
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        REAL MEANJ
        REAL A(NCMAX,NSMAX)
        REAL X(NCMAX),Y(NCMAX)
        REAL XMINS,XMAXS,YMINS,YMAXS,DX,DY
        REAL XMIN,XMAX,YMIN,YMAX
        REAL XC,YC
        REAL XFIT(NFITMAX),YFIT(NFITMAX)
        REAL XV1,XV2,YV1,YV2
        REAL XX0(1),YY0(1)
        CHARACTER*1 COPC,CH,CMARKOK
        LOGICAL IFCHAN(NCMAX)
        LOGICAL LEXIT
        LOGICAL LCOLOR(MAX_ID_RED)
C
        COMMON/BLK0/NSCAN,NCHAN
        COMMON/BLK1/A
        COMMON/BLKLIMITS/XMIN,XMAX,YMIN,YMAX
        COMMON/BLKFIT1/NPT
        COMMON/BLKFIT2/XFIT,YFIT
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        WRITE(*,101)'* Plotting Y-cut: '
        DO J=1,NCHAN
          IFCHAN(J)=.FALSE.
        END DO
        LEXIT=.FALSE.
        DO WHILE(.NOT.LEXIT)
          WRITE(*,100)'Channel region (0,0=EXIT) '
          CALL READ2I('0,0',J1,J2)
          IF((J1.EQ.0).AND.(J2.EQ.0))THEN
            NJ=0
            DO J=1,NCHAN
              IF(IFCHAN(J)) NJ=NJ+1
            END DO
            IF(NJ.EQ.0)THEN
              WRITE(*,101)'ERROR: no. of channels = 0. Try again.'
            ELSE
              LEXIT=.TRUE.
            END IF
          ELSEIF((J1.LT.1).OR.(J2.GT.NCHAN).OR.(J1.GT.J2))THEN
            WRITE(*,101)'ERROR: out of range. Try again.'
          ELSE
            DO J=J1,J2
              IFCHAN(J)=.TRUE.
            END DO
          END IF
        END DO
C------------------------------------------------------------------------------
        DO I=1,NSCAN
          X(I)=REAL(I)
          Y(I)=0.
          DO J=1,NCHAN
            IF(IFCHAN(J)) Y(I)=Y(I)+A(J,I)
          END DO
          Y(I)=Y(I)/REAL(NJ)
        END DO
C
        I1=1
        I2=NSCAN
C
10      XMINS=REAL(I1)
        XMAXS=REAL(I2)
        DX=XMAXS-XMINS
        XMINS=XMINS-DX/20.
        XMAXS=XMAXS+DX/20.
C
        CALL FINDMML(NSCAN,I1,I2,Y,YMINS,YMAXS)
        DY=YMAXS-YMINS
        YMINS=YMINS-DY/20.
        YMAXS=YMAXS+DY/20.
C
        CALL BUTTQBR(XV1,XV2,YV1,YV2)
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGVPORT(0.40,XV2,YV1,YV2)
ccc          CALL PGVPORT(0.40,X2VPORT,Y3VPORT,Y4VPORT)
          CALL PGWINDOW(XMINS,XMAXS,YMINS,YMAXS)
          CALL PGSCI(0)
          CALL PGRECT(XMINS,XMAXS,YMINS,YMAXS)
          CALL PGVPORT(0.43,XV2,YV1+0.05,YV2)
ccc          CALL PGVPORT(0.43,X2VPORT,Y3VPORT+0.05,Y4VPORT)
          CALL PGSCI(1)
          CALL PGSCH(0.7)
          CALL PGSLW(3)
          CALL PGBOX('BCNTS',0.0,0,'BCNTS',0.0,0)
          CALL PGSCH(1.0)
          CALL PGMTEXT('B',2.5,0.5,0.5,'scan')
          CALL PGMTEXT('L',2.0,0.5,0.5,'no. counts')
          IF(LCOLOR(ITERM)) CALL PGSCI(2)
          CALL PGBIN(NSCAN,X,Y,.TRUE.)
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
          CALL PGSLW(1)
        END DO
C------------------------------------------------------------------------------
ccc20      WRITE(*,100)'[z]oom, [w]hole, e[x]it, [m]ark (z/w/x/m) '
        WRITE(*,100)'[z]oom, [w]hole, e[x]it, [m]ark (z/w/x/m) '
        COPC(1:1)=READC('x','zwxm')
C
        IF(COPC.EQ.'x')THEN
          RETURN
        ELSEIF(COPC.EQ.'z')THEN
          WRITE(*,100)'Press mouse button...'
          IF(LCOLOR(1)) CALL PGSCI(5)
          CALL RPGBAND(6,0,0.,0.,XC,YC,CH)
          IF(LCOLOR(1)) CALL PGSCI(1)
          I1=NINT(XC)
          IF(I1.LT.1) I1=1
          IF(I1.GT.NSCAN) I1=NSCAN
          WRITE(*,100)'again...'
          IF(LCOLOR(1)) CALL PGSCI(5)
          CALL RPGBAND(4,0,REAL(I1),0.,XC,YC,CH)
          IF(LCOLOR(1)) CALL PGSCI(1)
          I2=NINT(XC)
          IF(I2.LT.1) I2=1
          IF(I2.GT.NSCAN) I2=NSCAN
          WRITE(*,101)'OK!'
          IF(I1.GT.I2)THEN
            I0=I1
            I1=I2
            I2=I0
          END IF
          GOTO 10
        ELSEIF(COPC.EQ.'w')THEN
          I1=1
          I2=NSCAN
          GOTO 10
        ELSEIF(COPC.EQ.'m')THEN
          WRITE(*,100)'Press mouse button...'
          CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
          WRITE(*,101)' thanks!'
          MEANJ=0
          DO J=1,NCHAN
            IF(IFCHAN(J)) MEANJ=MEANJ+REAL(J)
          END DO
          MEANJ=MEANJ/REAL(NJ)
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL BUTTQPR(XV1,XV2,YV1,YV2)
            CALL PGVPORT(XV1,XV2,YV1,YV2)
            CALL PGWINDOW(XMIN,XMAX,YMIN,YMAX)
            IF(LCOLOR(ITERM)) CALL PGSCI(2)
            CALL PGSCH(1.5)
            XX0(1)=MEANJ
            YY0(1)=XC
            CALL PGPOINT(1,XX0,YY0,5)
            CALL PGSCH(1.0)
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END DO
          WRITE(*,100)'Is the mark location ok (y/n) '
          CMARKOK(1:1)=READC('y','yn')
          IF(CMARKOK.EQ.'n')THEN
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              CALL PGSCI(0)
              CALL PGSCH(1.5)
              XX0(1)=MEANJ
              YY0(1)=XC
              CALL PGPOINT(1,XX0,YY0,5)
              CALL PGSCH(1.0)
              CALL PGSCI(1)
            END DO
          ELSE
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              IF(LCOLOR(ITERM)) CALL PGSCI(3)
              CALL PGSCH(1.5)
              XX0(1)=MEANJ
              YY0(1)=XC
              CALL PGPOINT(1,XX0,YY0,5)
              CALL PGSCH(1.0)
              IF(LCOLOR(ITERM)) CALL PGSCI(1)
            END DO
            NPT=NPT+1
            XFIT(NPT)=MEANJ
            YFIT(NPT)=XC
          END IF
          GOTO 10
        END IF
C------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
        END
C
C******************************************************************************
C
C Realiza la interpolacion en la direccion Y
        SUBROUTINE INTERP_Y
        IMPLICIT NONE
C
        INCLUDE 'redlib.inc'
        INCLUDE 'futils.inc'
        INTEGER READI
        INTEGER READILIM
C
        INTEGER I,J,K
        INTEGER NS_C1,NS_C2                                       !Central band
        INTEGER NS_A1,NS_A2                                         !Lower band
        INTEGER NS_B1,NS_B2                                         !Upper band
        INTEGER NCEN,NBANDL,NSKIPL,NBANDU,NSKIPU
        INTEGER NF
        INTEGER NDEG,NTERMS
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        REAL A(NCMAX,NSMAX)
        REAL ERR(NCMAX,NSMAX)
        REAL AA(NCMAX,NSMAX)
        REAL XMIN,XMAX,YMIN,YMAX
        REAL XPOL(NCMAX),YPOL(NCMAX)
        REAL YLIM1(NCMAX),YLIM2(NCMAX)
        REAL XF(NSMAX),YF(NSMAX)
        REAL CHISQR,B(20),XX,YY
        CHARACTER*1 CERR,CBOK,CFOK,CCONT
        LOGICAL LERROR
        LOGICAL LCOLOR(MAX_ID_RED)
        LOGICAL LOUT
C
        COMMON/BLK0/NSCAN,NCHAN
        COMMON/BLK1/A
        COMMON/BLK5/CERR
        COMMON/BLK6/ERR
        COMMON/BLKLIMITS/XMIN,XMAX,YMIN,YMAX
        COMMON/BLKPOL/XPOL,YPOL
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
10      WRITE(*,101)'* Define interpolation bands (in scans):'
        WRITE(*,100)'Central band width (odd please)'
        NCEN=READI('@')
        WRITE(*,100)'Lower band width'
        NBANDL=READI('@')
        WRITE(*,100)'Scans to skip between Central band and Lower band '
        NSKIPL=READI('0')
        WRITE(*,100)'Upper band width'
        NBANDU=READI('@')
        WRITE(*,100)'Scans to skip between Central band and Upper band '
        NSKIPU=READI('0')
C------------------------------------------------------------------------------
        DO ITERM=NTERM,1,-1
ccc          IF(ITERM.EQ.1)THEN
ccc          CALL RPGENV(XMIN,XMAX,YMIN,YMAX,0,-2)
ccc          ELSE
ccc          CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,-2)
ccc          END IF
          CALL PGWINDOW(XMIN,XMAX,YMIN,YMAX)
        END DO
C------------------------------------------------------------------------------
        DO J=1,NCHAN
          NS_C1=NINT(YPOL(J))-NCEN/2
          NS_C2=NINT(YPOL(J))+NCEN/2
          YLIM1(J)=NS_C1
          YLIM2(J)=NS_C2
        END DO
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM)) CALL PGSCI(2)
          CALL PGLINE(NCHAN,XPOL,YLIM1)
          CALL PGLINE(NCHAN,XPOL,YLIM2)
        END DO
        DO J=1,NCHAN
          NS_C1=NINT(YPOL(J))-NCEN/2
          NS_A2=NS_C1-NSKIPL-1
          NS_A1=NS_A2-NBANDL+1
          YLIM1(J)=NS_A1
          YLIM2(J)=NS_A2
        END DO
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM)) CALL PGSCI(3)
          CALL PGSLS(4)
          CALL PGLINE(NCHAN,XPOL,YLIM1)
          CALL PGLINE(NCHAN,XPOL,YLIM2)
          CALL PGSLS(1)
        END DO
        DO J=1,NCHAN
          NS_C2=NINT(YPOL(J))+NCEN/2
          NS_B1=NS_C2+NSKIPU+1
          NS_B2=NS_B1+NBANDU-1
          YLIM1(J)=NS_B1
          YLIM2(J)=NS_B2
        END DO
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM)) CALL PGSCI(4)
          CALL PGSLS(4)
          CALL PGLINE(NCHAN,XPOL,YLIM1)
          CALL PGLINE(NCHAN,XPOL,YLIM2)
          CALL PGSLS(1)
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
        END DO
C------------------------------------------------------------------------------
C Testeamos que no hay problemas con los limites
        LERROR=.FALSE.
        DO J=1,NCHAN
          NS_C1=NINT(YPOL(J))-NCEN/2
          NS_C2=NINT(YPOL(J))+NCEN/2
          NS_A2=NS_C1-NSKIPL-1
          NS_A1=NS_A2-NBANDL+1
          NS_B1=NS_C2+NSKIPU+1
          NS_B2=NS_B1+NBANDU-1
          IF(NS_B1.GE.NSCAN) LERROR=.TRUE.
          IF(NS_A2.LE.1) LERROR=.TRUE.
        END DO
        IF(LERROR)THEN
          WRITE(*,101)'WARNING: band limits out of frame.'
          WRITE(*,100)'Do you want to continue anyway (y/n) '
          CCONT(1:1)=READC('y','yn')
          IF(CCONT.EQ.'n') RETURN
        END IF
C------------------------------------------------------------------------------
        WRITE(*,100)'Are the band limits right (y/n) '
        CBOK(1:1)=READC('y','yn')
        IF(CBOK.EQ.'n')THEN
          CALL PLOT_A
          GOTO 10
        END IF
C
        DO J=1,NCHAN
          DO I=1,NSCAN
            AA(J,I)=A(J,I)
          END DO
        END DO
C
        WRITE(*,100)'Polynomial degree'
        NDEG=READILIM('@',0,19)
        NTERMS=NDEG+1
C------------------------------------------------------------------------------
        LOUT=.FALSE.
        WRITE(*,100)'Thinking...'
        DO J=1,NCHAN
          NS_C1=NINT(YPOL(J))-NCEN/2
          NS_C2=NINT(YPOL(J))+NCEN/2
          NS_A2=NS_C1-NSKIPL-1
          NS_A1=NS_A2-NBANDL+1
          NS_B1=NS_C2+NSKIPU+1
          NS_B2=NS_B1+NBANDU-1
          NF=0
          DO I=NS_A1,NS_A2
            IF((I.GE.1).AND.(I.LE.NSCAN))THEN
              NF=NF+1
              XF(NF)=REAL(I)
              YF(NF)=A(J,I)
            END IF
          END DO
          DO I=NS_B1,NS_B2
            IF((I.GE.1).AND.(I.LE.NSCAN))THEN
              NF=NF+1
              XF(NF)=REAL(I)
              YF(NF)=A(J,I)
            END IF
          END DO
          IF(NF.LT.NTERMS)THEN
            IF(.NOT.LOUT)THEN
              WRITE(*,'(A,I6,A)')'Channel #',J,'   --> Insuficient '//
     +         'no. of points to fit.'
              LOUT=.TRUE.
            END IF
          ELSE
            CALL POLFIT(XF,YF,YF,NF,NTERMS,0,B,CHISQR)
            DO I=NS_C1,NS_C2
              IF((I.GE.1).AND.(I.LE.NSCAN))THEN
                XX=REAL(I)
                YY=B(NTERMS)
                DO K=NTERMS-1,1,-1
                  YY=YY*XX+B(K)
                END DO
                A(J,I)=YY
              END IF
            END DO
          END IF
        END DO
        WRITE(*,101)'  ..OK!'
C
        CALL PLOT_A
C------------------------------------------------------------------------------
        WRITE(*,100)'Is the fit right (y/n) '
        CFOK(1:1)=READC('y','yn')
        IF(CFOK.EQ.'n')THEN
          DO J=1,NCHAN
            DO I=1,NSCAN
              A(J,I)=AA(J,I)
            END DO
          END DO
          RETURN
        END IF
C------------------------------------------------------------------------------
        IF(CERR.EQ.'y')THEN
          WRITE(*,100)'Thinking (fitting errors)...'
          DO J=1,NCHAN
            NS_C1=NINT(YPOL(J))-NCEN/2
            NS_C2=NINT(YPOL(J))+NCEN/2
            NS_A2=NS_C1-NSKIPL-1
            NS_A1=NS_A2-NBANDL+1
            NS_B1=NS_C2+NSKIPU+1
            NS_B2=NS_B1+NBANDU-1
            IF(NS_A1.LT.1) NS_A1=1
            IF(NS_B2.GT.NSCAN) NS_B2=NSCAN
            NF=0
            DO I=NS_A1,NS_A2
              NF=NF+1
              XF(NF)=REAL(I)
              YF(NF)=ERR(J,I)
            END DO
            DO I=NS_B1,NS_B2
              NF=NF+1
              XF(NF)=REAL(I)
              YF(NF)=ERR(J,I)
            END DO
            IF(NF.LT.NTERMS)THEN
              WRITE(*,'(A,I6,A)')'Channel #',J,'   --> Insuficient '//
     +         'no. of points to fit.'
            ELSE
              CALL POLFIT(XF,YF,YF,NF,NTERMS,0,B,CHISQR)
              DO I=NS_C1,NS_C2
                XX=REAL(I)
                YY=B(NTERMS)
                DO K=NTERMS-1,1,-1
                  YY=YY*XX+B(K)
                END DO
                ERR(J,I)=YY
              END DO
            END IF
          END DO
          WRITE(*,101)'  ..OK!'
        END IF
C------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
        END
