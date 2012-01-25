C Version 10-March-2005                                         file: sdistor.f
C @ ncl & fjg
C------------------------------------------------------------------------------
Comment
C
C Program: sdistor
C Classification: distortion
C Description: Corrects S-distortion.
C
Comment
C
C Calcula y corrige la distorsion-S en una imagen. Realiza tambien el trabajo
C con los errores si se solicita.
C En esta nueva version, el programa realiza la correcion utilizando un
C sistema que mejora el efecto de pixelizacion de la imagen. Para ello, despues
C de una primera correcion de la distorsion, el programa calcula una
C redistribucion de la senhal en cada pixel (ajustando un polinomio local de
C segundo grado), como si cada pixel se subdividiera en infinitos pixels mas
C pequenhos. Utilizando estos ajustes locales se realiza la correccion
C definitiva.
C
        PROGRAM SDISTOR
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER TRUELEN
        INTEGER READI
        INTEGER READILIM
        REAL READF
C
        INTEGER NMAX_REG
        PARAMETER(NMAX_REG=100)
C
        INTEGER I,J,K,L,IMAX
        INTEGER NS0,NS1,NS2,NW,NS0N
        INTEGER NC1,NC2,NCC1,NCC2
        INTEGER NBIN
        INTEGER NREG,NF,NFSIG,NFSIG_OLD,NFINI,NCAU
        INTEGER NTERMS,NCOLOR
        INTEGER I1,I2,II,NBEG,NEND
        INTEGER PLOTNS
        INTEGER NCP1,NCP2,IREG,NREG1(NMAX_REG),NREG2(NMAX_REG),NCHFIT
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        INTEGER NITERTOT
        REAL S(NCMAX,NSMAX),ES(NCMAX,NSMAX),S0(NCMAX,NSMAX)
        REAL SS(NSMAX)
        REAL XF(NCMAX),YF(NCMAX),XFSIG(NCMAX),YFSIG(NCMAX)
        REAL XFF(NCMAX),YFF(NCMAX)
        REAL XP(NCMAX),YP(NCMAX),YPFIN(NCMAX)
        REAL XT(NSMAX),YT(NSMAX)
        REAL XCAU(NSMAX),YCAU(NSMAX)
        REAL SIGMAS(NCMAX),AMPS(NCMAX)
        REAL YCAUMAX
        REAL AMP,XBAR,SIGMA,YOFF
        REAL AMP0,XBAR0,SIGMA0,YOFF0
        REAL X,POL,DYP,DYC
        REAL A(20),CHISQR
        REAL Z1,Z2,Z3
        REAL XMIN,XMAX,YMIN,YMAX,DX,DY
        REAL XMIN0,XMAX0,YMIN0,YMAX0
        REAL YMAXF
        REAL TSIGMA,SCANRADIUS
        REAL XPL,FNS
        REAL CCC1(NSMAX),CCC2(NSMAX),CCC3(NSMAX)         !coef. ajustes locales
        REAL FACTOR1,FACTOR2,X1,X2
        CHARACTER*1 CBIN,CPCAU,COTHER,CERR,CSEARCH,CSAVE,C2D,CIPOL
        CHARACTER*50 CDUMMY
        CHARACTER*75 INFILE,ERRFILE,OUTFILE
        CHARACTER*1 FITYPE,CSCALE,CSCALEY,FITCTE
        LOGICAL IFCHAN(NCMAX)
        LOGICAL FIRSTCAU,LITER,FITCT
        LOGICAL LCOLOR(MAX_ID_RED)
        LOGICAL L2D,LIPOL
C
        COMMON/BLKCAU1/NCAU
        COMMON/BLKCAU2/XCAU,YCAU
        COMMON/BLKFIT1/NF
        COMMON/BLKFIT2/XF,YF
        COMMON/BLKFIT3/S
        COMMON/BLKGFIT/AMP0,XBAR0,SIGMA0,YOFF0
        COMMON/BLKPLOT1/NTERM,IDN
        COMMON/BLKPLOT2/LCOLOR
        COMMON/BLKPLOT3/XMIN,XMAX,YMIN,YMAX
C------------------------------------------------------------------------------
        THISPROGRAM='sdistor'
        CALL WELCOME('10-March-2005')
C------------------------------------------------------------------------------
        NCC1=0                                    !evita warning de compilacion
        NCC2=0                                    !evita warning de compilacion
        IREG=0                                    !evita warning de compilacion
        FIRSTCAU=.FALSE.                          !evita warning de compilacion
        FITCT=.FALSE.                             !evita warning de compilacion
        NBIN=1
C
        WRITE(*,100)'Work with error images (y/n) '
        CERR(1:1)=READC('y','yn')
C
        NCOLOR=2
        LITER=.FALSE.
        PLOTNS=4
C
        CALL PIDEGTER(NTERM,IDN,LCOLOR)
C
        WRITE(*,100)'Input file name......'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        IF(CERR.EQ.'y') CALL GUESSEF(INFILE,ERRFILE)
        IF(TRUELEN(OBJECT).GT.0)THEN
          INFILE=INFILE(1:TRUELEN(INFILE))//' ['//
     +     OBJECT(1:TRUELEN(OBJECT))//']'
        END IF
        DO I=1,NSCAN
          READ(20) (S(J,I),J=1,NCHAN)
        END DO
        CLOSE(20)
        DO I=1,NSCAN
          DO J=1,NCHAN
            S0(J,I)=S(J,I)                        !duplicamos la imagen inicial
          END DO
        END DO
        IF(CERR.EQ.'y')THEN
          WRITE(*,100)'Input error file name '
          ERRFILE=INFILEX(21,ERRFILE,NSCAN,NCHAN,STWV,DISP,21,.TRUE.) !...match
          DO I=1,NSCAN
            READ(21) (ES(J,I),J=1,NCHAN)
          END DO
          CLOSE(21)
        END IF
C------------------------------------------------------------------------------
400     WRITE(*,101)'Perform fit to:'
        WRITE(*,101)'    1 : Simmetric central region'
        WRITE(*,101)'    2 : More general region'
        WRITE(*,101)'    3 : Use user-defined polynomial for distortion'
        WRITE(*,100)'Option (1/2/3) '
        FITYPE(1:1)=READC('1','123')
        IF(FITYPE.EQ.'3') GOTO 33
C dibujamos corte en la direccion espacial para buscar el maximo
        WRITE(*,100)'Channel region to plot averaged spatial profile '
        WRITE(CDUMMY,'(A2,I10)')'1,',NCHAN
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        CALL READ2I(CDUMMY(1:L),NC1,NC2)
        IMAX=0
        YMAXF=0.0
        DO I=1,NSCAN
          YT(I)=0.
          XT(I)=REAL(I)
          DO J=NC1,NC2
            YT(I)=YT(I)+S(J,I)/REAL(NC2-NC1+1)
          END DO
          IF(YT(I).GT.YMAXF)THEN
            YMAXF=YT(I)
            IMAX=I
          END IF
        END DO
        XMIN=XT(1)
        XMAX=XT(NSCAN)
        CALL FINDMM(NSCAN,YT,YMIN,YMAX)
        DX=XMAX-XMIN
        DY=YMAX-YMIN
        XMIN0=XMIN-DX/30.
        XMAX0=XMAX+DX/30.
        YMIN0=YMIN-DY/30.
        YMAX0=YMAX+DY/30.
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGSCH(1.3)
        END DO
50      DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGENV(XMIN0,XMAX0,YMIN0,YMAX0,0,0)
          CALL PGIDEN_RED
          CALL PGBIN(NSCAN,XT,YT,.TRUE.)
          CALL PGLABEL('scan','No. counts','file: '//INFILE)
        END DO
C------------------------------------------------------------------------------
        IF(FITYPE.EQ.'2') THEN
          WRITE(*,100)'Change x scale (y/n) '
          CSCALE(1:1)=READC('y','yn')
          IF(CSCALE.EQ.'y') THEN
51          WRITE(*,100)'Enter Xmin,Xmax '
            CALL READ2I('0,0',NCP1,NCP2)
            IF((NCP1.LT.1).OR.(NCP1.GT.NCP2).OR.(NCP2.GT.NSCAN))THEN
              WRITE(*,101)'ERROR: Invalid entry. Try again.'
              GOTO 51
            END IF
            XMIN0=REAL(NCP1)
            XMAX0=REAL(NCP2)
            WRITE(*,100)'Y scale: [a]utomatic, [s]ame, [k]eyboard '//
     +       '(a/s/k) '
            CSCALEY(1:1)=READC('s','ask')
            IF(CSCALEY.EQ.'a')THEN
              YMIN0=YT(NCP1)
              YMAX0=YMIN0
              IF(NCP2.GT.NCP1)THEN
                DO I=NCP1+1,NCP2
                  IF(YT(I).LT.YMIN0) YMIN0=YT(I)
                  IF(YT(I).GT.YMAX0) YMAX0=YT(I)
                END DO
              END IF
              DY=YMAX0-YMIN0
              YMIN0=YMIN0-DY/30.
              YMAX0=YMAX0+DY/30.
            ELSEIF(CSCALEY.EQ.'k')THEN
              WRITE(CDUMMY,*) YMIN0
              WRITE(*,100)'Ymin '
              YMIN0=READF(CDUMMY)
              WRITE(CDUMMY,*) YMAX0
              WRITE(*,100)'Ymax '
              YMAX0=READF(CDUMMY)
            END IF
            GOTO 50
          END IF
          IREG=0
          WRITE(*,*)
52        WRITE(*,100)'Enter scan region to fit (0,0=exit) '
          CALL READ2I('0,0',NS1,NS2)
          IF(NS1.EQ.0.AND.NS2.EQ.0) GOTO 53
          IF((NS1.LT.1).OR.(NS2.GT.NSCAN).OR.NS1.GT.NS2)THEN
            WRITE(*,101)'ERROR: Invalid entry. Try again.'
            GOTO 52
          END IF
          IREG=IREG+1
          IF(IREG.GT.NMAX_REG)THEN
            WRITE(*,101)'FATAL ERROR: number or regions .gt. NMAX_REG'
            STOP
          END IF
          NREG1(IREG)=NS1
          NREG2(IREG)=NS2
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(2)
            CALL PGSLS(2)
            CALL PGMOVE(REAL(NS1),YMIN0)
            CALL PGDRAW(REAL(NS1),YMAX0)
            CALL PGMOVE(REAL(NS2),YMIN0)
            CALL PGDRAW(REAL(NS2),YMAX0)
            CALL PGSLS(1)
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END DO
          GOTO 52
53        IF(IREG.EQ.0) GOTO 52
          NCAU=0
          DO L=1,IREG
            DO I=NREG1(L),NREG2(L)
              NCAU=NCAU+1
            END DO
          END DO
          IF(NCAU.LT.3)THEN
            WRITE(*,101)'ERROR: Insufficient number of scans for fit.'
            WRITE(*,101)'       Redefine scan region to be employed.'
            IREG=0
            GOTO 52
          END IF
          NS0=(NREG1(1)+NREG2(IREG))/2.       !valor aproximado de scan central
C..............................................................................
        ELSE
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(2)
            CALL PGRECT(XT(IMAX)-.5,XT(IMAX)+.5,0.0,YMAXF)
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END DO
          WRITE(*,110)'>>> Maximum is located at scan: ',IMAX
C
3         WRITE(*,100)'Search for a different maximum (y/n) '
          CSEARCH(1:1)=READC('n','yn')
          IF(CSEARCH.EQ.'y')THEN
            WRITE(*,101)'Define scan region: '
            WRITE(*,100)'1st scan '
            NS1=READILIM('@',1,NSCAN)
            WRITE(*,100)'2nd scan '
            NS2=READILIM('@',NS1,NSCAN)
            IMAX=0
            YMAXF=0.0
            DO I=NS1,NS2
              IF(YT(I).GT.YMAXF)THEN
                YMAXF=YT(I)
                IMAX=I
              END IF
            END DO
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              IF(LCOLOR(ITERM)) CALL PGSCI(2)
              CALL PGRECT(XT(IMAX)-.5,XT(IMAX)+.5,0.0,YMAXF)
              IF(LCOLOR(ITERM)) CALL PGSCI(1)
            END DO
            WRITE(*,110)'>>> Maximum is located at scan: ',IMAX
            GOTO 3
          END IF
C
4         WRITE(*,100)'Center scan '
          WRITE(CDUMMY,*) IMAX
          NS0=READILIM(CDUMMY,1,NSCAN)
          IMAX=NS0
5         WRITE(*,100)'No. of scans at each side to find/fit maximum '
          NW=READI('3')
          IF(NW.LT.1)THEN
            WRITE(*,101)'ERROR: sorry, minimum value is 1. Try again.'
            GOTO 5
          END IF
          NS1=NS0-NW
          NS2=NS0+NW
          IF((NS1.LT.1).OR.(NS2.GT.NSCAN))THEN
            WRITE(*,101)'ERROR: range out of frame. Try again.'
            GOTO 4
          END IF
        END IF
C------------------------------------------------------------------------------
        DO J=1,NCHAN
          IFCHAN(J)=.FALSE.
        END DO
9       NREG=0
        WRITE(*,101)'Enter channels to be fitted (0,0=EXIT):'
10      WRITE(*,100)'Channel region '
        CALL READ2I('0,0',NC1,NC2)
        IF((NC1.EQ.0).AND.(NC2.EQ.0)) GOTO 12
        IF((NC1.LT.1).OR.(NC2.GT.NCHAN).OR.(NC2.LT.NC1))THEN
          WRITE(*,101)'ERROR: Invalid entry. Try again.'
          GOTO 10
        END IF
        NCC1=NC1
        NCC2=NC2
        DO J=NC1,NC2
          IFCHAN(J)=.TRUE.
        END DO
        NREG=NREG+1
        GOTO 10
C
12      IF(NREG.EQ.0)THEN
          WRITE(*,101)'ERROR: No. of channels for fit = 0. Try again.'
          GOTO 9
        ELSEIF(NREG.EQ.1)THEN
          WRITE(*,100)'Binning (y/n) '
          CBIN(1:1)=READC('n','yn')
          IF(CBIN.EQ.'y')THEN
            WRITE(*,100)'Bin width'
            WRITE(CDUMMY,*) NBIN
            NBIN=READI(CDUMMY)
          ELSE
            NBIN=1
          END IF
        ELSE
          NBIN=1
        END IF
C
        WRITE(*,100)'No. of scans around center scan to be plotted '
        WRITE(CDUMMY,*) PLOTNS
        PLOTNS=READILIM(CDUMMY,1,NSCAN)
C
        WRITE(*,100)'Plot fits to Cauchy functions for each '//
     +   'channel (y/n) '
        CPCAU(1:1)=READC('n','yn')
        IF(CPCAU.EQ.'y') FIRSTCAU=.TRUE.
C------------------------------------------------------------------------------
C inicializamos polinomio final de distorsion
        DO J=1,NCHAN
          YPFIN(J)=0.
        END DO
C------------------------------------------------------------------------------
C Ajustamos a los canales correspondientes una funcion de Cauchy para
C encontrar los maximos. El tratamiento es diferente dependiendo de si
C se hace binning o no.
        NCHFIT=0
        DO J=1,NCHAN
          IF(IFCHAN(J)) NCHFIT=NCHFIT+1
        END DO
C hacemos un ajuste inicial al perfil promedio (sin corregir de distorsion)
        IF(FITYPE.EQ.'1') THEN
          DO I=NS1,NS2
            XCAU(I-NS1+1)=REAL(I-NS0)
            YCAU(I-NS1+1)=0.
            DO J=1,NCHAN
              IF(IFCHAN(J)) YCAU(I-NS1+1)=YCAU(I-NS1+1)+S(J,I)
            END DO
            YCAU(I-NS1+1)=YCAU(I-NS1+1)/REAL(NCHFIT)    
          END DO
          FNS=1.
          NCAU=NS2-NS1+1
        ELSE
          WRITE(*,100) 'Central scan '
          WRITE(CDUMMY,*) NS0
          NS0=READILIM(CDUMMY,1,NSCAN)
          NS1=NS0-NREG1(1)
          NS2=NREG2(IREG)-NS0
          FNS=REAL(NS1)
          IF(REAL(NS2).GT.FNS) FNS=REAL(NS2)
ccc       type*,'fns=',fns
          NCAU=0
          DO L=1,IREG
            DO I=NREG1(L),NREG2(L)
              NCAU=NCAU+1
              XCAU(NCAU)=REAL(I-NS0)/FNS
              YCAU(NCAU)=0.
              DO J=1,NCHAN
                IF(IFCHAN(J)) YCAU(NCAU)=YCAU(NCAU)+S(J,I)
              END DO
              YCAU(NCAU)=YCAU(NCAU)/REAL(NCHFIT)    
            END DO
          END DO
        END IF
C       Total fit (not used afterwards)
        YCAUMAX=YCAU(1)
        DO I=2,NCAU
          IF(YCAU(I).GT.YCAUMAX) YCAUMAX=YCAU(I)
        END DO
        IF(YCAUMAX.EQ.0.) STOP 'FATAL ERROR: division by zero.'
        DO I=1,NCAU                                               !normalizamos
          YCAU(I)=YCAU(I)/YCAUMAX
        END DO
        FITCT=.FALSE.                   !determinamos si ajustados Cauchy + cte
        IF(FITYPE.EQ.'2') THEN
          IF(NCAU.GE.4)THEN
            WRITE(*,100)'Add cte to Cauchy function (y/n) '
            FITCTE(1:1)=READC('n','yn')
            IF(FITCTE.EQ.'y') FITCT=.TRUE.
          END IF
        END IF
        IF(FITCT) THEN
          CALL GFIT2(AMP,XBAR,SIGMA,YOFF,0)  !ajustamos funcion de Cauchy + cte
          AMP0=AMP
          XBAR0=XBAR
          SIGMA0=SIGMA
          YOFF0=YOFF
        ELSE
          CALL GFIT(AMP,XBAR,SIGMA)    !ajustamos funcion de Cauchy
          YOFF=0.
        END IF
ccc     WRITE(*,101)'AMP,XBAR,SIGMA,YOFF'
ccc     WRITE(*,*)AMP,XBAR,SIGMA,YOFF
        DO I=1,NCMAX
          XP(I)=REAL(I-1)/REAL(NCMAX-1)*(XMAX-XMIN)+XMIN
          XPL=(XP(I)-REAL(NS0))/FNS
          YP(I)=AMP/(SIGMA*SIGMA+(XPL-XBAR)*(XPL-XBAR))+YOFF
          YP(I)=YP(I)*YCAUMAX
        END DO
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM)) CALL PGSCI(3)
          CALL PGLINE(1000,XP,YP)
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
        END DO
        NS0N=INT(XBAR*FNS+REAL(NS0)+0.5)
        WRITE(*,100)'>>> Original center scan: '
        WRITE(*,*) NS0
        WRITE(*,100)'>>> Fitted center scan..: '
        WRITE(*,*)NS0N
        IF(NS0.NE.NS0N) THEN       !recalculamos las distancias al scan central
          NS0=NS0N
          IF(FITYPE.EQ.'1') THEN
            DO I=NS1,NS2
              XCAU(I-NS1+1)=REAL(I-NS0)
            END DO
          ELSE
            NCAU=0
            DO L=1,IREG
              DO I=NREG1(L),NREG2(L)
                NCAU=NCAU+1
                XCAU(NCAU)=REAL(I-NS0)/FNS
              END DO
            END DO
          END IF
        END IF
C------------------------------------------------------------------------------
ccc     WRITE(*,100)'Press RETURN to continue...'
ccc     READ(*,*)
C------------------------------------------------------------------------------
20      WRITE(*,101)'Working...'
C---------------------------------------  1) SIN BINNING
        IF(NBIN.EQ.1)THEN
          K=0
          DO J=1,NCHAN
            IF(IFCHAN(J)) THEN    !solo canales seleccionados
              IF(FITYPE.EQ.'1') THEN
                DO I=NS1,NS2
                  YCAU(I-NS1+1)=S(J,I)
                END DO
              ELSE
                NCAU=0
                DO L=1,IREG
                  DO I=NREG1(L),NREG2(L)
                    NCAU=NCAU+1
                    YCAU(NCAU)=S(J,I)
                  END DO
                END DO
              END IF
              YCAUMAX=YCAU(1)
              DO I=2,NCAU
                IF(YCAU(I).GT.YCAUMAX) YCAUMAX=YCAU(I)
              END DO
              IF(YCAUMAX.EQ.0.) STOP 'FATAL ERROR: division by zero.'
              DO I=1,NCAU
                YCAU(I)=YCAU(I)/YCAUMAX
              END DO
              IF(FITCT) THEN
                CALL GFIT2(AMP,XBAR,SIGMA,YOFF,1)      !funcion de Cauchy + cte
              ELSE
                CALL GFIT(AMP,XBAR,SIGMA)                    !funcion de Cauchy
              END IF
              IF(CPCAU.EQ.'y')THEN
                IF(FIRSTCAU)THEN
                  FIRSTCAU=.FALSE.
                ELSE
                  WRITE(*,100)'More plots (y/n) '
                  CPCAU(1:1)=READC('y','yn')
                END IF
                IF(CPCAU.EQ.'y')THEN
                  DO ITERM=NTERM,1,-1
                    CALL PGSLCT(IDN(ITERM))
                    CALL PLOTCAU(AMP,XBAR,SIGMA,YOFF,J,0,NS0,FNS,
     +               LCOLOR(ITERM))
                  END DO
                ELSE
                  WRITE(*,101)'Working...'
                END IF
              END IF
              K=K+1
              XF(K)=REAL(J)
              YF(K)=XBAR*FNS+REAL(NS0)
              SIGMAS(K)=SIGMA
              AMPS(K)=AMP/(SIGMA*SIGMA)
            END IF
          END DO
          NF=K                                  !No. puntos para ajustar pol.
C---------------------------------------  2) CON BINNING
        ELSE
          NF=(NCC2-NCC1+1)/NBIN
          IF(MOD(NCC2-NCC1+1,NBIN).NE.0) NF=NF+1
          DO K=1,NF
            I1=(K-1)*NBIN+NCC1
            I2=I1+NBIN-1
            IF(I2.GT.NCC2) I2=NCC2
            IF(FITYPE.EQ.'1') THEN
              DO I=NS1,NS2
                YCAU(I-NS1+1)=0.
                DO J=I1,I2
                  YCAU(I-NS1+1)=YCAU(I-NS1+1)+S(J,I)
                END DO
              END DO
            ELSE
              NCAU=0
              DO L=1,IREG
                DO I=NREG1(L),NREG2(L)
                  NCAU=NCAU+1
                  YCAU(NCAU)=0.
                  DO J=I1,I2
                    YCAU(NCAU)=YCAU(NCAU)+S(J,I)
                  END DO
                END DO
              END DO
            END IF
            YCAUMAX=YCAU(1)
            DO I=2,NCAU
              IF(YCAU(I).GT.YCAUMAX) YCAUMAX=YCAU(I)
            END DO
            IF(YCAUMAX.EQ.0.) STOP 'FATAL ERROR: division by zero.'
            DO I=1,NCAU                                 !normalizamos
              YCAU(I)=YCAU(I)/YCAUMAX
            END DO
            IF(FITCT) THEN
              CALL GFIT2(AMP,XBAR,SIGMA,YOFF,1)        !funcion de Cauchy + cte
            ELSE
              CALL GFIT(AMP,XBAR,SIGMA)                      !funcion de Cauchy
            END IF
            IF(CPCAU.EQ.'y')THEN
              IF(FIRSTCAU)THEN
                FIRSTCAU=.FALSE.
              ELSE
                WRITE(*,100)'More plots (y/n) '
                CPCAU(1:1)=READC('y','yn')
              END IF
              IF(CPCAU.EQ.'y')THEN
                DO ITERM=NTERM,1,-1
                  CALL PGSLCT(IDN(ITERM))
                  CALL PLOTCAU(AMP,XBAR,SIGMA,YOFF,I1,I2,NS0,FNS,
     +             LCOLOR(ITERM))
                END DO
              ELSE
                WRITE(*,101)'Working...'
              END IF
            END IF
            XF(K)=REAL(I1+I2)/2.
            YF(K)=XBAR*FNS+REAL(NS0)
            SIGMAS(K)=SIGMA
            AMPS(K)=AMP/(SIGMA*SIGMA)
          END DO
        END IF
C------------------------------------------------------------------------------
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGSUBP(1,3)
ccc       CALL PGSCH(2.2)
        END DO
        XMIN=1.
        XMAX=REAL(NCHAN)
        DX=XMAX-XMIN
        XMIN=XMIN-DX/30.
        XMAX=XMAX+DX/30.
C dibujamos los valores de sigma
        CALL FINDMM(NF,SIGMAS,YMIN,YMAX)
        DY=YMAX-YMIN
        YMIN=YMIN-DY/30.
        YMAX=YMAX+DY/30.
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,0)
          CALL PGPOINT(NF,XF,SIGMAS,1)
          CALL PGLABEL('Channel','\\gs',' ')
        END DO
C dibujamos los valores de amp
        CALL FINDMM(NF,AMPS,YMIN,YMAX)
        DY=YMAX-YMIN
        YMIN=YMIN-DY/30.
        YMAX=YMAX+DY/30.
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,0)
          CALL PGPOINT(NF,XF,AMPS,1)
          CALL PGLABEL('Channel','A/\\gs\\u2',' ')
        END DO
C dibujamos los valores centrales
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PLOTFIT(NCHAN,NSCAN,NS0,PLOTNS,INFILE,LCOLOR(ITERM),1)
        END DO
        WRITE(*,100)'Press RETURN to continue...'
        READ(*,*)
C redibujamos centros y calculamos ajuste polinomico
        NFINI=NF
30      DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
ccc       CALL PGSCH(1.3)
          CALL PGSUBP(1,1)
ccc     type*,'linea 30, NF=',NF
          CALL PLOTFIT(NCHAN,NSCAN,NS0,PLOTNS,INFILE,LCOLOR(ITERM),2)
        END DO
        DYC=0.
        DO I=1,NF
          DYC=DYC+ABS(YF(I)-REAL(NS0))
        END DO
        DYC=DYC/REAL(NF)
        WRITE(*,100)'Mean deviation around central scan (scans): '
        WRITE(*,*)DYC
31      IF(LITER)THEN
          WRITE(*,100)'Polynomial degree (<20,-1=NO MORE FITS,'//
     +     '-3=QUIT)'
          NTERMS=READI('@')
          IF((NTERMS.GT.19).OR.(NTERMS.LT.-3)) GOTO 30
          IF(NTERMS.EQ.-2) GOTO 31
        ELSE
          WRITE(*,100)'Polynomial degree (<20,-2=RESTART,-3=QUIT)'
          NTERMS=READI('@')
          IF((NTERMS.GT.19).OR.(NTERMS.LT.-3)) GOTO 30
          IF(NTERMS.EQ.-1) GOTO 31
        END IF
        IF(NTERMS.EQ.-1)THEN
          GOTO 40
        ELSEIF(NTERMS.EQ.-2)THEN
          GOTO 400
        ELSEIF(NTERMS.EQ.-3)THEN
          CALL PGEND
          STOP
        END IF
        NTERMS=NTERMS+1
C eliminamos puntos muy alejados del centro
        WRITE(*,100)'No. of SCANS from center scan to exclude points '
        WRITE(CDUMMY,*)NW
        SCANRADIUS=READF(CDUMMY)
        NFSIG=0
ccc     type*,'antes de eliminar puntos lejanos, NF=',NF
        DO I=1,NF
          IF(ABS(YF(I)-REAL(NS0)).LE.SCANRADIUS)THEN
            NFSIG=NFSIG+1
            XFSIG(NFSIG)=XF(I)
            YFSIG(NFSIG)=YF(I)
          END IF
        END DO
        WRITE(*,110)'No. of points too far from the center scan: ',
     +   NF-NFSIG
        NF=NFSIG
ccc     type*,'despues de eliminar puntos lejanos, NF=',NF
        DO I=1,NFSIG
          XFF(I)=XFSIG(I)
          YFF(I)=YFSIG(I)
        END DO
C determinamos el numero de veces sigma para eliminar puntos
        WRITE(*,100)'Times SIGMA to exclude points '
        TSIGMA=READF('3.0')
        CALL POLFIT(XFF,YFF,YFF,NF,NTERMS,0,A,CHISQR)
C calculamos SIGMA alrededor del ajuste
        NFSIG_OLD=NF
32      DYP=0.
        DO I=1,NF
          X=XFF(I)
          POL=A(NTERMS)
          DO K=NTERMS-1,1,-1
            POL=POL*X+A(K)
          END DO
          DYP=DYP+(YFF(I)-POL)*(YFF(I)-POL)
        END DO
        DYP=SQRT(DYP/REAL(NF-1))
C repetimos ajuste eliminando los puntos muy alejados del ajuste
        NFSIG=0
        DO I=1,NF
          X=XFF(I)
          POL=A(NTERMS)
          DO K=NTERMS-1,1,-1
            POL=POL*X+A(K)
          END DO
          IF(ABS(YFF(I)-POL).LE.TSIGMA*DYP)THEN
            NFSIG=NFSIG+1
            XFSIG(NFSIG)=XFF(I)
            YFSIG(NFSIG)=YFF(I)
          ELSE
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              IF(LCOLOR(ITERM)) CALL PGSCI(2)
              CALL PGPOINT(1,XFF(I),YFF(I),5)
              IF(LCOLOR(ITERM)) CALL PGSCI(1)
            END DO
          END IF
        END DO
C si hemos eliminado puntos, repetimos el ajuste
        IF(NFSIG_OLD.NE.NFSIG)THEN
          CALL POLFIT(XFSIG,YFSIG,YFSIG,NFSIG,NTERMS,0,A,CHISQR)
          NFSIG_OLD=NFSIG
          GOTO 32
        END IF
C mostramos los coeficientes finales del ajuste
        WRITE(*,110)'No. of points removed from fit: ',NF-NFSIG
        WRITE(*,101)'Coefficients from fit:'
        DO I=1,NTERMS
          IF(I.LT.11)THEN
            WRITE(*,'(A,I1,A,$)')'> a(',I-1,') = '
          ELSE
            WRITE(*,'(A,I2,A,$)')'> a(',I-1,')= '
          END IF
          WRITE(*,*)A(I)
        END DO
        DYP=0.
        DO I=1,NF
          X=XFF(I)
          POL=A(NTERMS)
          DO K=NTERMS-1,1,-1
            POL=POL*X+A(K)
          END DO
          DYP=DYP+(YFF(I)-POL)*(YFF(I)-POL)
        END DO
        DYP=SQRT(DYP/REAL(NF-1))
        WRITE(*,100)'Mean dispersion around polynomial  : '
        WRITE(*,*)DYP
        WRITE(*,*)
C dibujamos el ajuste realizado
        DO I=1,NCHAN
          XP(I)=REAL(I)
          YP(I)=A(NTERMS)
          DO K=NTERMS-1,1,-1
            YP(I)=YP(I)*XP(I)+A(K)
          END DO
        END DO
        NCOLOR=NCOLOR+1
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM)) CALL PGSCI(NCOLOR)
          CALL PGLINE(NCHAN,XP,YP)
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
        END DO
        WRITE(*,100)'Change polynomial degree (y/n) '
        COTHER(1:1)=READC('n','yn')
        IF(COTHER.EQ.'y')THEN
          NF=NFINI
          GOTO 30
        END IF
C------------------------------------------------------------------------------
C realizamos la correccion de distorsion S (duplicamos el codigo de COSDIS3)
        DO J=1,NCHAN
          YP(J)=YP(J)-REAL(NS0)
          YPFIN(J)=YPFIN(J)+YP(J)           !acumulamos los ajustes polinomicos
        END DO
C
33      IF(FITYPE.EQ.'3')THEN
          WRITE(*,100) 'Polynomial degree '
          NTERMS=READILIM('1',1,19)
          NTERMS=NTERMS+1
          DO K=1,NTERMS
            WRITE(*,'(A2,I2.2,A2,$)') 'a(',K,')='
            A(K)=READF('@')
          END DO
          DO J=1,NCHAN
            XP(J)=REAL(J)
            YP(J)=A(NTERMS)
            DO K=NTERMS-1,1,-1
              YP(J)=YP(J)*XP(J)+A(K)
            END DO
            YPFIN(J)=YP(J)
          END DO
        END IF
C
        DO J=1,NCHAN
          DO I=1,NSCAN
            SS(I)=0.
          END DO
          Z3=ABS(YP(J))
          Z2=Z3-AINT(Z3)
          Z1=1.-Z2
          I1=INT(YP(J))
          I2=I1+INT(SIGN(1.,YP(J)))
          II=1-I2
          NBEG=MAX0(1,II)
          II=NSCAN-I2
          NEND=MIN0(NSCAN,II)
          DO I=NBEG,NEND
            SS(I)=Z1*S(J,I+I1)+Z2*S(J,I+I2)
          END DO
          DO I=1,NSCAN
            S(J,I)=SS(I)
          END DO
        END DO
        IF(FITYPE.EQ.'3') GOTO 40 !....................no iteramos en este caso
C------------------------------------------------------------------------------
        LITER=.TRUE.
        GOTO 20
C------------------------------------------------------------------------------
40      CONTINUE
        WRITE(*,100)'Save input file after initial correction (y/n) '
        CSAVE(1:1)=READC('n','yn')
        IF(CSAVE.EQ.'y')THEN
          WRITE(*,100)'Output file name'
          OUTFILE=OUTFILEX(30,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
          DO I=1,NSCAN
            WRITE(30) (S(J,I),J=1,NCHAN)
          END DO
          CLOSE(30)
        END IF
C------------------------------------------------------------------------------
C dibujamos los ajuste locales en un canal seleccionado
60      WRITE(*,*)
        WRITE(*,101)'* Plotting fits to local polynomials:'
        WRITE(*,100)'Channel (0=EXIT)'
        J=READILIM('@',0,NCHAN)
        IF(J.EQ.0) GOTO 70
        WRITE(*,100)'No. of iterations to link local fits '
        NITERTOT=READILIM('0',0,100)
        DO I=1,NSCAN
          XT(I)=REAL(I)
          YT(I)=S0(J,I)
        END DO
        IF(NITERTOT.GT.0)THEN
          WRITE(*,100)'Prevent sign change in second order derivates '
          WRITE(*,100)'(y/n) '
          C2D(1:1)=READC('n','yn')
          L2D=(C2D.EQ.'y')
          WRITE(*,100)'Plot intermediate polynomials (y/n) '
          CIPOL(1:1)=READC('y','yn')
          LIPOL=(CIPOL.EQ.'y')
        ELSE
          L2D=.FALSE.
          LIPOL=.TRUE.
        END IF
C dibujamos corte espacial en el canal J
        XMIN=1.
        XMAX=REAL(NSCAN)
        DX=XMAX-XMIN
        XMIN=XMIN-DX/30.
        XMAX=XMAX+DX/30.
        CALL FINDMM(NSCAN,YT,YMIN,YMAX)
        DY=YMAX-YMIN
        YMIN=YMIN-DY/30.
        YMAX=YMAX+DY/30.
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGSLW(3)
          CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,0)
          CALL PGIDEN_RED
          CALL PGBIN(NSCAN,XT,YT,.TRUE.)
          CALL PGLABEL('scan','No. counts','File: '//INFILE)
          CALL PGSLW(1)
        END DO
C cambiamos limites en X
        WRITE(*,100)'X zoom: ns1,ns2'
        CALL READ2I('@',NS1,NS2)
        XMIN=REAL(NS1)
        XMAX=REAL(NS2)
        DX=XMAX-XMIN
        XMIN=XMIN-DX/30.
        XMAX=XMAX+DX/30.
        CALL FINDMML(NSCAN,NS1,NS2,YT,YMIN,YMAX)
        DY=YMAX-YMIN
        YMIN=YMIN-DY/10.
        YMAX=YMAX+DY/10.
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGSLW(3)
          CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,0)
          CALL PGIDEN_RED
          CALL PGBIN(NSCAN,XT,YT,.TRUE.)
          CALL PGLABEL('scan','No. counts','File: '//INFILE)
          CALL PGSLW(1)
        END DO
C dibujamos ajustes locales
        CALL LIGHTSCAN(NSCAN,YT,NITERTOT,CCC1,CCC2,CCC3,.TRUE.,L2D,
     +   LIPOL)
        GOTO 60
C------------------------------------------------------------------------------
70      WRITE(*,100)'No. of iterations to link local fits '
        NITERTOT=READILIM('0',0,100)
        IF(NITERTOT.GT.0)THEN
          WRITE(*,100)'Prevent sign change in second order derivates '
          WRITE(*,100)'(y/n) '
          C2D(1:1)=READC('n','yn')
          L2D=(C2D.EQ.'y')
          LIPOL=.FALSE.
        ELSE
          L2D=.FALSE.
          LIPOL=.FALSE.
        END IF
C corregimos imagen inicial
        WRITE(*,100)'Working with initial frame...'
        DO I=1,NSCAN
          XT(I)=REAL(I)
        END DO
        DO J=1,NCHAN
          IF(MOD(J,100).EQ.0)THEN
            WRITE(*,*)
            WRITE(*,'(A,I6,A,$)')'channell: ',J,'...'
          END IF
          DO I=1,NSCAN
            YT(I)=S0(J,I)
          END DO
          CALL LIGHTSCAN(NSCAN,YT,NITERTOT,CCC1,CCC2,CCC3,.FALSE.,L2D,
     +     LIPOL)
          DO I=1,NSCAN
            SS(I)=0.
          END DO
          Z3=ABS(YPFIN(J))
          Z2=Z3-AINT(Z3)
          Z1=1.-Z2
          I1=INT(YPFIN(J))
          I2=I1+INT(SIGN(1.,YPFIN(J)))
          II=1-I2
          NBEG=MAX0(1,II)
          II=NSCAN-I2
          NEND=MIN0(NSCAN,II)
          DO I=NBEG,NEND
            IF(I1.LE.I2)THEN
              IF(Z1.EQ.0.0)THEN
                FACTOR1=0.
              ELSE
                X1=+0.5-Z1
                X2=+0.5
                FACTOR1=CCC3(I+I1)*(X2*X2*X2-X1*X1*X1)/3.+
     +                  CCC2(I+I1)*(X2*X2-X1*X1)/2.+
     +                  CCC1(I+I1)*(X2-X1)
                FACTOR1=FACTOR1/(X2-X1)
              END IF
              IF(Z2.EQ.0.0)THEN
                FACTOR2=0.
              ELSE
                X1=-0.5
                X2=Z2-0.5
                FACTOR2=CCC3(I+I2)*(X2*X2*X2-X1*X1*X1)/3.+
     +                  CCC2(I+I2)*(X2*X2-X1*X1)/2.+
     +                  CCC1(I+I2)*(X2-X1)
                FACTOR2=FACTOR2/(X2-X1)
              END IF
            ELSE
              IF(Z1.EQ.0.0)THEN
                FACTOR1=0.
              ELSE
                X1=-0.5
                X2=Z1-0.5
                FACTOR1=CCC3(I+I1)*(X2*X2*X2-X1*X1*X1)/3.+
     +                  CCC2(I+I1)*(X2*X2-X1*X1)/2.+
     +                  CCC1(I+I1)*(X2-X1)
                FACTOR1=FACTOR1/(X2-X1)
              END IF
              IF(Z2.EQ.0.0)THEN
                FACTOR2=0.
              ELSE
                X1=+0.5-Z2
                X2=+0.5
                FACTOR2=CCC3(I+I2)*(X2*X2*X2-X1*X1*X1)/3.+
     +                  CCC2(I+I2)*(X2*X2-X1*X1)/2.+
     +                  CCC1(I+I2)*(X2-X1)
                FACTOR2=FACTOR2/(X2-X1)
              END IF
            END IF
            SS(I)=Z1*FACTOR1+Z2*FACTOR2
          END DO
          DO I=1,NSCAN
            S(J,I)=SS(I)
          END DO
        END DO
        WRITE(*,101)'   ...OK!'
C
        IF(CERR.EQ.'n') GOTO 90
C..............................................................................
C corregimos imagen de errores
        WRITE(*,100)'Working with error frame...'
        DO J=1,NCHAN
          IF(MOD(J,100).EQ.0)THEN
            WRITE(*,*)
            WRITE(*,'(A,I6,A,$)')'channell: ',J,'...'
          END IF
          DO I=1,NSCAN
            YT(I)=ES(J,I)
          END DO
          CALL LIGHTSCAN(NSCAN,YT,NITERTOT,CCC1,CCC2,CCC3,.FALSE.,L2D,
     +     LIPOL)
          DO I=1,NSCAN
            SS(I)=0.
          END DO
          Z3=ABS(YPFIN(J))
          Z2=Z3-AINT(Z3)
          Z1=1.-Z2
          I1=INT(YPFIN(J))
          I2=I1+INT(SIGN(1.,YPFIN(J)))
          II=1-I2
          NBEG=MAX0(1,II)
          II=NSCAN-I2
          NEND=MIN0(NSCAN,II)
          DO I=NBEG,NEND
            IF(I1.LE.I2)THEN
              IF(Z1.EQ.0.0)THEN
                FACTOR1=0.
              ELSE
                X1=+0.5-Z1
                X2=+0.5
                FACTOR1=CCC3(I+I1)*(X2*X2*X2-X1*X1*X1)/3.+
     +                  CCC2(I+I1)*(X2*X2-X1*X1)/2.+
     +                  CCC1(I+I1)*(X2-X1)
                FACTOR1=FACTOR1/(X2-X1)
              END IF
              IF(Z2.EQ.0.0)THEN
                FACTOR2=0.
              ELSE
                X1=-0.5
                X2=Z2-0.5
                FACTOR2=CCC3(I+I2)*(X2*X2*X2-X1*X1*X1)/3.+
     +                  CCC2(I+I2)*(X2*X2-X1*X1)/2.+
     +                  CCC1(I+I2)*(X2-X1)
                FACTOR2=FACTOR2/(X2-X1)
              END IF
            ELSE
              IF(Z1.EQ.0.0)THEN
                FACTOR1=0.
              ELSE
                X1=-0.5
                X2=Z1-0.5
                FACTOR1=CCC3(I+I1)*(X2*X2*X2-X1*X1*X1)/3.+
     +                  CCC2(I+I1)*(X2*X2-X1*X1)/2.+
     +                  CCC1(I+I1)*(X2-X1)
                FACTOR1=FACTOR1/(X2-X1)
              END IF
              IF(Z2.EQ.0.0)THEN
                FACTOR2=0.
              ELSE
                X1=+0.5-Z2
                X2=+0.5
                FACTOR2=CCC3(I+I2)*(X2*X2*X2-X1*X1*X1)/3.+
     +                  CCC2(I+I2)*(X2*X2-X1*X1)/2.+
     +                  CCC1(I+I2)*(X2-X1)
                FACTOR2=FACTOR2/(X2-X1)
              END IF
            END IF
            SS(I)=Z1*FACTOR1+Z2*FACTOR2
          END DO
          DO I=1,NSCAN
            ES(J,I)=SS(I)
          END DO
        END DO
        WRITE(*,101)'   ...OK!'
C------------------------------------------------------------------------------
C Salvamos fichero(s)
90      WRITE(*,100)'Output file name.......(corrected from '//
     +   'S-distortion)'
        OUTFILE=OUTFILEX(30,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        DO I=1,NSCAN
          WRITE(30) (S(J,I),J=1,NCHAN)
        END DO
        CLOSE(30)
        IF(CERR.EQ.'y')THEN
          WRITE(*,100)'Output error file name (corrected from '//
     +     'S-distortion) '
          CALL GUESSEF(OUTFILE,ERRFILE)
          OUTFILE=OUTFILEX(31,ERRFILE,NSCAN,NCHAN,STWV,DISP,1,.TRUE.)
          DO I=1,NSCAN
            WRITE(31) (ES(J,I),J=1,NCHAN)
          END DO
          CLOSE(31)
        END IF
C------------------------------------------------------------------------------
        CALL PGEND
C
        STOP
100     FORMAT(A,$)
101     FORMAT(A)
110     FORMAT(A,I6)
        END
C
C******************************************************************************
C
        SUBROUTINE GFIT(AMP,XBAR,SIGMA)
        IMPLICIT NONE
C
        INTEGER NDIM,NEVAL
        REAL XX0(3),DXX0(3),XXF(3),DXXF(3)
        REAL YRMSTOL
        EXTERNAL FUNKSDISTOR
        REAL FUNKSDISTOR
        REAL AMP,XBAR,SIGMA
C------------------------------------------------------------------------------
        NDIM=3
        YRMSTOL=1.E-6
        XX0(1)=1.
        XX0(2)=0.
        XX0(3)=2.
        DXX0(1)=0.1
        DXX0(2)=0.1
        DXX0(3)=0.1
C
        CALL DOWNHILL(NDIM,XX0,DXX0,FUNKSDISTOR,1.0,0.5,2.0,YRMSTOL,
     +   XXF,DXXF,NEVAL)
C
        AMP=XXF(1)
        XBAR=XXF(2)
        SIGMA=XXF(3)
        END
C
C******************************************************************************
C
        REAL FUNCTION FUNKSDISTOR(X)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
C
        INTEGER I
        INTEGER NCAU
        REAL X(3)
        REAL XCAU(NSMAX),YCAU(NSMAX)
        DOUBLE PRECISION SM,X1,Z1
        DOUBLE PRECISION AMP,XBAR,SIGMA
C
        COMMON/BLKCAU1/NCAU
        COMMON/BLKCAU2/XCAU,YCAU
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        SM=0.D0
        AMP=DBLE(X(1))
        XBAR=DBLE(X(2))
        SIGMA=DBLE(X(3))
C
        IF(X(1).LE.0.0D+00) GO TO 10
        IF(X(3).LE.0.0D+00) GO TO 10
C
        DO I=1,NCAU
          X1=DBLE(XCAU(I))
          Z1=AMP/(SIGMA**2.D0+(X1-XBAR)**2.D0)
          SM=SM+(DBLE(YCAU(I))-Z1)**2.D0
        END DO
C
        FUNKSDISTOR=REAL(SM)
        RETURN
  10    FUNKSDISTOR=1.E+10
        END
C
C
C******************************************************************************
C
        SUBROUTINE GFIT2(AMP,XBAR,SIGMA,YOFF,IGFIT)
        IMPLICIT NONE
C
        INTEGER NDIM,NEVAL,IGFIT
        REAL YRMSTOL
        REAL XX0(4),DXX0(4),XXF(4),DXXF(4)
        EXTERNAL FUNKSDISTOR2
        REAL FUNKSDISTOR2
        REAL AMP,XBAR,SIGMA,YOFF
        REAL AMP0,XBAR0,SIGMA0,YOFF0

        COMMON/BLKGFIT/AMP0,XBAR0,SIGMA0,YOFF0

C------------------------------------------------------------------------------
        NDIM=4
        YRMSTOL=1.E-6
C
        IF(IGFIT.EQ.0) THEN
          XX0(1)=1.
          XX0(2)=0.
          XX0(3)=0.1
          XX0(4)=0.
          DXX0(1)=0.1
          DXX0(2)=0.1
          DXX0(3)=0.1
          DXX0(4)=0.1
        ELSE
          XX0(1)=AMP0
          XX0(2)=XBAR0
          XX0(3)=SIGMA0
          XX0(4)=YOFF0
          DXX0(1)=0.1
          DXX0(2)=0.1
          DXX0(3)=0.1
          DXX0(4)=0.1
        END IF
C
        CALL DOWNHILL(NDIM,XX0,DXX0,FUNKSDISTOR2,1.0,0.5,2.0,YRMSTOL,
     +   XXF,DXXF,NEVAL)
C
        AMP=XXF(1)
        XBAR=XXF(2)
        SIGMA=XXF(3)
        YOFF=XXF(4)
        END
C
C******************************************************************************
C
        REAL FUNCTION FUNKSDISTOR2(X)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
C
        INTEGER I
        INTEGER NCAU
        REAL X(4)
        REAL XCAU(NSMAX),YCAU(NSMAX)
        DOUBLE PRECISION SM,X1,Z1
        DOUBLE PRECISION AMP,XBAR,SIGMA,YOFF
C
        COMMON/BLKCAU1/NCAU
        COMMON/BLKCAU2/XCAU,YCAU
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        SM=0.D0
        AMP=DBLE(X(1))
        XBAR=DBLE(X(2))
        SIGMA=DBLE(X(3))
        YOFF=DBLE(X(4))
C
        IF(X(1).LE.0.0D+00) GO TO 10
        IF(X(3).LE.0.0D+00) GO TO 10
C
        DO I=1,NCAU
          X1=DBLE(XCAU(I))
          Z1=AMP/(SIGMA**2.D0+(X1-XBAR)**2.D0)+YOFF
          SM=SM+(DBLE(YCAU(I))-Z1)**2.D0
        END DO
C
        FUNKSDISTOR2=REAL(SM)
        RETURN
  10    FUNKSDISTOR2=1.E+10
        END
C
C******************************************************************************
C
        SUBROUTINE PLOTCAU(AMP,XBAR,SIGMA,YOFF,NC1,NC2,N0,FNS,LCOLOR)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        REAL AMP,XBAR,SIGMA,YOFF
C
        INTEGER I,NC1,NC2,N0
        INTEGER NCAU,L
        REAL XCAU(NSMAX),YCAU(NSMAX)
        REAL XMIN,XMAX,YMIN,YMAX,DX,DY
        REAL XP(1000),YP(1000)
        REAL FNS,XCAU2(NSMAX),XPL
        CHARACTER*50 CDUMMY
        LOGICAL LCOLOR
C
        COMMON/BLKCAU1/NCAU
        COMMON/BLKCAU2/XCAU,YCAU
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        XCAU2(1)=XCAU(1)*FNS
        XMIN=XCAU2(1)
        XMAX=XMIN
        YMAX=YCAU(1)
        YMIN=YMAX
        DO I=2,NCAU
          XCAU2(I)=XCAU(I)*FNS
          IF(XCAU2(I).LT.XMIN) XMIN=XCAU2(I)
          IF(XCAU2(I).GT.XMAX) XMAX=XCAU2(I)
          IF(YCAU(I).LT.YMIN) YMIN=YCAU(I)
          IF(YCAU(I).GT.YMAX) YMAX=YCAU(I)
        END DO
        IF(YMAX.LT.AMP/(SIGMA*SIGMA)+YOFF) YMAX=AMP/(SIGMA*SIGMA)+YOFF
        DX=XMAX-XMIN
        DY=YMAX-YMIN
        XMIN=XMIN-DX/50.
        XMAX=XMAX+DX/50.
        YMIN=YMIN-DY/50.
        YMAX=YMAX+DY/50.
        CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,0)
        CALL PGIDEN_RED
        CALL PGBIN(NCAU,XCAU2,YCAU,.TRUE.)
        IF(NC2.EQ.0)THEN
          WRITE(CDUMMY,110)'Channel #',NC1
          CALL PGLABEL(' ','normalized No. of counts',CDUMMY)
        ELSE
          WRITE(CDUMMY,'(I10,A1,I10)')NC1,'-',NC2
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          CALL PGLABEL(' ','normalized No. of counts',
     +     'Added Channels: '//CDUMMY(1:L))
        END IF
        WRITE(CDUMMY,*)N0
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        CALL PGLABEL('scan (relative to scan#'//CDUMMY(1:L)//')',
     +   ' ',' ')
        IF(LCOLOR) CALL PGSCI(2)
        DO I=1,1000
          XP(I)=(REAL(I-1)/999.*(XMAX-XMIN)+XMIN)
          XPL=XP(I)/FNS
          YP(I)=AMP/(SIGMA*SIGMA+(XPL-XBAR)*(XPL-XBAR))+YOFF
        END DO
        CALL PGLINE(1000,XP,YP)
        IF(LCOLOR) CALL PGSCI(1)
C
110     FORMAT(A,I6)
        END
C
C******************************************************************************
C INUM=1: dibuja imagen en escale de grises
C INUM=2: no dibuja imagen en escale de grises
        SUBROUTINE PLOTFIT(NCHAN,NSCAN,NS0,PLOTNS,INFILE,LCOLOR,INUM)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INTEGER TRUELEN
C
        INTEGER NS0,PLOTNS,INUM
        CHARACTER*(*) INFILE
        LOGICAL LCOLOR
C
        INTEGER I,J,L
        INTEGER NF
        REAL S(NCMAX,NSMAX)
        REAL XF(NCMAX),YF(NCMAX)
        REAL XMIN,XMAX,DX
        REAL TR(6),FG,BG
        CHARACTER*255 LOCALGLABEL
C
        COMMON/BLKFIT1/NF
        COMMON/BLKFIT2/XF,YF
        COMMON/BLKFIT3/S
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
        XMIN=1.
        XMAX=REAL(NCHAN)
        DX=XMAX-XMIN
        XMIN=XMIN-DX/30.
        XMAX=XMAX+DX/30.
        CALL PGENV(XMIN,XMAX,REAL(NS0-PLOTNS),REAL(NS0+PLOTNS),0,0)
        CALL PGIDEN_RED
        LOCALGLABEL='S-distortion, file='
        L=TRUELEN(INFILE)
        LOCALGLABEL(19+1:19+L)=INFILE(1:L)
        CALL PGLABEL('Channel','Scan',LOCALGLABEL(1:19+L))
        IF(INUM.EQ.1)THEN
          FG=S(1,NS0-PLOTNS)
          BG=FG
          DO I=NS0-PLOTNS,NS0+PLOTNS
            DO J=1,NCHAN
              IF(BG.GT.S(J,I)) BG=S(J,I)
              IF(FG.LT.S(J,I)) FG=S(J,I)
            END DO
          END DO
          CALL PGGRAY(S,NCMAX,NSMAX,1,NCHAN,NS0-PLOTNS,NS0+PLOTNS,
     +     FG,BG,TR)
        END IF
        IF(LCOLOR)CALL PGSCI(2)
        CALL PGSLS(2)
        DO I=NS0-PLOTNS,NS0+PLOTNS
          CALL PGMOVE(0.,REAL(I)+.5)
          CALL PGDRAW(REAL(NCHAN+1),REAL(I)+.5)
        END DO
        CALL PGSLS(1)
        IF(INUM.EQ.1)THEN
C dibujamos un hueco para ver mejor en el PostScript
          CALL PGSLW(3)
          CALL PGSCI(0)
          CALL PGLINE(NF,XF,YF)
          CALL PGSCI(1)
          CALL PGSLW(1)
        END IF
C
        IF(LCOLOR)CALL PGSCI(4)
        CALL PGPOINT(NF,XF,YF,1)
        IF(LCOLOR)CALL PGSCI(1)
        END
C
C******************************************************************************
C Distribuye la luz en cada scan utilizando ajustes locales a polinomios de
C segundo grado.
C Si LPLOT=.TRUE. dibujamos.
C Si LPLOT=.TRUE. y LIPOL=.TRUE. dibujamos polinomios intermedios
C Si LPLOT=.TRUE. y LIPOL=.FALSE. solo dibujamos polinomio final
C
        SUBROUTINE LIGHTSCAN(NS,YT,NITERTOT,CCC1,CCC2,CCC3,LPLOT,L2D,
     +   LIPOL)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INTEGER NS
        REAL YT(NS)
        INTEGER NITERTOT
        REAL CCC1(NS),CCC2(NS),CCC3(NS)
        LOGICAL LPLOT,L2D,LIPOL
C
        INTEGER NLOCALITERMAX
        PARAMETER (NLOCALITERMAX=20)
C
        INTEGER I,I1,I2,L
        INTEGER NITER
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        INTEGER SIGNO(NSMAX),SIGNO0
        INTEGER NLOCALITER
        REAL YF(3)
        REAL XP(NSMAX),YP(NSMAX)
        REAL X,POL
        REAL X1,X2,Y1,Y2
        REAL MEANPOL(NSMAX)
        REAL Y1POL(NSMAX),Y2POL(NSMAX)
        REAL PESO1(NSMAX),PESO2(NSMAX)
        REAL XMIN,XMAX,YMIN,YMAX,DY
        LOGICAL LCOLOR(MAX_ID_RED)
        LOGICAL LREPEAT
C
        COMMON/BLKPLOT1/NTERM,IDN
        COMMON/BLKPLOT2/LCOLOR
        COMMON/BLKPLOT3/XMIN,XMAX,YMIN,YMAX
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
        DY=0.0 !evita un warning de compilacion
C
        NITER=0
        IF(NS.LT.3)THEN
          WRITE(*,101)'FATAL ERROR: in subroutine LIGHTSCAN'
          WRITE(*,101)'>>> NSCAN .LT. 3!'
          STOP
        END IF
C Ajustamos localmente polinomios de segundo grado: calculamos los polinomios
C suponiendo que estos deben pasar por el centro de los pixels con un valor
C en ordenadas igual al valor del pixel.
        IF(LPLOT.AND.LIPOL)THEN
          DO I=1,NS
            IF(I.EQ.1)THEN
              I1=1
              I2=3
              X1=-1.5
              X2=-0.5
            ELSEIF(I.EQ.NS)THEN
              I2=NS
              I1=NS-2
              X1=+0.5
              X2=+1.5
            ELSE
              I1=I-1
              I2=I+1
              X1=-0.5
              X2=+0.5
            END IF
            DO L=I1,I2
              YF(L-I1+1)=YT(L)
            END DO
            CCC1(I)=YF(2)
            CCC2(I)=(YF(3)-YF(1))/2.
            CCC3(I)=(YF(3)+YF(1))/2.-YF(2)
            MEANPOL(I)=CCC3(I)*(X2*X2*X2-X1*X1*X1)/3.+
     +                 CCC2(I)*(X2*X2-X1*X1)/2.+
     +                 CCC1(I)*(X2-X1)
            IF(NINT(SIGN(1.,YT(I))).NE.NINT(SIGN(1.,MEANPOL(I))))THEN
              MEANPOL(I)=-MEANPOL(I)        !evita problemas de cambio de signo
            END IF
          END DO
C dibujamos
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(4)
            DO I=1,NS
              DO L=1,30    !dividimos cada pixel en 30 subunidades para dibujar
                X=-0.5+REAL(L-1)/29.
                XP(L)=X+REAL(I)
                IF(I.EQ.1) X=X-1.0
                IF(I.EQ.NS) X=X+1.0
                POL=CCC3(I)*X*X+CCC2(I)*X+CCC1(I)
                YP(L)=POL
              END DO
              CALL PGLINE(30,XP,YP)
            END DO
          END DO
          WRITE(*,100)'Type RETURN to continue...'
          READ(*,*)
        END IF
C..............................................................................
C polinomios iniciales escalados al numero de cuentas de cada pixel
        IF(LPLOT.AND.LIPOL)THEN
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(2)
            DO I=1,NS
              DO L=1,30    !dividimos cada pixel en 30 subunidades para dibujar
                X=-0.5+REAL(L-1)/29.
                XP(L)=X+REAL(I)
                IF(I.EQ.1) X=X-1.0
                IF(I.EQ.NS) X=X+1.0
                POL=CCC3(I)*X*X+CCC2(I)*X+CCC1(I)
                YP(L)=POL*YT(I)/MEANPOL(I)
              END DO
              CALL PGLINE(30,XP,YP)      !senhal normalizada al polinomio local
            END DO
          END DO
          WRITE(*,100)'Type RETURN to continue...'
          READ(*,*)
        END IF
C------------------------------------------------------------------------------
C calculamos polinomios locales imponiendo conservacion de area en cada 3
C pixels consecutivos
        DO I=1,NS
          IF(I.EQ.1)THEN
            I1=1
            I2=3
          ELSEIF(I.EQ.NS)THEN
            I2=NS
            I1=NS-2
          ELSE
            I1=I-1
            I2=I+1
          END IF
          DO L=I1,I2
            YF(L-I1+1)=YT(L)
          END DO
          CCC3(I)=(YF(3)+YF(1))/2.-YF(2)
          CCC2(I)=(YF(3)-YF(1))/2.
          CCC1(I)=YF(2)-CCC3(I)/12.
          SIGNO(I)=NINT(SIGN(1.,CCC3(I)))
        END DO
C dibujamos
        IF(LPLOT)THEN
          DY=YMAX-YMIN
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(7)
            DO I=1,NS
              DO L=1,30    !dividimos cada pixel en 30 subunidades para dibujar
                X=-0.5+REAL(L-1)/29.
                XP(L)=X+REAL(I)
                IF(I.EQ.1) X=X-1.0
                IF(I.EQ.NS) X=X+1.0
                POL=CCC3(I)*X*X+CCC2(I)*X+CCC1(I)
                YP(L)=POL
              END DO
              CALL PGSLS(4)
              CALL PGLINE(30,XP,YP)
              CALL PGSLS(1)
              if((real(i).gt.xmin).and.(real(i).lt.xmax))then
                if(signo(i).gt.0)then
                  call pgptext(real(i),ymin+DY/40.,0.,0.5,'+')
                elseif(signo(i).lt.0)then
                  call pgptext(real(i),ymin+DY/40.,0.,0.5,'-')
                else
                  call pgptext(real(i),ymin+DY/40.,0.,0.5,'0')
                end if
              end if
            END DO
          END DO
          WRITE(*,100)'Type RETURN to continue...'
          READ(*,*)
        END IF
C..............................................................................
C si no iteramos, volvemos al programa principal sin forzar union de 
C polinomios ni igualdad de las derivadas
        IF(NITERTOT.EQ.0)THEN
          IF(LPLOT)THEN
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              IF(LCOLOR(ITERM)) CALL PGSCI(3)
              CALL PGSLW(3)
              DO I=1,NS
                DO L=1,30  !dividimos cada pixel en 30 subunidades para dibujar
                  X=-0.5+REAL(L-1)/29.
                  XP(L)=X+REAL(I)
                  IF(I.EQ.1) X=X-1.0
                  IF(I.EQ.NS) X=X+1.0
                  POL=CCC3(I)*X*X+CCC2(I)*X+CCC1(I)
                  YP(L)=POL
                END DO
                CALL PGLINE(30,XP,YP)
              END DO
              CALL PGSLW(1)
              IF(LCOLOR(ITERM)) CALL PGSCI(1)
            END DO
          END IF
          RETURN
        END IF
C..............................................................................
C polinomios anteriores empalmados en los bordes de cada pixel
        DO I=1,NS
          X1=-0.5
          IF(I.EQ.1) X1=X1-1.0
          IF(I.EQ.NS) X1=X1+1.0
          Y1=CCC3(I)*X1*X1+CCC2(I)*X1+CCC1(I)
          Y1POL(I)=Y1
          X2=+0.5
          IF(I.EQ.1) X2=X2-1.0
          IF(I.EQ.NS) X2=X2+1.0
          Y2=CCC3(I)*X2*X2+CCC2(I)*X2+CCC1(I)
          Y2POL(I)=Y2
        END DO
C promediamos los valores en las fronteras de cada pixel, salvo en los bordes
C de la imagen y cuando se viole la conservacion de la derivada segunda 
C inicial (siempre y cuando esta conservacion sea requerida); la influencia de
C los pixels cuya derivada segunda hay que conservar se transmite tambien a los
C pixels vecinos
10      CONTINUE
        DO I=1,NS
          PESO1(I)=0.5
          PESO2(I)=0.5
        END DO
        PESO1(1)=1.0
        PESO2(NS)=1.0
        NLOCALITER=0
C
12      LREPEAT=.FALSE.
        NLOCALITER=NLOCALITER+1
        DO I=1,NS
          IF(I.EQ.1)THEN
            Y1=Y1POL(I)
          ELSE
            Y1=Y2POL(I-1)*PESO2(I-1)+Y1POL(I)*PESO1(I)
          END IF
          IF(I.EQ.NS)THEN
            Y2=Y2POL(I)
          ELSE
            Y2=Y2POL(I)*PESO2(I)+Y1POL(I+1)*PESO1(I+1)
          END IF
          IF(I.EQ.1)THEN
            CCC1(I)=1.75*Y1+3.75*Y2-4.5*YT(I)
            CCC2(I)=5.0*Y1+7.0*Y2-12.0*YT(I)
            CCC3(I)=3.0*(Y1+Y2)-6.0*YT(I)
          ELSEIF(I.EQ.NS)THEN
            CCC1(I)=3.75*Y1+1.75*Y2-4.5*YT(I)
            CCC2(I)=-7.0*Y1-5.0*Y2+12.0*YT(I)
            CCC3(I)=3.0*(Y1+Y2)-6.0*YT(I)
          ELSE
            CCC1(I)=1.5*YT(I)-0.25*(Y1+Y2)
            CCC2(I)=Y2-Y1
            CCC3(I)=3.0*(Y1+Y2)-6.0*YT(I)
          END IF
c si la derivada segunda cambia de signo, aumentamos los pesos del pixel I
C en detrimento de los pixels vecinos
          IF(L2D)THEN                 !hay que conservar las derivadas segundas
            SIGNO0=NINT(SIGN(1.,CCC3(I)))
            IF(SIGNO0.NE.SIGNO(I))THEN
              LREPEAT=.TRUE.
              IF(I.EQ.1)THEN
                PESO1(I+1)=PESO1(I+1)/2.
                PESO2(I)=1.-PESO1(I+1)
              ELSEIF(I.EQ.NS)THEN
                PESO2(I-1)=PESO2(I-1)/2.
                PESO1(I)=1.-PESO2(I-1)
              ELSE
                PESO2(I-1)=PESO2(I-1)/2.
                PESO1(I)=1.-PESO2(I-1)
                PESO1(I+1)=PESO1(I+1)/2.
                PESO2(I)=1.-PESO1(I+1)
ccc           type*,i,nlocaliter,
ccc     +         peso2(i-1)+peso1(i),peso2(i)+peso1(i+1)
              END IF
            END IF
          END IF
        END DO
        IF((LREPEAT).AND.(NLOCALITER.LE.NLOCALITERMAX)) GOTO 12
C dibujamos
        IF(LPLOT.AND.LIPOL)THEN
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(14)
            DO I=1,NS
              DO L=1,30    !dividimos cada pixel en 30 subunidades para dibujar
                X=-0.5+REAL(L-1)/29.
                XP(L)=REAL(I)+X
                IF(I.EQ.1) X=X-1.0
                IF(I.EQ.NS) X=X+1.0
                POL=CCC3(I)*X*X+CCC2(I)*X+CCC1(I)
                YP(L)=POL
              END DO
              CALL PGLINE(30,XP,YP)      !senhal normalizada al polinomio local
            END DO
          END DO
        END IF
C..............................................................................
C dibujamos polinomios anteriores con las derivadas igualadas en los 
C bordes de cada pixel
        DO I=1,NS
          X1=-0.5
          IF(I.EQ.1) X1=X1-1.0
          IF(I.EQ.NS) X1=X1+1.0
          Y1=2*CCC3(I)*X1+CCC2(I)
          Y1POL(I)=Y1
          X2=+0.5
          IF(I.EQ.1) X2=X2-1.0
          IF(I.EQ.NS) X2=X2+1.0
          Y2=2*CCC3(I)*X2+CCC2(I)
          Y2POL(I)=Y2
        END DO
C promediamos los valores en las fronteras de cada pixel, salvo en los bordes
C de la imagen y cuando se viole la conservacion de la derivada segunda 
C inicial (siempre y cuando esta conservacion sea requerida)
        DO I=1,NS
          PESO1(I)=0.5
          PESO2(I)=0.5
        END DO
        PESO1(1)=1.0
        PESO2(NS)=1.0
        NLOCALITER=0
C
14      LREPEAT=.FALSE.
        NLOCALITER=NLOCALITER+1
        DO I=1,NS
          IF(I.EQ.1)THEN
            Y1=Y1POL(I)
          ELSE
            Y1=Y2POL(I-1)*PESO2(I-1)+Y1POL(I)*PESO1(I)
          END IF
          IF(I.EQ.NS)THEN
            Y2=Y2POL(I)
          ELSE
            Y2=Y2POL(I)*PESO2(I)+Y1POL(I+1)*PESO1(I+1)
          END IF
          IF(I.EQ.1)THEN
            CCC1(I)=YT(I)+Y1/24.0+23.0*Y2/24.0
            CCC2(I)=1.5*Y2-0.5*Y1
            CCC3(I)=(Y2-Y1)/2.0
          ELSEIF(I.EQ.NS)THEN
            CCC1(I)=YT(I)-23.0*Y1/24.0-Y2/24.0
            CCC2(I)=1.5*Y1-0.5*Y2
            CCC3(I)=(Y2-Y1)/2.0
          ELSE
            CCC1(I)=YT(I)-(Y2-Y1)/24.0
            CCC2(I)=(Y2+Y1)/2.0
            CCC3(I)=(Y2-Y1)/2.0
          END IF
          IF(L2D)THEN                 !hay que conservar las derivadas segundas
            SIGNO0=NINT(SIGN(1.,CCC3(I)))
            IF(SIGNO0.NE.SIGNO(I))THEN
              LREPEAT=.TRUE.
              IF(I.EQ.1)THEN
                PESO1(I+1)=PESO1(I+1)/2.
                PESO2(I)=1.-PESO1(I+1)
              ELSEIF(I.EQ.NS)THEN
                PESO2(I-1)=PESO2(I-1)/2.
                PESO1(I)=1.-PESO2(I-1)
              ELSE
                PESO2(I-1)=PESO2(I-1)/2.
                PESO1(I)=1.-PESO2(I-1)
                PESO1(I+1)=PESO1(I+1)/2.
                PESO2(I)=1.-PESO1(I+1)
              END IF
            END IF
          END IF
        END DO
        IF((LREPEAT).AND.(NLOCALITER.LE.NLOCALITERMAX)) GOTO 14
C dibujamos
        IF(LPLOT.AND.LIPOL)THEN
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(14)
            DO I=1,NS
              DO L=1,30    !dividimos cada pixel en 30 subunidades para dibujar
                X=-0.5+REAL(L-1)/29.
                XP(L)=REAL(I)+X
                IF(I.EQ.1) X=X-1.0
                IF(I.EQ.NS) X=X+1.0
                POL=CCC3(I)*X*X+CCC2(I)*X+CCC1(I)
                YP(L)=POL
              END DO
              CALL PGLINE(30,XP,YP)      !senhal normalizada al polinomio local
            END DO
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END DO
        END IF
C..............................................................................
C si hay que seguir iterando, lo hacemos
        NITER=NITER+1
        IF(NITER.LT.NITERTOT)THEN
          DO I=1,NS
            X1=-0.5
            IF(I.EQ.1) X1=X1-1.0
            IF(I.EQ.NS) X1=X1+1.0
            Y1POL(I)=CCC3(I)*X1*X1+CCC2(I)*X1+CCC1(I)
            X2=+0.5
            IF(I.EQ.1) X2=X2-1.0
            IF(I.EQ.NS) X2=X2+1.0
            Y2POL(I)=CCC3(I)*X2*X2+CCC2(I)*X2+CCC1(I)
          END DO
          GOTO 10
        END IF
C..............................................................................
C si salimos y estamos dibujando, mostramos el signo de las derivadas segundas
C y el polinomio final
        IF(LPLOT)THEN
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(3)
            DO I=1,NS
              SIGNO0=NINT(SIGN(1.,CCC3(I)))
              IF((SIGNO0.NE.SIGNO(I)).AND.(LCOLOR(ITERM))) CALL PGSCI(2)
              if((real(i).gt.xmin).and.(real(i).lt.xmax))then
                if(signo0.gt.0)then
                  call pgptext(real(i),ymax-dy/15.,0.,0.5,'+')
                elseif(signo0.lt.0)then
                  call pgptext(real(i),ymax-dy/15.,0.,0.5,'-')
                else
                  call pgptext(real(i),ymax-dy/15.,0.,0.5,'0')
                end if
ccc     type*,i,signo0,signo(i),ccc3(i),ccc2(i),ccc1(i)
              end if
              IF((SIGNO0.NE.SIGNO(I)).AND.(LCOLOR(ITERM))) CALL PGSCI(3)
            END DO
            IF(LCOLOR(ITERM)) CALL PGSCI(3)
            DO I=1,NS
              DO L=1,30    !dividimos cada pixel en 30 subunidades para dibujar
                X=-0.5+REAL(L-1)/29.
                XP(L)=X+REAL(I)
                IF(I.EQ.1) X=X-1.0
                IF(I.EQ.NS) X=X+1.0
                POL=CCC3(I)*X*X+CCC2(I)*X+CCC1(I)
                YP(L)=POL
              END DO
              CALL PGLINE(30,XP,YP)      !senhal normalizada al polinomio local
            END DO
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END DO
        END IF
C
100     FORMAT(A,$)
101     FORMAT(A)
        END
