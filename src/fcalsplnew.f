C------------------------------------------------------------------------------
C Version 21-January-2000                                    file: fcalsplnew.f
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
C Program: fcalsplnew
C Classification: flux calibration
C Description: Generates a flux calibration spectrum using splines, allowing an
C absolute flux calibration.
C
Comment
C
C Ajusta curva de respuesta del detector usando splines
C 
        PROGRAM FCALSPLNEW
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER TRUELEN
        INTEGER READI
        INTEGER READILIM
        REAL READF
C
        INTEGER NPMAX           !No. maximo de puntos tabulados en fichero
        PARAMETER (NPMAX=10000)
        INTEGER NDMAX           !No. maximo de knots (ver splfit.f)
        PARAMETER (NDMAX=20)
C constantes
        REAL C                                           !speed of light (km/s)
        PARAMETER (C=2.9979246E+5)
C
        INTEGER I,J,K,L
        INTEGER NPT,NFIT,NCHANBIN
        INTEGER NBIN
        INTEGER NCBEG,NCEND,ND,ND0
        INTEGER NC0,NC1,NC2
        INTEGER NCC1,NCC2
        INTEGER NS1,NS2
        INTEGER NCOLOR
        INTEGER IFB             !funcion que calcula el numero de binning
        INTEGER NCHANXX,NSCANXX
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        REAL S0(NCMAX),S(NCMAX),X(NCMAX),T(NCMAX),TW(NCMAX),Y(NCMAX)
        REAL T2(NCMAX)
        REAL SYNTHETIC(NCMAX)
        REAL XBIN(NCMAX),YBIN(NCMAX)
        REAL XF(NCMAX),YF(NCMAX)
        REAL XSPL(NCMAX),YSPL(NCMAX),YSPL_SAVE(NCMAX)
        REAL YOVER(NCMAX)
        REAL WV(NPMAX),FV(NPMAX)
        REAL XD(NDMAX)
        REAL XC,YC,SIGMAF
        REAL YRMSTOL
        REAL XMIN,XMAX,YMIN,YMAX,DX,DY,XMINL,XMAXL
        REAL FS,FT,FCTE,FRES
        REAL YPS(NCMAX),YPT(NCMAX),XP(NCMAX)
        REAL STWVXX,DISPXX
        REAL VSIGMA,VSIGMAVEC(NCMAX),RADVEL,RCVEL,RCVEL1
        CHARACTER*1 CCONT,CH,COPC,CKMOD,CREMOV,CCHB
        CHARACTER*1 CINIT,CINITMODE,CPSYN,CRENOR
        CHARACTER*1 CFORMAT
        CHARACTER*50 CDUMMY
        CHARACTER*75 DATFILE,TABFILE,OUTFILE,OVERFILE
        LOGICAL FINDTABFILE
        LOGICAL LFIT(NCMAX),LSYNTHETIC(NCMAX)
        LOGICAL LCOLOR(MAX_ID_RED)
C
        COMMON/BLKNBIN/NBIN
        COMMON/BLKYMM/YMIN,YMAX
        COMMON/BLKFILES/DATFILE,TABFILE
        COMMON/BLKEJEL/STWV,DISP
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
C------------------------------------------------------------------------------
        THISPROGRAM='fcalsplnew'
        CALL WELCOME('21-January-2000')
C
        NCOLOR=5
        CKMOD=' '
        CREMOV='n'
        YRMSTOL=1E-6
        CALL PIDEGTER(NTERM,IDN,LCOLOR)
C
        DO J=1,NCMAX
          S(J)=0.
          X(J)=REAL(J)
        END DO
C
C cargamos espectro
        WRITE(*,100) 'Input file name'
        DATFILE=INFILEX(14,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        IF(TRUELEN(OBJECT).GT.0)THEN
          DATFILE=DATFILE(1:TRUELEN(DATFILE))//' ['//
     +     OBJECT(1:TRUELEN(OBJECT))//']'
        END IF
        IF(STWV.EQ.0.)THEN
          WRITE(*,101) 'FATAL ERROR: STWV = 0.'
          CLOSE(14)
          STOP
        END IF
        IF(DISP.EQ.0.)THEN
          WRITE(*,101) 'FATAL ERROR: DISP = 0.'
          CLOSE(14)
          STOP
        END IF
        IF(NSCAN.GT.1)THEN
          WRITE(*,101) 'WARNING: this file contains more than 1 '//
     +     'spectrum.'
          WRITE(*,100) 'Do you want to continue (y/n) '
          CCONT(1:1)=READC('y','yn')
          IF(CCONT.EQ.'y')THEN
10          WRITE(*,100) 'Spectra range (1st & last scan)'
            CALL READ2I('@',NS1,NS2)
            IF((NS1.LT.1).OR.(NS2.GT.NSCAN).OR.(NS1.GT.NS2))THEN
              WRITE(*,101)'ERROR: numbers out of range. Try again.'
              GOTO 10
            END IF
          ELSE
            CALL PGEND
            STOP
          END IF
        ELSE
          NS1=1
          NS2=1
        END IF
        WRITE(*,100) 'Reading file...'
        IF(NS1.GT.1)THEN                !skip first NS1-1 scans
          DO I=1,NS1-1
            READ(14) (S0(J),J=1,NCHAN)
          END DO
        END IF
        DO I=NS1,NS2
          READ(14) (S0(J),J=1,NCHAN)
          DO J=1,NCHAN
            S(J)=S(J)+S0(J)
          END DO
        END DO
        IF(NS2.GT.NS1)THEN
          DO J=1,NCHAN
            S(J)=S(J)/REAL(NS2-NS1+1)
          END DO
        END IF
        CLOSE(14)
        WRITE(*,101) '  ...OK!'
C
        WRITE(*,100) 'Radial velocity of this spectrum (km/s) '
        RADVEL=READF('0.0')
        RCVEL=RADVEL/C
        RCVEL1=1.+RCVEL
        RCVEL1=RCVEL1/SQRT(1.-RCVEL*RCVEL)              !correccion relativista
C
        WRITE(*,100) 'Sigma (km/s) to broad this spectrum (0=NONE) '
        VSIGMA=READF('0.0')
        IF(VSIGMA.GT.0.0)THEN
          DO J=1,NCHAN
            VSIGMAVEC(J)=VSIGMA
          END DO
          CALL BROADEN(S,S0,NCHAN,STWV,DISP,VSIGMAVEC,.FALSE.)
        END IF
C
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PLOT(NCHAN,X,S,0,0,.FALSE.,0,0)
          IF(VSIGMA.GT.0.0)THEN
            CALL PGSCI(2)
            CALL PGBIN(NCHAN,X,S0,.TRUE.)
            CALL PGSCI(1)
          END IF
          CALL PGLABEL('Channel','counts/sec','file: '//DATFILE)
        END DO
C
        DO J=1,NCHAN
          S(J)=S0(J)
        END DO
C cargamos tabla con longitudes de onda y flujos
        WRITE(*,100) 'Flux table file name'
        IF(FINDTABFILE(OBJECT,TABFILE))THEN
          WRITE(*,100)' '
          TABFILE=INFILEX(12,TABFILE,0,0,.0,.0,3,.FALSE.)
        ELSE
          TABFILE=INFILEX(12,'@',0,0,.0,.0,3,.FALSE.)
        END IF
        WRITE(*,100) 'Reading file...'
        READ(12,*) NPT
        IF(NPT.GT.NPMAX)THEN
          WRITE(*,101) 'FATAL ERROR: No. of points too large in '//
     +     'file: '//TABFILE(1:TRUELEN(TABFILE))
          CLOSE(12)
          STOP
        END IF
        DO I=1,NPT
          READ(12,*) WV(I),FV(I)
        END DO
        CLOSE(12)
        WRITE(*,110) '      No. of points read: ',NPT
C
        WRITE(*,101) 'Indicate table format:'
        WRITE(*,101) '(1) Angs, erg cm^-2 s^-1 A^-1'
        WRITE(*,101) '(2) Angs, erg cm^-2 s^-1 Hz^-1'
        WRITE(*,100) 'Option (1/2) '
        CFORMAT(1:1)=READC('1','12')
        IF(CFORMAT.EQ.'2')THEN
C Si el fichero corresponde a los datos de Massey et al. 1988, ApJ, 328, 315,
C tenemos WV en angstroms y FV en magnitudes. Pasamos las magnitudes a 
C f_nu (erg cm^-2 s^-1 Hz^-1) usando la ecuación de la página 333, y a
C continuación calculamos f_lambda (erg cm^-2 s^-1 A^-1) usando 
C f_lambda = f_nu * c/ldo^2, donde c es la velocidad de la luz en Angs/s, y
C ldo es la longitud de onda en Angs. Como hemos definido C como la velocidad
C de la luz en km/s, tenemos un factor 10**13 al usar c en Angs/s.
          DO I=1,NPT
            FV(I)=10**(-0.4*(FV(I)+48.59)+13)*C/(WV(I)*WV(I))
          END DO
        END IF
C
        DO I=1,NPT
          WV(I)=ALOG10(WV(I))
          FV(I)=ALOG10(FV(I))
        END DO
C
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PLOT(NPT,WV,FV,0,0,.FALSE.,0,0)
          CALL PGLABEL('log[wavelength (\A)]',
     +     'log[flux (erg s\u-1\d cm\u-2\d \A\u-1\d)]',
     +     'file: '//TABFILE(1:TRUELEN(TABFILE)))
        END DO
C
        IF(10.**WV(1)*RCVEL1.GT.STWV)THEN
          WRITE(*,101) 'FATAL ERROR: range covered by the table in '//
     +     'the BLUE is insufficient.'
          CALL PGEND
          STOP
        END IF
        IF(10.**WV(NPT)*RCVEL1.LT.STWV+REAL(NCHAN-1)*DISP)THEN
          WRITE(*,101) 'FATAL ERROR: range covered by the table in '//
     +     'the RED is insufficient.'
          CALL PGEND
          STOP
        END IF
C pintamos el rango que nos interesa
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM)) CALL PGSCI(2)
          CALL PGMOVE(ALOG10(STWV),YMIN)
          CALL PGDRAW(ALOG10(STWV),YMAX)
          CALL PGMOVE(ALOG10(STWV+REAL(NCHAN-1)*DISP),YMIN)
          CALL PGDRAW(ALOG10(STWV+REAL(NCHAN-1)*DISP),YMAX)
          CALL PGRECT(ALOG10(STWV),ALOG10(STWV+REAL(NCHAN-1)*DISP),
     +     YMIN,YMIN+(YMAX-YMIN)/40.)
          CALL PGRECT(ALOG10(STWV),ALOG10(STWV+REAL(NCHAN-1)*DISP),
     +     YMAX,YMAX-(YMAX-YMIN)/40.)
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
        END DO
C recuperamos escala lineal en ambos ejes
        DO I=1,NPT
          WV(I)=10.**WV(I)
          FV(I)=10.**FV(I)
        END DO
C dibujamos el rango de interes en escala lineal
        WRITE(*,100) 'Press <CR> to continue...'
        READ(*,*)
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PLOT(NPT,WV,FV,
     +     NINT(STWV-5.*DISP),NINT(STWV+REAL(NCHAN+4)*DISP),.FALSE.,0,0)
          CALL PGLABEL('wavelength (\A)',
     +     'flux (erg s\u-1\d cm\u-2\d \A\u-1\d)',
     +     'file: '//TABFILE(1:TRUELEN(TABFILE)))
        END DO
C rebineamos tabla en el rango que nos interesa
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM)) CALL PGSCI(2)
          CALL PGLINE(NPT,WV,FV)
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
        END DO
        CALL REBINING(WV,FV,NPT,TW,T,NCHAN,STWV,DISP)
        CALL RVREBIN(RADVEL,NCHAN,T,T2,STWV,DISP)
        DO J=1,NCHAN
          T(J)=T2(J)
        END DO
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM)) CALL PGSCI(3)
          CALL PGLINE(NCHAN,TW,T)
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
        END DO
C------------------------------------------------------------------------------
C IMPORTANTE: el flujo en la tabla esta en erg s^-1 cm^-2 A^-1. Sin embargo,
C el espectro esta en cuentas pixel^-1. En principio habria que convertir el
C flujo por Angs en flujo por pixel. Sin embargo, al dividir posteriormente 
C un espectro observado por la curva respuesta, se vuelve a hacer este cambio
C de unidades pero en sentido contrario. Por tanto, no es necesario hacerlo
C ahora pero si conviene recordar cuales son las unidades correctas.
C------------------------------------------------------------------------------
        DO I=1,NCHAN
          XP(I)=REAL(I)
        END DO
C calculamos la media de S
        FS=0.
        DO I=1,NCHAN
          FS=FS+S(I)
        END DO
        FS=FS/REAL(NCHAN)
C calculamos la media de T
        FT=0.
        DO I=1,NCHAN
          FT=FT+T(I)
        END DO
        FT=FT/REAL(NCHAN)
C creamos los espectros a dibujar normalizados
        DO I=1,NCHAN
          YPS(I)=S(I)/FS
        END DO
        DO I=1,NCHAN
          YPT(I)=T(I)/FT
        END DO
C calculamos la curva respuesta
        DO I=1,NCHAN
          IF(YPT(I).GT.0.)THEN
            Y(I)=YPS(I)/YPT(I)
          ELSE
            Y(I)=1.0
          END IF
        END DO
        print*,'ft,fs,fs/ft:',ft,fs,fs/ft
C calculamos la media de Y
ccc     FRES=0.
ccc     DO I=1,NCHAN
ccc       FRES=FRES+Y(I)
ccc     END DO
ccc     FRES=FRES/REAL(NCHAN)
ccc     print*,'fres:',fres
C renormalizamos Y
ccc     DO I=1,NCHAN
ccc       Y(I)=Y(I)/FRES
ccc     END DO
C
21      WRITE(*,100)'Enter First and Last Channel '
        WRITE(CDUMMY,'(I1,A1,I10)')1,',',NCHAN
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        CALL READ2I(CDUMMY(1:L),NCBEG,NCEND)
        IF((NCBEG.LT.1).OR.(NCEND.GT.NCHAN).OR.
     +   (NCEND.LT.NCBEG))THEN
          WRITE(*,101)'ERROR: numbers out of range. Try again.'
          GOTO 21
        END IF
C calculamos limites de los plots
        XMIN=REAL(NCBEG)
        XMAX=REAL(NCEND)
        YMIN=1.E30
        YMAX=-1.E30
        DO I=NCBEG,NCEND
          IF(YPS(I).LT.YMIN) YMIN=YPS(I)
          IF(YPT(I).LT.YMIN) YMIN=YPT(I)
          IF(Y(I).LT.YMIN) YMIN=Y(I)
          IF(YPS(I).GT.YMAX) YMAX=YPS(I)
          IF(YPT(I).GT.YMAX) YMAX=YPT(I)
          IF(Y(I).GT.YMAX) YMAX=Y(I)
        END DO
        DX=XMAX-XMIN
        XMIN=XMIN-DX/50.
        XMAX=XMAX+DX/50.
        DY=YMAX-YMIN
        YMIN=YMIN-DY/20.
        YMAX=YMAX+DY/20.
        XMINL=(XMIN-1.)*DISP+STWV
        XMAXL=(XMAX-1.)*DISP+STWV
C
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,-2)
          CALL PGIDEN_RED
          CALL PGSWIN(XMINL,XMAXL,YMIN,YMAX)
          CALL PGBOX('CMTS',0.,0,'BCNST',0.,0)
          CALL PGSWIN(XMIN,XMAX,YMIN,YMAX)
          CALL PGBOX('BNTS',0.,0,' ',0.,0)
          IF(LCOLOR(ITERM)) CALL PGSCI(3)
          CALL PGSCH(0.7)
          CALL PGBIN(NCHAN,XP,YPS,.TRUE.)
          WRITE(CDUMMY,*) FS
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          CALL PGMTEXT('T',3.0,0.,0.,
     +     DATFILE(1:TRUELEN(DATFILE))//'/'//CDUMMY(1:L))
          CALL PGMTEXT('T',4.5,0.,0.,
     +     '(count s\u-1\d pixel\u-1\d)')
          IF(LCOLOR(ITERM)) CALL PGSCI(4)
          CALL PGBIN(NCHAN,XP,YPT,.TRUE.)
          WRITE(CDUMMY,*) FT
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          CALL PGMTEXT('T',3.0,1.,1.,
     +     TABFILE(1:TRUELEN(TABFILE))//'/'//CDUMMY(1:L))
          CALL PGMTEXT('T',4.5,1.,1.,
     +     '(erg s\u-1\d cm\u-2\d \A\u-1\d)')
          IF(LCOLOR(ITERM)) CALL PGSCI(5)
          WRITE(CDUMMY,*) FS/FT
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          CALL PGMTEXT('T',3.0,.5,.5,
     +     'Response curve'//'/'//CDUMMY(1:L))
          CALL PGMTEXT('T',4.5,.5,.5,
     +     '(count erg\u-1\d cm\u2\d \A pixel\u-1\d)')
          CALL PGBIN(NCHAN,XP,Y,.TRUE.)
          CALL PGSCH(1.0)
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
          CALL PGLABEL('Channel','Relative Flux',CHAR(32))
        END DO
C------------------------------------------------------------------------------
22      WRITE(*,100)'Bin width '
        NBIN=READILIM('@',1,NCHAN)
C si se ha solicitado, hacemos el binning
        IF(NBIN.NE.1)THEN
          I=0
          J=0
          DO WHILE (I.LT.NCHAN)
            J=J+1
            K=0
            XBIN(J)=0.
            YBIN(J)=0.
            DO WHILE ((K.LT.NBIN).AND.(I.LT.NCHAN))
              I=I+1
              K=K+1
              XBIN(J)=XBIN(J)+X(I)
              YBIN(J)=YBIN(J)+Y(I)
            END DO
            XBIN(J)=XBIN(J)/REAL(K)
            YBIN(J)=YBIN(J)/REAL(K)
          END DO
          NCHANBIN=J
        ELSE
          NCHANBIN=NCHAN
          DO I=1,NCHAN
            XBIN(I)=X(I)
            YBIN(I)=Y(I)
          END DO
        END IF
C redibujamos con el nuevo binning
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,-2)
          CALL PGIDEN_RED
          CALL PGSWIN(XMINL,XMAXL,YMIN,YMAX)
          CALL PGBOX('CMTS',0.,0,'BCNST',0.,0)
          CALL PGSWIN(XMIN,XMAX,YMIN,YMAX)
          CALL PGBOX('BNTS',0.,0,' ',0.,0)
          IF(LCOLOR(ITERM)) CALL PGSCI(3)
ccc       CALL PGSCH(0.7)
ccc       CALL PGBIN(NCHAN,XP,YPS,.TRUE.)
ccc       CALL PGMTEXT('T',3.0,0.,0.,DATFILE)
ccc       IF(LCOLOR(ITERM)) CALL PGSCI(4)
ccc       CALL PGBIN(NCHAN,XP,YPT,.TRUE.)
ccc       CALL PGMTEXT('T',3.0,1.,1.,TABFILE)
ccc       IF(LCOLOR(ITERM)) CALL PGSCI(5)
ccc       CALL PGMTEXT('T',3.0,.5,.5,'Response curve')
ccc       CALL PGBIN(NCHAN,XP,Y,.TRUE.)
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
          CALL PGBIN(NCHANBIN,XBIN,YBIN,.TRUE.)
          CALL PGSCH(1.0)
          IF(NBIN.NE.1)THEN
            WRITE(CDUMMY,*)NBIN
            CALL RMBLANK(CDUMMY,CDUMMY,L)
            CDUMMY='Channels (binning='//CDUMMY(1:L)//')'
          ELSE
            CDUMMY='Channels'
          END IF
          CALL PGLABEL(CDUMMY,'Relative Flux',CHAR(32))
        END DO
C
        WRITE(*,100)'Change binning (y/n) '
        CCHB(1:1)=READC('n','yn')
        IF(CCHB.EQ.'y')GOTO 22
C
        DO I=1,NCHANBIN
          LFIT(I)=.FALSE.
        END DO
C
        DO I=IFB(NCBEG),IFB(NCEND)
          LFIT(I)=.TRUE.
        END DO
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM)) CALL PGSCI(2)
          CALL PLOT(NCHANBIN,XBIN,YBIN,IFB(NCBEG),IFB(NCEND),.TRUE.,0,0)
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
        END DO
C
ccc35      WRITE(*,*)
        WRITE(*,*)
        WRITE(*,100)'Remove regions with [m]ouse, [k]eyboard, [n]one '//
     +   ' (m/k/n) '
        CREMOV(1:1)=READC(CREMOV,'mkn')
        IF(CREMOV.EQ.'n') GOTO 40
36      IF(CREMOV.EQ.'m')THEN
          WRITE(*,101)'Remove region with mouse (q=X=EXIT):'
          WRITE(*,100)'Press mouse (point #1)...'
          CALL PGBAND(6,0,0.,0.,XC,YC,CH)
          IF((CH.EQ.'q').OR.(CH.EQ.'Q').OR.(CH.EQ.'X').OR.
     +     (CH.EQ.'x'))THEN
            WRITE(*,*)
            GOTO 40
          END IF
          NC1=NINT(XC)
          WRITE(*,110)'     position: ',NC1
          WRITE(*,100)'Press mouse (point #2)...'
          CALL PGBAND(4,0,REAL(NC1),0.,XC,YC,CH)
          IF((CH.EQ.'q').OR.(CH.EQ.'Q').OR.(CH.EQ.'X').OR.
     +     (CH.EQ.'x'))THEN
            WRITE(*,*)
            GOTO 40
          END IF
          NC2=NINT(XC)
          WRITE(*,110)'     position: ',NC2
        ELSE
          WRITE(*,100)'Channel region (0,0=EXIT) '
          CALL READ2I('0,0',NC1,NC2)
          IF((NC1.EQ.0).AND.(NC2.EQ.0)) GOTO 40
        END IF
        IF(NC1.GT.NC2)THEN
          NC0=NC1
          NC1=NC2
          NC2=NC0
        END IF
        IF(NC1.LT.1) NC1=1
        IF(NC2.GT.NCHAN) NC2=NCHAN
        DO I=IFB(NC1),IFB(NC2)
          LFIT(I)=.FALSE.
        END DO
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGSLS(2)
          CALL PGMOVE(REAL(NC1),YMIN)
          CALL PGDRAW(REAL(NC1),YMAX)
          CALL PGMOVE(REAL(NC2),YMIN)
          CALL PGDRAW(REAL(NC2),YMAX)
          CALL PGSLS(1)
          CALL PGRECT(REAL(NC1),REAL(NC2),YMIN,YMIN+(YMAX-YMIN)/40.)
          CALL PGRECT(REAL(NC1),REAL(NC2),YMAX,YMAX-(YMAX-YMIN)/40.)
          CALL PLOT(NCHANBIN,XBIN,YBIN,IFB(NC1),IFB(NC2),.TRUE.,0,0)
        END DO
        GOTO 36
C
40      NFIT=0
        DO I=1,NCHANBIN
          IF(LFIT(I))THEN
            NFIT=NFIT+1
            XF(NFIT)=XBIN(I)
            YF(NFIT)=YBIN(I)
          END IF
        END DO
C
        DO I=1,NCHANBIN
          IF(LFIT(I)) GOTO 41
        END DO
41      NC1=NINT(XBIN(I)-REAL(NBIN)/2.+0.5)  !limite inferior del binning
        DO I=NCHANBIN,1,-1
          IF(LFIT(I)) GOTO 42
        END DO
42      NC2=NINT(XBIN(I)+REAL(NBIN)/2.-0.5)  !limite superior del binning
        IF(NC2.GT.NCHAN) NC2=NCHAN           !puede ser mayor que NCHAN
        IF(NCBEG.LT.NC1) NCBEG=NC1
        IF(NCEND.GT.NC2) NCEND=NC2
C
        ND=2
        XD(1)=REAL(NCBEG)
        XD(2)=REAL(NCEND)
45      WRITE(*,*)
        WRITE(*,110)'Current number of Knots: ',ND
        DO I=1,ND
          WRITE(*,'(A,I2,A,I6)')'Knot #',I,'  located at: ',
     +     NINT(XD(I))
        END DO
        WRITE(CDUMMY,'(I10,A1,I10)')NCBEG,',',NCEND
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,100)'Enter new Knots by mouse/keyboard/automatic'//
     +   ' (m/k/a) '
        IF(CKMOD.EQ.' ')THEN
          CKMOD(1:1)=READC('a','mka')
        ELSE
          CKMOD(1:1)=READC(CKMOD,'mka')
        END IF
        IF(CKMOD.EQ.'m')THEN
          ND0=0
          WRITE(*,101)'Press mouse between '//CDUMMY(1:L)//
     +     ' (excluding '//CDUMMY(1:L)//')'
46        WRITE(*,'(A,I2,A,$)')'New Knot #',ND0+1,
     +     '   Press mouse (q=X=EXIT)...'
          CALL PGBAND(6,0,0.,0.,XC,YC,CH)
          WRITE(*,110)'    position: ',NINT(XC)
          IF((CH.EQ.'q').OR.(CH.EQ.'Q').OR.(CH.EQ.'X').OR.
     +     (CH.EQ.'x'))THEN
            GOTO 47
          END IF
          XD(ND+ND0+1)=ANINT(XC)
          IF((NINT(XD(ND+ND0+1)).LE.NCBEG).OR.
     +     (NINT(XD(ND+ND0+1)).GE.NCEND))THEN
            WRITE(*,101)'ERROR: knot position out of range. Try again.'
            GOTO 46
          END IF
          ND0=ND0+1
          IF(ND+ND0.EQ.18)THEN
            WRITE(*,101)'WARNING: Maximum number of Knots is 18.'
            WRITE(*,101)'No more Knots allowed.'
            GOTO 47
          END IF
          GOTO 46
47        CONTINUE
        ELSE
48        WRITE(*,100)'No. of new Knots between '//CDUMMY(1:L)//
     +       ' (excluding '//CDUMMY(1:L)//')'
          ND0=READI('@')
          IF(ND+ND0.GT.18)THEN
            WRITE(*,101)'ERROR: Maximum number of knots must be < 18.'//
     +       ' Try again.'
            GOTO 48
          END IF
          IF(CKMOD.EQ.'a')THEN
            DO I=ND+1,ND+ND0
              XD(I)=REAL(NCBEG)+REAL(I-ND)*REAL(NCEND-NCBEG)/REAL(ND0+1)
            END DO
          ELSE
            DO I=ND+1,ND+ND0
49            WRITE(*,'(A,I2.2,$)')'New Knot #',I-ND
              XD(I)=READF('@')
              IF((NINT(XD(I)).LE.NCBEG).OR.(NINT(XD(I)).GE.NCEND))THEN
                WRITE(*,101)'ERROR: knot position out of range. '//
     +           'Try again.'
                GOTO 49
              END IF
              DO L=1,I-1
                IF(XD(L).EQ.XD(I))THEN
                  WRITE(*,101)'ERROR: this Knot already exist. '//
     +             'Try again.'
                  GOTO 49
                END IF
              END DO
            END DO
          END IF
        END IF
        ND=ND+ND0
        WRITE(*,110)'Total number of Knots: ',ND
        CALL ORDENA1F(ND,XD)
        IF(CKMOD.EQ.'a')THEN
          DO I=1,ND
            WRITE(CDUMMY,*)XD(I)
            CALL RMBLANK(CDUMMY,CDUMMY,L)
            WRITE(*,'(A,I2,A)')'Knot #',I,' located at '//CDUMMY(1:L)
          END DO
        END IF
50      WRITE(CDUMMY,*)YRMSTOL
        WRITE(*,100)'YRMSTOL for DOWNHILL '
        YRMSTOL=READF(CDUMMY)
        NCOLOR=NCOLOR+1
        IF(NCOLOR.EQ.15) NCOLOR=6
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM)) CALL PGSCI(NCOLOR)
        END DO
        CALL SPLFIT(NFIT,XF,YF,ND,XD,YRMSTOL,NCHAN,XSPL,YSPL,
     +   1.,REAL(NCHAN),SIGMAF,.TRUE.)
        WRITE(CDUMMY,*)SIGMAF
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,101)'Sigma of the fit: '//CDUMMY(1:L)
        WRITE(*,*)
        WRITE(*,101)'>>> Final Knot locations:'
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
          DO I=1,ND
            CALL PGMOVE(XD(I),YSPL(NINT(XD(I)))+(YMAX-YMIN)/30.)
            CALL PGDRAW(XD(I),YSPL(NINT(XD(I)))-(YMAX-YMIN)/30.)
          END DO
        END DO
        DO I=1,ND
          WRITE(CDUMMY,*)XD(I)
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          WRITE(*,'(A,I2,A)')'Knot #',I,'  ---> channel: '//CDUMMY(1:L)
        END DO
        WRITE(*,*)
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          DO I=1,ND
            WRITE(CDUMMY,*)I
            CALL RMBLANK(CDUMMY,CDUMMY,L)
            CALL PGPTEXT(XD(I),YSPL(NINT(XD(I)))+(YMAX-YMIN)/24.,0.,0.5,
     +       CDUMMY(1:L))
          END DO
        END DO
55      WRITE(*,*)
        WRITE(*,101) '(1) Repeat whole fit'
        WRITE(*,101) '(2) Add Knots'
        IF(ND.GT.2)THEN
          WRITE(*,101) '(3) List/Remove Knots'
        END IF
        WRITE(*,101) '(f) Fit'
        WRITE(*,101) '(r) Replot'
        WRITE(*,101) '(o) Overplot curve from other file'
        WRITE(*,101) '(g) Cut & paste --> synthetic spectrum'
        WRITE(*,101) '(s) Save normalized & binned response curve '//
     +   'before fit'
        WRITE(*,101) '(w) Save normalized tabulated spectrum'
        WRITE(*,101) '(z) Save synthetic spectrum'
        WRITE(*,101) '(x) Save fit'
        WRITE(*,101) '(0) QUIT'
        WRITE(*,100) 'Option '
        COPC(1:1)=READC('r','x0123frogswz')
        CALL CHLOWER(COPC)
C..............................................................................
        IF(COPC.EQ.'x')THEN
          WRITE(*,100)'Renormalize fitted spectrum (y/n) '
          CRENOR(1:1)=READC('y','yn')
          IF(CRENOR.EQ.'y')THEN
            WRITE(CDUMMY,'(A,I10)')'1,',NCHAN
            CALL RMBLANK(CDUMMY,CDUMMY,L)
            WRITE(*,100)'Channel region to obtain normalization factor '
            CALL READ2I(CDUMMY(1:L),NC1,NC2)
            FRES=0.
            DO J=NC1,NC2
              FRES=FRES+YSPL(J)
            END DO
            FRES=FRES/REAL(NC2-NC1+1)
          ELSE
            FRES=1.0
          END IF
          WRITE(*,100)'>>> Normalization factor: '
          WRITE(*,*) FRES
          DO J=1,NCHAN
            YSPL_SAVE(J)=YSPL(J)/FRES
          END DO
          WRITE(*,100)'Output file name'
          OUTFILE=OUTFILEX(15,'@',1,NCHAN,STWV,DISP,1,.FALSE.)
          WRITE(15) (YSPL_SAVE(J),J=1,NCHAN)
          CLOSE(15)
          GOTO 55
C..............................................................................
        ELSEIF(COPC.EQ.'0')THEN
          GOTO 80
C..............................................................................
        ELSEIF(COPC.EQ.'1')THEN
          GOTO 21
C..............................................................................
        ELSEIF(COPC.EQ.'2')THEN
          GOTO 45
C..............................................................................
        ELSEIF(COPC.EQ.'3')THEN
60        DO I=1,ND
            WRITE(CDUMMY,*)XD(I)
            CALL RMBLANK(CDUMMY,CDUMMY,L)
            WRITE(*,'(A,I2,A)')'Knot #',I,' located at '//CDUMMY(1:L)
          END DO
62        IF(ND.EQ.2)THEN
            WRITE(*,101)'No. of Knots = minimum value.'
            GOTO 55
          END IF
          WRITE(*,100)'Which Knot number are you removing (0=EXIT) '
          ND0=READILIM('0',0,ND)
          IF(ND0.EQ.0) GOTO 55
          IF((ND0.EQ.1).OR.(ND0.EQ.ND))THEN
            WRITE(*,101)'WARNING: 1st and Last Knot cannot be removed.'
            GOTO 62
          END IF
          IF(ND0.LT.ND)THEN
            DO I=ND0,ND-1
              XD(I)=XD(I+1)
            END DO
          END IF
          ND=ND-1
          GOTO 60
C..............................................................................
        ELSEIF(COPC.EQ.'f')THEN
          GOTO 50
C..............................................................................
        ELSEIF(COPC.EQ.'r')THEN
          XMIN=1.
          XMAX=REAL(NCHAN)
          YMIN=1.E30
          YMAX=-1.E30
          DO I=1,NCHAN
            IF(YPS(I).LT.YMIN) YMIN=YPS(I)
            IF(YPT(I).LT.YMIN) YMIN=YPT(I)
            IF(Y(I).LT.YMIN) YMIN=Y(I)
            IF(YPS(I).GT.YMAX) YMAX=YPS(I)
            IF(YPT(I).GT.YMAX) YMAX=YPT(I)
            IF(Y(I).GT.YMAX) YMAX=Y(I)
          END DO
          DX=XMAX-XMIN
          XMIN=XMIN-DX/50.
          XMAX=XMAX+DX/50.
          DY=YMAX-YMIN
          YMIN=YMIN-DY/20.
          YMAX=YMAX+DY/20.
          WRITE(*,100) 'YMIN'
          WRITE(CDUMMY,*) YMIN
          YMIN=READF(CDUMMY)
          WRITE(*,100) 'YMAX'
          WRITE(CDUMMY,*) YMAX
          YMAX=READF(CDUMMY)
          XMINL=(XMIN-1.)*DISP+STWV
          XMAXL=(XMAX-1.)*DISP+STWV
C
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,-2)
            CALL PGIDEN_RED
            CALL PGSWIN(XMINL,XMAXL,YMIN,YMAX)
            CALL PGBOX('CMTS',0.,0,'BCNST',0.,0)
            CALL PGSWIN(XMIN,XMAX,YMIN,YMAX)
            CALL PGBOX('BNTS',0.,0,' ',0.,0)
            CALL PGSCH(0.7)
            IF(LCOLOR(ITERM)) CALL PGSCI(3)
            CALL PGBIN(NCHAN,XP,YPS,.TRUE.)
            WRITE(CDUMMY,*) FS
            CALL RMBLANK(CDUMMY,CDUMMY,L)
            CALL PGMTEXT('T',3.0,0.,0.,
     +       DATFILE(1:TRUELEN(DATFILE))//'/'//CDUMMY(1:L))
            CALL PGMTEXT('T',4.5,0.,0.,
     +       '(count s\u-1\d pixel\u-1\d)')
            IF(LCOLOR(ITERM)) CALL PGSCI(4)
            CALL PGBIN(NCHAN,XP,YPT,.TRUE.)
            WRITE(CDUMMY,*) FT
            CALL RMBLANK(CDUMMY,CDUMMY,L)
            CALL PGMTEXT('T',3.0,1.,1.,
     +       TABFILE(1:TRUELEN(TABFILE))//'/'//CDUMMY(1:L))
            CALL PGMTEXT('T',4.5,1.,1.,
     +       '(erg s\u-1\d cm\u-2\d \A\u-1\d)')
            IF(LCOLOR(ITERM)) CALL PGSCI(5)
            WRITE(CDUMMY,*) FS/FT
            CALL RMBLANK(CDUMMY,CDUMMY,L)
            CALL PGMTEXT('T',3.0,.5,.5,
     +       'Response curve'//'/'//CDUMMY(1:L))
            CALL PGMTEXT('T',4.5,.5,.5,
     +       '(count erg\u-1\d cm\u2\d \A pixel\u-1\d)')
            CALL PGBIN(NCHAN,XP,Y,.TRUE.)
            IF(LCOLOR(ITERM)) CALL PGSCI(NCOLOR)
            CALL PGBIN(NCHANBIN,XBIN,YBIN,.TRUE.)
            IF(NBIN.NE.1)THEN
              WRITE(CDUMMY,*)NBIN
              CALL RMBLANK(CDUMMY,CDUMMY,L)
              CDUMMY='Channels (binning='//CDUMMY(1:L)//')'
            ELSE
              CDUMMY='Channels'
            END IF
            CALL PGSCH(1.0)
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
            CALL PGLABEL(CDUMMY,'Relative Flux',CHAR(32))
            DO I=1,ND
              CALL PGMOVE(XD(I),YSPL(NINT(XD(I)))+(YMAX-YMIN)/30.)
              CALL PGDRAW(XD(I),YSPL(NINT(XD(I)))-(YMAX-YMIN)/30.)
            END DO
            DO I=1,ND
              WRITE(CDUMMY,*)I
              CALL RMBLANK(CDUMMY,CDUMMY,L)
              CALL PGPTEXT(XD(I),YSPL(NINT(XD(I)))+(YMAX-YMIN)/24.,
     +         0.,0.5,CDUMMY(1:L))
            END DO
            CALL PGBIN(NCHAN,XSPL,YSPL,.TRUE.)
          END DO
          GOTO 55
C..............................................................................
        ELSEIF(COPC.EQ.'o')THEN
          WRITE(*,100)'(Overplot) file name'
          OVERFILE=INFILEX(40,'@',NSCANXX,NCHANXX,STWVXX,DISPXX,1,
     +     .FALSE.)
          IF(NSCANXX.GT.1)THEN
            WRITE(*,101)'WARNING: this file contains more than'//
     +       ' 1 spectrum.'
          END IF
          IF(NCHANXX.NE.NCHAN)THEN
            WRITE(*,101)'ERROR: number of channels is different.'
            GOTO 55
          END IF
          IF(STWVXX.NE.STWV)THEN
            WRITE(*,101)'ERROR: STWV is different.'
            GOTO 55
          END IF
          IF(DISPXX.NE.DISP)THEN
            WRITE(*,101)'ERROR: DISP is different.'
            GOTO 55
          END IF
          READ(40) (YOVER(J),J=1,NCHANXX)
          CLOSE(40)
          NCOLOR=NCOLOR+1
          IF(NCOLOR.EQ.15) NCOLOR=6
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(NCOLOR)
            CALL PLOT(NCHAN,XSPL,YOVER,0,0,.TRUE.,0,0)
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END DO
          GOTO 55
C..............................................................................
        ELSEIF(COPC.EQ.'g')THEN
          WRITE(*,*)
          WRITE(*,101)'* Synthetic spectrum: '
          WRITE(*,100)'Initialize synthetic spectrum (y/n) '
          CINIT(1:1)=READC('n','yn')
          IF(CINIT.EQ.'y')THEN
            WRITE(*,101)'1 - Cte'
            WRITE(*,101)'2 - Normalized & binned response curve '//
     +       'before fit'
            WRITE(*,100)'Option (1/2) '
            CINITMODE(1:1)=READC('2','12')
            IF(CINITMODE.EQ.'1')THEN
              WRITE(*,100)'Cte '
              FCTE=READF('1.0')
              WRITE(*,100)'Initializing...'
              DO J=1,NCHAN
                SYNTHETIC(J)=FCTE
              END DO
              WRITE(*,101)'...OK!'
            ELSE
              WRITE(*,100)'Initializing...'
              IF(NBIN.EQ.1)THEN
                DO J=1,NCHAN
                  SYNTHETIC(J)=YBIN(J)
                END DO
              ELSE
                I=0
                J=0
                DO WHILE (I.LT.NCHAN)
                  J=J+1
                  K=0
                  DO WHILE ((K.LT.NBIN).AND.(I.LT.NCHAN))
                    I=I+1
                    K=K+1
                    SYNTHETIC(I)=YBIN(J)
                  END DO
                END DO
              END IF
              WRITE(*,101)'...OK!'
            END IF
          END IF
C
          DO J=1,NCHAN
            LSYNTHETIC(J)=.FALSE.
          END DO
C
          WRITE(*,101)'* Define regions from last fit to be stored '//
     +     'in synthetic spectrum:'
          NCC1=1
          NCC2=1
          DO WHILE((NCC1.NE.0).AND.(NCC2.NE.0))
            WRITE(*,100)'Channels (0,0=EXIT) '
            CALL READ2I('0,0',NCC1,NCC2)
            IF((NCC1.EQ.0).AND.(NCC2.EQ.0))THEN
            ELSEIF((NCC1.LT.1).OR.(NCC2.GT.NCHAN).OR.(NCC1.GT.NCC2))THEN
              WRITE(*,101)'ERROR: invalid entry. Try again.'
            ELSE
              DO J=NCC1,NCC2
                LSYNTHETIC(J)=.TRUE.
              END DO
            END IF
          END DO
          DO J=1,NCHAN
            IF(LSYNTHETIC(J)) SYNTHETIC(J)=YSPL(J)
          END DO
C
          WRITE(*,100)'Plot synthetic spectrum (y/n) '
          CPSYN(1:1)=READC('y','yn')
          IF(CPSYN.EQ.'y')THEN
            NCOLOR=NCOLOR+1
            IF(NCOLOR.EQ.15) NCOLOR=6
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              IF(LCOLOR(ITERM)) CALL PGSCI(NCOLOR)
              CALL PGBIN(NCHAN,XP,SYNTHETIC,.TRUE.)
              IF(LCOLOR(ITERM)) CALL PGSCI(1)
            END DO
          END IF
C
          GOTO 55
C..............................................................................
        ELSEIF(COPC.EQ.'s')THEN
          WRITE(*,100)'Output file name'
          OUTFILE=OUTFILEX(50,'@',1,NCHAN,STWV,DISP,1,.FALSE.)
          IF(NBIN.EQ.1)THEN
            DO J=1,NCHAN
              YOVER(J)=YBIN(J)*FS/FT
            END DO
            WRITE(50) (YOVER(J),J=1,NCHAN)
          ELSE
            I=0
            J=0
            DO WHILE (I.LT.NCHAN)
              J=J+1
              K=0
              DO WHILE ((K.LT.NBIN).AND.(I.LT.NCHAN))
                I=I+1
                K=K+1
                YOVER(I)=YBIN(J)*FS/FT
              END DO
            END DO
            WRITE(50) (YOVER(J),J=1,NCHAN)
          END IF
          CLOSE(50)
          GOTO 55
C..............................................................................
        ELSEIF(COPC.EQ.'z')THEN
          WRITE(*,100)'Output file name'
          OUTFILE=OUTFILEX(15,'@',1,NCHAN,STWV,DISP,1,.FALSE.)
          WRITE(15) (SYNTHETIC(J),J=1,NCHAN)
          CLOSE(15)
          GOTO 55
C..............................................................................
        ELSEIF(COPC.EQ.'w')THEN
          WRITE(*,100)'Output file name'
          OUTFILE=OUTFILEX(50,'@',1,NCHAN,STWV,DISP,1,.FALSE.)
          WRITE(50) (T(J),J=1,NCHAN)
          GOTO 55
C..............................................................................
        ELSE
          WRITE(*,101)'ERROR: option out of range. Try again.'
          GOTO 55
        END IF
C------------------------------------------------------------------------------
80      CALL PGEND
C------------------------------------------------------------------------------
        STOP
100     FORMAT(A,$)
101     FORMAT(A)
110     FORMAT(A,I6)
        END
C
C******************************************************************************
C Calcula el numero de binning al que corresponde el canal N. Para ello
C tiene en cuenta que el binning se ha realizado cada NBIN canales
        INTEGER FUNCTION IFB(N)
        IMPLICIT NONE
        INTEGER N
        INTEGER NBIN
        COMMON/BLKNBIN/NBIN
C------------------------------------------------------------------------------
        IFB=N/NBIN+1
        IF(MOD(N,NBIN).EQ.0) IFB=IFB-1
        END
C
C******************************************************************************
C Si LOVER=.TRUE. y MODE=1 dibujamos el eje X superior en Angstroms
C Si IYMIN=1 se fuerza a que YMIN=0 (siempre con LOVER=.FALSE.)
C Si IYMIN=2 se invierte el eje Y
        SUBROUTINE PLOT(N,X,Y,IXMIN,IXMAX,LOVER,MODE,IYMIN)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INTEGER N
        INTEGER IXMIN,IXMAX
        INTEGER I
        REAL X(N),Y(N)
        REAL XP(NCMAX),YP(NCMAX)
        REAL XMIN,XMAX,YMIN,YMAX,DX,DY,YYYY
        REAL XMINL,XMAXL
        LOGICAL LOVER
        INTEGER MODE,IYMIN
C
        COMMON/BLKYMM/YMIN,YMAX
        COMMON/BLKEJEL/STWV,DISP
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        IF(.NOT.LOVER)THEN
          IF((IXMIN.EQ.0).AND.(IXMAX.EQ.0))THEN
            CALL FINDMM(N,X,XMIN,XMAX)
            DX=XMAX-XMIN
            XMIN=XMIN-DX/50.
            XMAX=XMAX+DX/50.
            CALL FINDMM(N,Y,YMIN,YMAX)
            DY=YMAX-YMIN
            YMIN=YMIN-DY/10.
            YMAX=YMAX+DY/10.
          ELSE
            XMIN=REAL(IXMIN)
            XMAX=REAL(IXMAX)
            YMIN=1.E30
            YMAX=-YMIN
            DO I=1,N
              IF((X(I).GE.XMIN).AND.(X(I).LE.XMAX))THEN
                IF(Y(I).LT.YMIN) YMIN=Y(I)
                IF(Y(I).GT.YMAX) YMAX=Y(I)
              END IF
            END DO
            DY=YMAX-YMIN
            YMIN=YMIN-DY/10.
            YMAX=YMAX+DY/10.
          END IF
          IF(IYMIN.EQ.1) YMIN=0.
          IF(IYMIN.EQ.2)THEN
            YYYY=YMIN
            YMIN=YMAX
            YMAX=YYYY
          END IF
          IF(MODE.EQ.1)THEN
            CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,-2)
            CALL PGIDEN_RED
            XMINL=(XMIN-1.)*DISP+STWV
            XMAXL=(XMAX-1.)*DISP+STWV
            CALL PGSWIN(XMINL,XMAXL,YMIN,YMAX)
            CALL PGBOX('CMTS',0.,0,'BCNST',0.,0)
            CALL PGSWIN(XMIN,XMAX,YMIN,YMAX)
            CALL PGBOX('BNTS',0.,0,' ',0.,0)
          ELSE
            CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,0)
            CALL PGIDEN_RED
          END IF
          CALL PGBIN(N,X,Y,.TRUE.)
        ELSE
          IF((IXMIN.EQ.0).AND.(IXMAX.EQ.0))THEN
            CALL PGBIN(N,X,Y,.TRUE.)
          ELSE
            DO I=IXMIN,IXMAX
              XP(I-IXMIN+1)=X(I)
              YP(I-IXMIN+1)=Y(I)
            END DO
            CALL PGBIN(IXMAX-IXMIN+1,XP,YP,.TRUE.)
          END IF
        END IF
C
        END
C
C******************************************************************************
C Buscar, a partir del nombre del objeto OBJECT, el nombre del fichero
C TABFILE que puede contener la tabla con el flujo de la estrella. Si la
C subrutina no encuentra ningun fichero apropiado retorna .FALSE.
        LOGICAL FUNCTION FINDTABFILE(OBJECT,TABFILE)
        IMPLICIT NONE
        INTEGER TRUELEN
C
        CHARACTER*(*) OBJECT
        CHARACTER*(*) TABFILE
C
        INTEGER LMIN,LRED
        INTEGER K,I
        INTEGER LOBJ,LTAB
        CHARACTER*255 LOCALFILE,REDUCEMEDIR
        LOGICAL LOGFILE
C------------------------------------------------------------------------------
        FINDTABFILE=.FALSE.                !salvo que se demuestre lo contrario
C
        LOBJ=TRUELEN(OBJECT)
        TABFILE=OBJECT(1:LOBJ)
        CALL RMBLANK(TABFILE,TABFILE,LTAB)
        CALL CHLOWER(TABFILE(1:LTAB))
C
        CALL GETENV('reduceme_dir',REDUCEMEDIR)
        LRED=TRUELEN(REDUCEMEDIR)
C
        K=0                                 !eliminamos caracteres no esperados
        DO I=1,LTAB
          IF(INDEX('abcdefghijklmnopqrstuvwxyz1234567890',
     +     TABFILE(I:I)).EQ.0)THEN                !hay que eliminar el caracter
            IF(I.LT.LTAB)THEN
              TABFILE(I:LTAB-1)=TABFILE(I+1:LTAB)
            END IF
            K=K+1
          END IF
        END DO
        LTAB=LTAB-K
C
        LMIN=5                                      !exigimos un taman~o minimo
        IF(LTAB.LT.LMIN) RETURN
C
        K=LMIN
        DO WHILE(K.LE.LTAB)
          DO I=1,LTAB-K+1
            LOCALFILE(1:LRED)=REDUCEMEDIR(1:LRED)
            LOCALFILE(LRED+1:LRED+8)='/files/t'
            LOCALFILE(LRED+9:)=TABFILE(I:I+K-1)
            INQUIRE(FILE=LOCALFILE,EXIST=LOGFILE)
            IF(LOGFILE)THEN
              TABFILE=LOCALFILE
              FINDTABFILE=.TRUE.
              RETURN
            END IF
          END DO
          K=K+1
        END DO
C
        END
