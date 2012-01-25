C------------------------------------------------------------------------------
C Version 29-October-2008                                       file: vaucoul.f
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
C Program: vaucoul
C Classification: miscellany
C Description: Fits a de Vaucouleurs profile to a galaxy frame. With the
C previous fit, the program also computes a Sersic profile.
C
Comment
C
C Realiza un ajuste a un perfile de de Vaucouleurs mas una recta (cielo)
C
        PROGRAM VAUCOUL
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER TRUELEN
        INTEGER READILIM
        REAL READF
C
        INTEGER NB
        INTEGER I,J,K,L
        INTEGER N1,N2,N0
        INTEGER NADD,NFIT,NPLOT
        INTEGER NC1,NC2,NDC,NINTERV
        INTEGER NTERMS                          !grado de polinomio para POLFIT
        INTEGER NSTEP,NMEAN
        INTEGER NDIM,NEVAL
        INTEGER NCMAX_LOCAL,NSMAX_LOCAL
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        INTEGER NDEG_SKY
        REAL A(NCMAX,NSMAX),RES(NCMAX)
        REAL S(NCMAX),X(NCMAX)            !dimensionado al mayor de NCMAX,NSMAX
        REAL SS(NCMAX)
        REAL XMIN,XMAX,YMIN,YMAX,DX,DY
        REAL XMINV,XMAXV,YMINV,YMAXV,DXV,DYV
        REAL XMINV2,XMAXV2,YMINV2,YMAXV2                       !para intervalos
        REAL XC,YC
        REAL B(2),CHISQR   !coeficientes polinomios en POLFIT y otras variables
        REAL SERSIC_B(2)
        REAL BV(2)                             !ajuste recta en perfil r**(1/4)
        REAL BB(2),BV2(2)                     !igual que B y BV para intervalos
        REAL XFIT(NCMAX),YFIT(NCMAX)      !dimensionado al mayor de NCMAX,NSMAX
        REAL YFIT_SERSIC(NCMAX)
        REAL XP(NCMAX),YP(NCMAX)          !dimensionado al mayor de NCMAX,NSMAX
        REAL Y1,Y2,MEAN,SIGMA
        REAL MEAN_,SIGMA_
        REAL PESOGAL,PESOSKY
        REAL X0,RE,I0,MINRAD
        REAL SERSIC_I0,SERSIC_RE,SERSIC_BN,SERSIC_N
        REAL RE2,I02                         !igual que RE y I0 para intervalos
        EXTERNAL FUNKVAUC_0,FUNKVAUC_1
        REAL FUNKVAUC_0,FUNKVAUC_1
        EXTERNAL FUNKSERC_0,FUNKSERC_1
        REAL FUNKSERC_0,FUNKSERC_1
        REAL YRMSTOL,XX0(6),DXX0(6),XXF(6),DXXF(6)
        REAL STWV2,DISP2
        CHARACTER*1 CH,CREP,CMOUSE,CFLUX
        CHARACTER*50 CDUMMY,XLABEL,YLABEL
        CHARACTER*75 INFILE,FLUXFILE,OUTFILE
        LOGICAL IFCHAN(NCMAX)
        LOGICAL IFX(NCMAX)                          !region a ajustar el perfil
        LOGICAL IFXSKY(NCMAX)                        !region a ajustar el cielo
        LOGICAL IFFUNK1(NCMAX)                     !region a minimizar de cielo
        LOGICAL IFFUNK2(NCMAX)                   !region a minimizar de galaxia
        LOGICAL LFIRSTPLOT,LINIFIT,LREFINED
        LOGICAL LBEXIST
        LOGICAL LCOLOR(MAX_ID_RED)
C
        COMMON/BLKFUNK0/NSCAN
        COMMON/BLKFUNK1/X0
        COMMON/BLKFUNK2/S
        COMMON/BLKFUNK3/IFFUNK1,IFFUNK2
        COMMON/BLKFUNK4/PESOSKY,PESOGAL
C------------------------------------------------------------------------------
        THISPROGRAM='vaucoul'
        CALL WELCOME('20-September-2007')
C
        NSMAX_LOCAL=NSMAX
        NCMAX_LOCAL=NCMAX
        IF(NSMAX_LOCAL.GT.NCMAX_LOCAL)STOP 'FATAL ERROR: NSMAX.GT.NCMAX'
C------------------------------------------------------------------------------
C Valores iniciales
        LFIRSTPLOT=.TRUE.
        LINIFIT=.FALSE.
        LREFINED=.FALSE.
        CMOUSE='m'
        PESOGAL=1.
        PESOSKY=1.
C evita warnings de compilacion
        NDEG_SKY=0
        MINRAD=0.0
        SERSIC_B(1)=0.0
        SERSIC_B(2)=0.0
        SERSIC_RE=0.0
        SERSIC_N=0.0
        SERSIC_BN=0.0
C------------------------------------------------------------------------------
        CALL RPGBEGIN(NTERM,IDN,LCOLOR)
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGSLW(2)
        END DO
C------------------------------------------------------------------------------
        WRITE(*,100)'Input file name'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        DO I=1,NSCAN
          READ(20) (A(J,I),J=1,NCHAN)
        END DO
        CLOSE(20)
C------------------------------------------------------------------------------
        WRITE(*,101)'* Flux calibration:'
        WRITE(*,100)'Are you using a response curve (y/n) '
        CFLUX(1:1)=READC('n','yn')
        IF(CFLUX.EQ.'y')THEN
          WRITE(*,100)'Response curve file name'
          FLUXFILE=INFILEX(30,'../stand/curvresf',N1,N2,
     +     STWV2,DISP2,1,.FALSE.)
          IF((N2.NE.NCHAN).OR.(STWV2.NE.STWV).OR.(DISP2.NE.DISP))THEN
            WRITE(*,101)'FATAL ERROR: NCHAN/STWV/DISP is(are) '//
     +       'different.'
            STOP
          END IF
          IF(N1.GT.1)THEN
            WRITE(*,101)'WARNING: this file contains more than '//
     +       '1 spectrum.'
            WRITE(*,101)'The first spectrum will be employed as '//
     +       'response curve.'
          END IF
          READ(30) (RES(J),J=1,NCHAN)
          CLOSE(30)
          DO I=1,NSCAN
            DO J=1,NCHAN
              A(J,I)=A(J,I)/RES(J)
            END DO
          END DO
          WRITE(*,*)
        END IF
C------------------------------------------------------------------------------
        DO J=1,NCHAN
          IFCHAN(J)=.FALSE.
        END DO
        WRITE(*,101)'* Define channel region to be added:'
        WRITE(CDUMMY,*)NCHAN
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,101)'  (valid range is: 1,'//CDUMMY(1:L)//')'
10      WRITE(*,100)'Channels (0,0 = EXIT) '
        CALL READ2I('0,0',NC1,NC2)
        IF((NC1.EQ.0).AND.(NC2.EQ.0))THEN
        ELSE
          IF((NC1.LT.1).OR.(NC2.GT.NCHAN).OR.(NC1.GT.NC2))THEN
            WRITE(*,101)'ERROR: numbers out of range. Try again.'
          ELSE
            DO J=NC1,NC2
              IFCHAN(J)=.TRUE.
            END DO
          END IF
          GOTO 10
        END IF
        NADD=0
        DO J=1,NCHAN
          IF(IFCHAN(J)) NADD=NADD+1
        END DO
        IF(NADD.EQ.0)THEN
          WRITE(*,101)'ERROR: No. of channels added = 0!'
          GOTO 10
        ELSE
          WRITE(*,100)'No. of channels added: '
          WRITE(*,*)NADD
        END IF
        DO I=1,NSCAN
          S(I)=0.
          DO J=1,NCHAN
            IF(IFCHAN(J)) S(I)=S(I)+A(J,I)
          END DO
          S(I)=S(I)/REAL(NADD)
        END DO
C------------------------------------------------------------------------------
        CALL BUTTON(1,'[I]ni. fit',0)
        CALL BUTTON(2,'[W]hole',0)
        CALL BUTTON(3,'[Z]oom',0)
        CALL BUTTON(4,'[Y]-limits',0)
        CALL BUTTON(5,'postScript',0)
        CALL BUTTON(5,'postScript',3)
        CALL BUTTON(6,'[E]xit',0)
        CALL BUTTON(7,'[R]efine',0)
        CALL BUTTON(7,'[R]efine',3)
        CALL BUTTON(8,'[S]tatistics',0)
        CALL BUTTON(8,'[S]tatistics',3)
        CALL BUTTON(9,'Mouse',0)
        CALL BUTTON(10,'Wei[g]hts',0)
        CALL BUTTON(11,'Sa[v]e...',0)
        CALL BUTTON(11,'Sa[v]e...',3)
C------------------------------------------------------------------------------
C Dibujamos corte realizado
        DO I=1,NSCAN
          X(I)=REAL(I)
        END DO
        XMIN=1.
        XMAX=REAL(NSCAN)
        DX=XMAX-XMIN
        XMIN=XMIN-DX/50.
        XMAX=XMAX+DX/50.
C
        XLABEL='scan'
        YLABEL='No. of counts'
        N1=1
        N2=NSCAN
C
40      YMIN=S(N1)
        DO I=N1+1,N2
          IF(S(I).LT.YMIN) YMIN=S(I)
        END DO
        YMAX=YMIN
        DO I=N1,N2
          IF(S(I).GT.YMAX) YMAX=S(I)
        END DO
        DY=YMAX-YMIN
        YMIN=YMIN-DY/30.
        YMAX=YMAX+DY/30.
C
42      IF(LFIRSTPLOT)THEN
          LFIRSTPLOT=.FALSE.
        END IF
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
ccc       CALL PGSLW(2)
          IF(ITERM.EQ.1)THEN
            CALL RPGERASW(0.,1.,0.,0.80)
            CALL RPGENV(XMIN,XMAX,YMIN,YMAX,0,0)
          ELSE
            CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,0)
          END IF
ccc       CALL PGSLW(1)
          CALL PGIDEN_RED
          CALL PGBIN(NSCAN,X,S,.TRUE.)
          CALL PGLABEL(XLABEL,YLABEL,CHAR(32))
          CALL PGMTEXT('T',1.5,0.,0.,'File: '//INFILE)
          CALL PGMTEXT('T',1.5,1.,1.,OBJECT(1:TRUELEN(OBJECT)))
        END DO
C------------------------------------------------------------------------------
C Si existe ajuste inicial/refinado, lo dibujamos
        IF(LINIFIT)THEN
          CALL BUTTON(11,'Sa[v]e...',0)
          WRITE(CDUMMY,*)RE
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGSCH(0.8)
            CALL PGMTEXT('T',-2.,0.02,0.,'r\de\u: '//CDUMMY(1:L))
            CALL PGSCH(1.0)
          END DO
          DO I=1,NSCAN
            YFIT(I)=B(1)+B(2)*X(I)
            YFIT_SERSIC(I)=SERSIC_B(1)+SERSIC_B(2)*X(I)
          END DO
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGSLS(4)
            CALL PGLINE(NSCAN,X,YFIT)
            CALL PGSLS(2)
            CALL PGSCI(14)
            CALL PGLINE(NSCAN,X,YFIT_SERSIC)
            CALL PGSLS(1)
          END DO
          IF(LREFINED)THEN
            CALL BUTTON(8,'[S]tatistics',0)
          ELSE
            CALL BUTTON(7,'[R]efine',0)
          END IF
          DO I=1,NSCAN
            Y2=I0*10**(-3.33*((ABS(X(I)-X0)/RE)**0.25-1.0))
            YFIT(I)=B(1)+B(2)*X(I)+Y2
            Y2=I0*EXP(-SERSIC_BN*
     >       ((ABS(X(I)-X0)/SERSIC_RE)**(1./SERSIC_N)-1.0))
            YFIT_SERSIC(I)=B(1)+B(2)*X(I)+Y2
          END DO
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(3)
            CALL PGSLS(4)
            CALL PGLINE(NSCAN,X,YFIT)
            CALL PGSLS(2)
            CALL PGSCI(14)
            CALL PGLINE(NSCAN,X,YFIT_SERSIC)
            CALL PGSCI(3)
            IF(LCOLOR(ITERM)) CALL PGSCI(4)
            CALL PGSLS(2)
            CALL PGMOVE(X0-MINRAD,YMIN)
            CALL PGDRAW(X0-MINRAD,YMAX)
            CALL PGMOVE(X0+MINRAD,YMIN)
            CALL PGDRAW(X0+MINRAD,YMAX)
            CALL PGSLS(1)
            DO I=1,NSCAN
              IF(IFXSKY(I))THEN
                IF(LCOLOR(ITERM)) CALL PGSCI(5)
                CALL PGPOINT(1,X(I),S(I),17)
              ELSEIF(IFX(I))THEN
                IF(LCOLOR(ITERM)) CALL PGSCI(2)
                CALL PGPOINT(1,X(I),S(I),17)
              END IF
            END DO
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END DO
        END IF
        IF(LREFINED)THEN
          WRITE(*,*)
          WRITE(*,101)'* Statistics in the sky region: '
          MEAN=0.  !r**(1/4)
          MEAN_=0. !Sersic
          NMEAN=0
          DO I=1,NSCAN
            IF(IFXSKY(I))THEN
              NMEAN=NMEAN+1
              Y2=I0*10**(-3.33*((ABS(X(I)-X0)/RE)**0.25-1.0))
              MEAN=MEAN+Y2
              Y2=I0*EXP(-SERSIC_BN*
     >         ((ABS(X(I)-X0)/SERSIC_RE)**(1./SERSIC_N)-1.0))
              MEAN_=MEAN_+Y2
            END IF
          END DO
          MEAN=MEAN/REAL(NMEAN)
          MEAN_=MEAN_/REAL(NMEAN)
          SIGMA=0.  !r**(1/4)
          SIGMA_=0. !Sersic
          DO I=1,NSCAN
            IF(IFXSKY(I))THEN
              Y2=I0*10**(-3.33*((ABS(X(I)-X0)/RE)**0.25-1.0))
              SIGMA=SIGMA+(Y2-MEAN)*(Y2-MEAN)
              Y2=I0*EXP(-SERSIC_BN*
     >         ((ABS(X(I)-X0)/SERSIC_RE)**(1./SERSIC_N)-1.0))
              SIGMA_=SIGMA_+(Y2-MEAN_)*(Y2-MEAN_)
            END IF
          END DO
          SIGMA=SQRT(SIGMA/(NMEAN-1))
          SIGMA_=SQRT(SIGMA_/(NMEAN-1))
          WRITE(*,100) 'De Vaucouleurs> Mean no. of counts (galaxy): '
          WRITE(*,*)MEAN
          WRITE(*,100) 'De Vaucouleurs> Std. deviation.... (galaxy): '
          WRITE(*,*)SIGMA
          WRITE(*,100) 'Sersic profile> Mean no. of counts (galaxy): '
          WRITE(*,*)MEAN_
          WRITE(*,100) 'Sersic profile> Std. deviation.... (galaxy): '
          WRITE(*,*)SIGMA_
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGSCH(0.8)
            WRITE(CDUMMY,*)MEAN
            CALL RMBLANK(CDUMMY,CDUMMY,L)
            CALL PGMTEXT('T',-2.,.98,1.,
     +       'Mean N\dgalaxy\u(sky region): '//CDUMMY(1:L))
            WRITE(CDUMMY,*)SIGMA
            CALL RMBLANK(CDUMMY,CDUMMY,L)
            CALL PGMTEXT('T',-4.,.98,1.,
     +       'Standard Deviation: '//CDUMMY(1:L))
            CALL PGSCH(1.0)
          END DO
        END IF
C------------------------------------------------------------------------------
50      CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
        CALL IFBUTTON(XC,YC,NB)
        IF(NB.EQ.0)THEN
          CALL CHLOWER(CH)
          IF(CH.EQ.'i')THEN
            NB=1
          ELSEIF(CH.EQ.'w')THEN
            NB=2
          ELSEIF(CH.EQ.'z')THEN
            NB=3
          ELSEIF(CH.EQ.'y')THEN
            NB=4
          ELSEIF(CH.EQ.'e')THEN
            NB=6
          ELSEIF(CH.EQ.'r')THEN
            CALL BUTTQEX(7,LBEXIST)
            IF(LBEXIST) NB=7
          ELSEIF(CH.EQ.'s')THEN
            CALL BUTTQEX(8,LBEXIST)
            IF(LBEXIST) NB=8
          ELSEIF(CH.EQ.'g')THEN
            NB=10
          ELSEIF(CH.EQ.'v')THEN
            NB=11
          END IF
        END IF
C------------------------------------------------------------------------------
        IF(NB.EQ.0)THEN
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.1)THEN
          CALL BUTTON(1,'[I]ni. fit',5)
          IF(LINIFIT)THEN
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
ccc           CALL PGSLW(2)
              IF(ITERM.EQ.1)THEN
                CALL RPGERASW(0.,1.,0.,0.80)
                CALL RPGENV(XMIN,XMAX,YMIN,YMAX,0,0)
              ELSE
                CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,0)
              END IF
ccc           CALL PGSLW(1)
              CALL PGIDEN_RED
              CALL PGBIN(NSCAN,X,S,.TRUE.)
              CALL PGLABEL(XLABEL,YLABEL,CHAR(32))
              CALL PGMTEXT('T',1.5,0.,0.,'File: '//INFILE)
              CALL PGMTEXT('T',1.5,1.,1.,OBJECT(1:TRUELEN(OBJECT)))
            END DO
            LINIFIT=.FALSE.
          END IF
C..............................................................................
C cielo inicial
          WRITE(*,101)'* Select region for a first '//
     +     'approximation to the sky level:'
          DO I=1,NSCAN
            IFXSKY(I)=.FALSE.
          END DO
          WRITE(*,100)'Polynomial degree for sky level (0 or 1) '
          NDEG_SKY=READILIM('0',0,1)
60        WRITE(*,101)'Select region to be used:'
          IF(CMOUSE.EQ.'k')THEN
            WRITE(*,100)'Region (0,0=EXIT) '
            CALL READ2I('0,0',N1,N2)
            IF((N1.EQ.0).AND.(N2.EQ.0)) GOTO 62
            IF(N1.LT.1) N1=1
            IF(N1.GT.NSCAN) N1=NSCAN
            IF(N2.LT.1) N2=1
            IF(N2.GT.NSCAN) N2=NSCAN
            IF(N1.GT.N2)THEN        !ordenamos los extremos si estan invertidos
              N0=N1
              N1=N2
              N2=N0
            END IF
            GOTO 61
          END IF
          WRITE(*,100)'Press mouse button...'
          IF(LCOLOR(1)) CALL PGSCI(5)
          CALL RPGBAND(6,0,0.,0.,XC,YC,CH)
          IF(LCOLOR(1)) CALL PGSCI(1)
          IF(CH.EQ.'X') GOTO 62
          N1=NINT(XC)
          IF(N1.LT.1) N1=1
          IF(N1.GT.NSCAN) N1=NSCAN
          WRITE(*,*) N1
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(2)
            CALL PGMOVE(REAL(N1)-.5,S(N1))
            CALL PGDRAW(REAL(N1)+.5,S(N1))
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END DO
          WRITE(*,100)'Press mouse button...'
          IF(LCOLOR(1)) CALL PGSCI(5)
          CALL RPGBAND(4,0,REAL(N1),0.,XC,YC,CH)
          IF(LCOLOR(1)) CALL PGSCI(1)
          IF(CH.EQ.'X') GOTO 62
          N2=NINT(XC)
          IF(N2.LT.1) N2=1
          IF(N2.GT.NSCAN) N2=NSCAN
          WRITE(*,*) N2
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(2)
            CALL PGMOVE(REAL(N2)-.5,S(N2))
            CALL PGDRAW(REAL(N2)+.5,S(N2))
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END DO
          IF(N1.GT.N2)THEN          !ordenamos los extremos si estan invertidos
            N0=N1
            N1=N2
            N2=N0
          END IF
61        NFIT=0
          DO I=N1,N2
            IFXSKY(I)=.TRUE.
            NFIT=NFIT+1
            XFIT(NFIT)=X(I)
            YFIT(NFIT)=S(I)              !solo para dibujar region seleccionada
          END DO
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(5)
            CALL PGBIN(NFIT,XFIT,YFIT,.TRUE.)
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END DO
          GOTO 60
62        WRITE(*,*)
          NFIT=0
          DO I=1,NSCAN
            IF(IFXSKY(I))THEN
              NFIT=NFIT+1
              XFIT(NFIT)=X(I)
              YFIT(NFIT)=S(I)
            END IF
          END DO
          IF(NFIT.LT.2)THEN
            WRITE(*,101)'ERROR: No. of points for fit < 2'
            CALL BUTTON(1,'[I]ni. fit',0)
            GOTO 50
          END IF
          NTERMS=NDEG_SKY+1                                !ajustamos una recta
          CALL POLFIT(XFIT,YFIT,YFIT,NFIT,NTERMS,0,B,CHISQR)
          IF(NTERMS.EQ.1) B(2)=0.0
          SERSIC_B(1)=B(1)
          SERSIC_B(2)=B(2)
          WRITE(*,*)
          WRITE(*,101)'Coefficients from fit: (y=a0+a1*x)'
          DO K=1,NTERMS
            WRITE(*,'(A,I1,A,$)')'> a(',K-1,') : '
            WRITE(*,*)B(K)
          END DO
          DO I=1,NSCAN
            Y2=B(NTERMS)
            DO K=NTERMS-1,1,-1
              Y2=Y2*X(I)+B(K)
            END DO
            YFIT(I)=Y2
          END DO
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(3)
            CALL PGBIN(NSCAN,X,YFIT,.TRUE.)
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END DO
C..............................................................................
C perfil inicial de la galaxia ajustado a ley r**(1/4)
          WRITE(*,101)'* Select region for a first '//
     +     'approximation to the galaxy profile:'
          DO I=1,NSCAN
            IFX(I)=.FALSE.
          END DO
70        WRITE(*,101)'Select region to be used:'
          IF(CMOUSE.EQ.'k')THEN
            WRITE(*,100)'Region (0,0=EXIT) '
            CALL READ2I('0,0',N1,N2)
            IF((N1.EQ.0).AND.(N2.EQ.0)) GOTO 72
            IF(N1.LT.1) N1=1
            IF(N1.GT.NSCAN) N1=NSCAN
            IF(N2.LT.1) N2=1
            IF(N2.GT.NSCAN) N2=NSCAN
            IF(N1.GT.N2)THEN        !ordenamos los extremos si estan invertidos
              N0=N1
              N1=N2
              N2=N0
            END IF
            GOTO 71
          END IF
          WRITE(*,100)'Press mouse button...'
          IF(LCOLOR(1)) CALL PGSCI(5)
          CALL RPGBAND(6,0,0.,0.,XC,YC,CH)
          IF(LCOLOR(1)) CALL PGSCI(1)
          IF(CH.EQ.'X') GOTO 72
          N1=NINT(XC)
          IF(N1.LT.1) N1=1
          IF(N1.GT.NSCAN) N1=NSCAN
          WRITE(*,*) N1
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(2)
            CALL PGMOVE(REAL(N1)-.5,S(N1))
            CALL PGDRAW(REAL(N1)+.5,S(N1))
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END DO
          WRITE(*,100)'Press mouse button...'
          IF(LCOLOR(1)) CALL PGSCI(5)
          CALL RPGBAND(4,0,REAL(N1),0.,XC,YC,CH)
          IF(LCOLOR(1)) CALL PGSCI(1)
          IF(CH.EQ.'X') GOTO 72
          N2=NINT(XC)
          IF(N2.LT.1) N2=1
          IF(N2.GT.NSCAN) N2=NSCAN
          WRITE(*,*) N2
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(2)
            CALL PGMOVE(REAL(N2)-.5,S(N2))
            CALL PGDRAW(REAL(N2)+.5,S(N2))
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END DO
          IF(N1.GT.N2)THEN          !ordenamos los extremos si estan invertidos
            N0=N1
            N1=N2
            N2=N0
          END IF
71        NFIT=0
          DO I=N1,N2
            IFX(I)=.TRUE.
            NFIT=NFIT+1
            XFIT(NFIT)=X(I)
            YFIT(NFIT)=S(I)              !solo para dibujar region seleccionada
          END DO
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(2)
            CALL PGBIN(NFIT,XFIT,YFIT,.TRUE.)
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END DO
          GOTO 70
C
72        WRITE(*,*)
          WRITE(*,100)'X0'
          X0=READF('@')
C
73        NFIT=0
          DO I=1,NSCAN
            IF(IFX(I))THEN
              NFIT=NFIT+1
              XFIT(NFIT)=X(I)
              Y2=B(NTERMS)
              DO K=NTERMS-1,1,-1
                Y2=Y2*X(I)+B(K)
              END DO
              YFIT(NFIT)=S(I)-Y2               !(galaxia+cielo)-(cielo inicial)
            END IF
          END DO
          IF(NFIT.LT.5)THEN
            WRITE(*,101)'ERROR: No. of points for fit < 5'
            CALL BUTTON(1,'[I]ni. fit',0)
            GOTO 50
          END IF
C pasamos a escala logaritmica y r**(1/4)
          DO I=1,NFIT
            XP(I)=ABS(XFIT(I)-X0)**0.25
            IF(YFIT(I).LE.0.)THEN
              WRITE(*,101)'ERROR: log(I), I<0.0!'
              WRITE(*,100)'Press <CR> to continue..'
              READ(*,*)
              CALL BUTTON(1,'[I]ni. fit',0)
              GOTO 50
            END IF
            YP(I)=ALOG10(YFIT(I))
          END DO
          NPLOT=NFIT
          XMINV=XP(1)
          XMAXV=XMINV
          YMINV=YP(1)
          YMAXV=YMINV
          DO I=2,NPLOT
            IF(XP(I).LT.XMINV) XMINV=XP(I)
            IF(XP(I).GT.XMAXV) XMAXV=XP(I)
            IF(YP(I).LT.YMINV) YMINV=YP(I)
            IF(YP(I).GT.YMAXV) YMAXV=YP(I)
          END DO
          DXV=XMAXV-XMINV
          DYV=YMAXV-YMINV
          XMINV=XMINV-DXV/30.
          XMAXV=XMAXV+DXV/30.
          YMINV=YMINV-DYV/30.
          YMAXV=YMAXV+DYV/30.
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
ccc         CALL PGSLW(2)
            IF(ITERM.EQ.1)THEN
              CALL RPGERASW(0.,1.,0.,0.80)
              CALL RPGENV(XMINV,XMAXV,YMINV,YMAXV,0,0)
            ELSE
              CALL PGENV(XMINV,XMAXV,YMINV,YMAXV,0,0)
            END IF
ccc         CALL PGSLW(1)
            CALL PGIDEN_RED
            CALL PGPOINT(NPLOT,XP,YP,17)
            CALL PGLABEL('r\u1/4\d','log(I)',
     +       'fit to a de Vaucouleurs profile')
            CALL PGMTEXT('T',1.5,0.,0.,'File: '//INFILE)
            CALL PGMTEXT('T',1.5,1.,1.,OBJECT(1:TRUELEN(OBJECT)))
          END DO
C eliminamos puntos centrales para el ajuste a la recta
          IF(CMOUSE.EQ.'k')THEN
            WRITE(*,100)'Minimum r**(.25) value '
            XC=READF('@')
          ELSE
            WRITE(*,100)'* Press mouse to remove inner points...'
            IF(LCOLOR(1)) CALL PGSCI(5)
            CALL RPGBAND(7,0,0.,0.,XC,YC,CH)
            IF(LCOLOR(1)) CALL PGSCI(1)
            WRITE(*,*)
          END IF
          MINRAD=XC*XC*XC*XC
          NFIT=0
          DO I=1,NPLOT
            IF(XP(I).GE.(MINRAD)**.25)THEN
              NFIT=NFIT+1
              XFIT(NFIT)=XP(I)
              YFIT(NFIT)=YP(I)
            END IF
          END DO
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(2)
            CALL PGPOINT(NFIT,XFIT,YFIT,17)
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END DO
C ajustamos recta
          IF(NFIT.LT.2)THEN
            WRITE(*,101)'ERROR: No. of points for fit < 2'
            CALL BUTTON(1,'[I]ni. fit',0)
            GOTO 50
          END IF
          NTERMS=2                                         !ajustamos una recta
          CALL POLFIT(XFIT,YFIT,YFIT,NFIT,NTERMS,0,BV,CHISQR)
          WRITE(*,*)
          WRITE(*,101)'Coefficients from fit: (y=a0+a1*x)'
          DO K=1,NTERMS
            WRITE(*,'(A,I1,A,$)')'> a(',K-1,') : '
            WRITE(*,*)BV(K)
          END DO
          DO I=1,NSCAN
            Y2=BV(NTERMS)
            DO K=NTERMS-1,1,-1
              Y2=Y2*X(I)+BV(K)
            END DO
            YFIT(I)=Y2
          END DO
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(3)
            CALL PGLINE(NSCAN,X,YFIT)
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END DO
C mostramos parametros iniciales del ajuste al perfil r**(1/4)
          RE=(-3.33/BV(2))*(-3.33/BV(2))*(-3.33/BV(2))*(-3.33/BV(2))
          I0=10.0**(BV(1)-3.33)
          WRITE(*,101)'Coefficients: '//
     +     'y=Io*10**{-3.33*[(abs(x-x0)/Re)**0.25-1]}'
          WRITE(*,100)'> Io : '
          WRITE(*,*)I0
          WRITE(*,100)'> x0 : '
          WRITE(*,*)X0
          WRITE(*,100)'> Re : '
          WRITE(*,*)RE
          WRITE(*,*)
C..............................................................................
C Repetimos el proceso a intervalos. Se realiza en dos pasos. En el primer
C paso hacemos los calculos y en el segundo realizamos los dibujos.
          WRITE(*,100)'Repeat fit in small channel regions (y/n) '
          CREP(1:1)=READC('n','yn')
          IF(CREP.EQ.'y')THEN
            WRITE(*,100)'First channel '
            NC1=READILIM('1',1,NCHAN)
            WRITE(*,100)'Last channel '
            WRITE(CDUMMY,*)NCHAN
            NC2=READILIM(CDUMMY,NC1,NCHAN)
            WRITE(*,100)'Increment (channels) '
            NDC=READILIM('@',1,NC2-NC1+1)
            NINTERV=(NC2-NC1+1)/NDC
            IF(MOD(NC2-NC1+1,NDC).NE.0) NINTERV=NINTERV+1
            NSTEP=1
            XMINV2=1.E6
            XMAXV2=-1.E6
            YMINV2=1.E6
            YMAXV2=-1.E6
80          IF(NSTEP.EQ.2)THEN
              DXV=XMAXV2-XMINV2
              DYV=YMAXV2-YMINV2
              XMINV2=XMINV2-DXV/30.
              XMAXV2=XMAXV2+DXV/30.
              YMINV2=YMINV2-DYV/30.
              YMAXV2=YMAXV2+DYV/30.
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
ccc             CALL PGSLW(2)
                IF(ITERM.EQ.1)THEN
                  CALL RPGERASW(0.,1.,0.,0.80)
                  CALL RPGENV(XMINV2,XMAXV2,YMINV2,YMAXV2,0,0)
                ELSE
                  CALL PGENV(XMINV2,XMAXV2,YMINV2,YMAXV2,0,0)
                END IF
ccc             CALL PGSLW(1)
                CALL PGIDEN_RED
                CALL PGLABEL('r\u1/4\d','log(I)',
     +           'fit to a de Vaucouleurs profile')
              END DO
            END IF
            DO L=1,NINTERV                       !bucle en numero de intervalos
              N1=NC1+(L-1)*NDC
              N2=NC1+L*NDC
              IF(N2.GT.NC2) N2=NC2
              DO I=1,NSCAN
                SS(I)=0.
                DO J=N1,N2
                  SS(I)=SS(I)+A(J,I)
                END DO
                SS(I)=SS(I)/REAL(N2-N1+1)
              END DO
              NFIT=0
              DO I=1,NSCAN
                IF(IFXSKY(I))THEN
                  NFIT=NFIT+1
                  XFIT(NFIT)=X(I)
                  YFIT(NFIT)=SS(I)
                END IF
              END DO
              NTERMS=NDEG_SKY+1                            !ajustamos una recta
              CALL POLFIT(XFIT,YFIT,YFIT,NFIT,NTERMS,0,BB,CHISQR)
              IF(NTERMS.EQ.2) BB(2)=0.0
              NFIT=0
              DO I=1,NSCAN
                IF(IFX(I))THEN
                  NFIT=NFIT+1
                  XFIT(NFIT)=X(I)
                  Y2=BB(NTERMS)
                  DO K=NTERMS-1,1,-1
                    Y2=Y2*X(I)+BB(K)
                  END DO
                  YFIT(NFIT)=SS(I)-Y2          !(galaxia+cielo)-(cielo inicial)
                END IF
              END DO
              DO I=1,NFIT
                XP(I)=ABS(XFIT(I)-X0)**0.25
                YP(I)=ALOG10(YFIT(I))
              END DO
              NPLOT=NFIT
              NFIT=0
              DO I=1,NPLOT
                IF(XP(I).GE.(MINRAD)**.25)THEN
                  NFIT=NFIT+1
                  XFIT(NFIT)=XP(I)
                  YFIT(NFIT)=YP(I)
                END IF
              END DO
              IF(NSTEP.EQ.1)THEN
                DO I=1,NFIT
                  IF(XMINV2.GT.XFIT(I)) XMINV2=XFIT(I)
                  IF(XMAXV2.LT.XFIT(I)) XMAXV2=XFIT(I)
                  IF(YMINV2.GT.YFIT(I)) YMINV2=YFIT(I)
                  IF(YMAXV2.LT.YFIT(I)) YMAXV2=YFIT(I)
                END DO
                GOTO 82
              END IF
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                IF(LCOLOR(ITERM)) CALL PGSCI(L)
                CALL PGPOINT(NFIT,XFIT,YFIT,L)
                IF(LCOLOR(ITERM)) CALL PGSCI(1)
              END DO
              NTERMS=2                                     !ajustamos una recta
              CALL POLFIT(XFIT,YFIT,YFIT,NFIT,NTERMS,0,BV2,CHISQR)
              DO I=1,NSCAN
                Y2=BV2(NTERMS)
                DO K=NTERMS-1,1,-1
                  Y2=Y2*X(I)+BV2(K)
                END DO
                YFIT(I)=Y2
              END DO
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                IF(LCOLOR(ITERM)) CALL PGSCI(L)
                CALL PGLINE(NSCAN,X,YFIT)
                IF(LCOLOR(ITERM)) CALL PGSCI(1)
              END DO
              RE2=(-3.33/BV2(2))*(-3.33/BV2(2))*
     +         (-3.33/BV2(2))*(-3.33/BV2(2))
              I02=10.0**(BV2(1)-3.33)
              WRITE(*,*)L,RE2,I02
82            CONTINUE
            END DO                                 !fin en numero de intervalos
            IF(NSTEP.EQ.1)THEN
              NSTEP=2
              GOTO 80
            END IF
          END IF
C..............................................................................
          WRITE(*,100)'Repeat last fit (y/n) '
          CREP(1:1)=READC('n','yn')
          IF(CREP.EQ.'y') GOTO 73
C..............................................................................
C hemos terminado
          LINIFIT=.TRUE.
          N1=1
          N2=NSCAN
          CALL BUTTON(1,'[I]ni. fit',0)
          GOTO 40
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.2)THEN
          CALL BUTTON(2,'[W]hole',5)
          XMIN=1.
          XMAX=REAL(NSCAN)
          XMIN=XMIN-DX/50.
          XMAX=XMAX+DX/50.
          N1=1
          N2=NSCAN
          CALL BUTTON(2,'[W]hole',0)
          GOTO 40
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.3)THEN
          CALL BUTTON(3,'[Z]oom',5)
          IF(CMOUSE.EQ.'k')THEN
            WRITE(*,100)'1st scan '
            N1=READILIM('@',1,NSCAN)
          ELSE
            WRITE(*,100)'Press mouse button...'
            IF(LCOLOR(1)) CALL PGSCI(5)
            CALL RPGBAND(6,0,0.,0.,XC,YC,CH)
            IF(LCOLOR(1)) CALL PGSCI(1)
            N1=NINT(XC)
            WRITE(*,*) N1
          END IF
          IF(CMOUSE.EQ.'k')THEN
            WRITE(*,100)'2nd scan '
            N2=READILIM('@',1,NSCAN)
          ELSE
            WRITE(*,100)'Press mouse button...'
            IF(LCOLOR(1)) CALL PGSCI(5)
            CALL RPGBAND(4,0,REAL(N1),0.,XC,YC,CH)
            IF(LCOLOR(1)) CALL PGSCI(1)
            N2=NINT(XC)
            WRITE(*,*) N2
          END IF
          IF(N1.GT.N2)THEN
            N0=N1
            N1=N2
            N2=N0
          END IF
          IF(N1.LT.1) N1=1
          IF(N2.GT.NSCAN) N2=NSCAN
          XMIN=REAL(N1)
          XMAX=REAL(N2)
          CALL BUTTON(3,'[Z]oom',0)
          GOTO 40
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.4)THEN
          CALL BUTTON(4,'[Y]-limits',5)
          WRITE(CDUMMY,*)YMIN
          WRITE(*,100)'Ymin: '
          YMIN=READF(CDUMMY)
          WRITE(CDUMMY,*)YMAX
          WRITE(*,100)'Ymax: '
          YMAX=READF(CDUMMY)
          CALL BUTTON(4,'[Y]-limits',0)
          GOTO 42
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.5)THEN
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.6)THEN
          CALL BUTTON(6,'[E]xit',5)
          CALL PGEND
          STOP
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.7)THEN
          CALL BUTTON(7,'[R]efine',5)
          LREFINED=.TRUE.
C..............................................................................
          WRITE(*,*)
          WRITE(*,101)'* NOTE: galaxy and sky regions will be fitted'
          DO I=1,NSCAN
            IFFUNK1(I)=IFXSKY(I)
            IFFUNK2(I)=(IFX(I).AND.(ABS(REAL(I)-X0).GE.MINRAD))
          END DO
          WRITE(*,100)'YRMSTOL for DOWNHILL '
          YRMSTOL=READF('1.E-6')
C..............................................................................
C
          IF(NDEG_SKY.EQ.0)THEN
            !perfil de de Vaucouleurs
            NDIM=3
            XX0(1)=B(1)
            XX0(2)=I0
            XX0(3)=RE
            DO I=1,3
              IF(XX0(I).NE.0.0)THEN
                DXX0(I)=0.05*XX0(I)
              ELSE
                DXX0(I)=1.0
              END IF
            END DO
            CALL DOWNHILL(NDIM,XX0,DXX0,FUNKVAUC_0,
     +       1.0,0.5,2.0,YRMSTOL,XXF,DXXF,NEVAL)
            B(1)=XXF(1)
            I0=XXF(2)
            RE=XXF(3)
            WRITE(*,*)
            WRITE(*,101) '* De Vaucouleurs profile'
            WRITE(*,101) '========================'
            WRITE(*,100)'No. of iterations (DOWNHILL): '
            WRITE(*,*)NEVAL
            WRITE(*,101)'Sky fit y=a0: '
            WRITE(*,100)'a0 (+ rmsDOWNHILL): '
            WRITE(*,*)XXF(1),DXXF(1)
            WRITE(*,100)'Io (+ rmsDOWNHILL): '
            WRITE(*,*)XXF(2),DXXF(2)
            WRITE(*,100)'Re (+ rmsDOWNHILL): '
            WRITE(*,*)XXF(3),DXXF(3)
            !perfil de Sersic (usando como valores iniciales los calculados 
            !en el perfil de de Vaucouleurs
            NDIM=5
            XX0(1)=B(1)
            XX0(2)=I0
            XX0(3)=RE
            XX0(4)=7.67 !valor para r**(1/4)
            XX0(5)=4.00 !valor para r**(1/4)
            DXX0(1)=DXXF(1)
            DXX0(2)=DXXF(2)
            DXX0(3)=DXXF(3)
            DO I=4,5
              DXX0(I)=0.05*XX0(I)
            END DO
            !tenemos que modificar un poco los valores anteriores porque
            !de lo contrario sale exactamente lo mismo (la solución no
            !cambia porque ya hemos encontrado el minimo para el perfil
            !de de Vaucouleurs
            DO I=1,5
              XX0(I)=XX0(I)*1.05 !modificamos un 5%
              IF(DXX0(I).LE.0.0) DXX0(I)=1.0
            END DO
            CALL DOWNHILL(NDIM,XX0,DXX0,FUNKSERC_0,1.0,0.5,2.0,YRMSTOL,
     +       XXF,DXXF,NEVAL)
            SERSIC_B(1)=XXF(1)
            SERSIC_B(2)=0.0
            SERSIC_I0=XXF(2)
            SERSIC_RE=XXF(3)
            SERSIC_BN=XXF(4)
            SERSIC_N=XXF(5)
            WRITE(*,*)
            WRITE(*,101) '* Sersic profile'
            WRITE(*,101) '================'
            WRITE(*,100)'No. of iterations (DOWNHILL): '
            WRITE(*,*)NEVAL
            WRITE(*,101)'Sky fit y=a0: '
            WRITE(*,100)'a0 (+ rmsDOWNHILL): '
            WRITE(*,*)XXF(1),DXXF(1)
            WRITE(*,100)'Io (+ rmsDOWNHILL): '
            WRITE(*,*)XXF(2),DXXF(2)
            WRITE(*,100)'bn (+ rmsDOWNHILL): '
            WRITE(*,*)XXF(4),DXXF(4)
            WRITE(*,100)'Re (+ rmsDOWNHILL): '
            WRITE(*,*)XXF(3),DXXF(3)
            WRITE(*,100)'n  (+ rmsDOWNHILL): '
            WRITE(*,*)XXF(5),DXXF(5)
          ELSE
            !perfil de de Vaucouleurs
            NDIM=4
            XX0(1)=B(1)
            XX0(2)=B(2)
            XX0(3)=I0
            XX0(4)=RE
            DO I=1,4
              IF(XX0(I).NE.0.0)THEN
                DXX0(I)=0.05*XX0(I)
              ELSE
                DXX0(I)=1.0
              END IF
            END DO
            CALL DOWNHILL(NDIM,XX0,DXX0,FUNKVAUC_1,
     +       1.0,0.5,2.0,YRMSTOL,XXF,DXXF,NEVAL)
            B(1)=XXF(1)
            B(2)=XXF(2)
            I0=XXF(3)
            RE=XXF(4)
            WRITE(*,100)'No. of iterations (DOWNHILL): '
            WRITE(*,*)NEVAL
            WRITE(*,101)'Sky fit y=a0+a1*x: '
            WRITE(*,100)'a0 (+ rmsDOWNHILL): '
            WRITE(*,*)XXF(1),DXXF(1)
            WRITE(*,100)'a1 (+ rmsDOWNHILL): '
            WRITE(*,*)XXF(2),DXXF(2)
            WRITE(*,100)'Io (+ rmsDOWNHILL): '
            WRITE(*,*)XXF(3),DXXF(3)
            WRITE(*,100)'Re (+ rmsDOWNHILL): '
            WRITE(*,*)XXF(4),DXXF(4)
            !perfil de Sersic (usando como valores iniciales los calculados 
            !en el perfil de de Vaucouleurs
            NDIM=6
            XX0(1)=B(1)
            XX0(2)=B(2)
            XX0(3)=I0
            XX0(4)=RE
            XX0(5)=7.67 !valor para r**(1/4)
            XX0(6)=4.00 !valor para r**(1/4)
            DXX0(1)=DXXF(1)
            DXX0(2)=DXXF(2)
            DXX0(3)=DXXF(3)
            DXX0(4)=DXXF(4)
            DO I=5,6
              DXX0(I)=0.05*XX0(I)
            END DO
            !tenemos que modificar un poco los valores anteriores porque
            !de lo contrario sale exactamente lo mismo (la solución no
            !cambia porque ya hemos encontrado el minimo para el perfil
            !de de Vaucouleurs
            DO I=1,6
              XX0(I)=XX0(I)*1.05 !modificamos un 5%
              IF(DXX0(I).LE.0.0) DXX0(I)=1.0
            END DO
            CALL DOWNHILL(NDIM,XX0,DXX0,FUNKSERC_1,1.0,0.5,2.0,YRMSTOL,
     +       XXF,DXXF,NEVAL)
            SERSIC_B(1)=XXF(1)
            SERSIC_B(2)=XXF(2)
            SERSIC_I0=XXF(3)
            SERSIC_RE=XXF(4)
            SERSIC_BN=XXF(5)
            SERSIC_N=XXF(6)
            WRITE(*,*)
            WRITE(*,101) '* Sersic profile'
            WRITE(*,101) '================'
            WRITE(*,100)'No. of iterations (DOWNHILL): '
            WRITE(*,*)NEVAL
            WRITE(*,101)'Sky fit y=a0+a1*x: '
            WRITE(*,100)'a0 (+ rmsDOWNHILL): '
            WRITE(*,*)XXF(1),DXXF(1)
            WRITE(*,100)'a1 (+ rmsDOWNHILL): '
            WRITE(*,*)XXF(2),DXXF(2)
            WRITE(*,100)'Io (+ rmsDOWNHILL): '
            WRITE(*,*)XXF(3),DXXF(3)
            WRITE(*,100)'bn (+ rmsDOWNHILL): '
            WRITE(*,*)XXF(5),DXXF(5)
            WRITE(*,100)'Re (+ rmsDOWNHILL): '
            WRITE(*,*)XXF(4),DXXF(4)
            WRITE(*,100)'n  (+ rmsDOWNHILL): '
            WRITE(*,*)XXF(6),DXXF(6)
          END IF
C..............................................................................
          CALL BUTTON(7,'[R]efine',0)
!         N1=1
!         N2=NSCAN
!         GOTO 40
          GOTO 42
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.8)THEN
          CALL BUTTON(8,'[S]tatistics',5)
C..............................................................................
          WRITE(*,101)'* Statistics in the sky region: '
          MEAN=0.
          NMEAN=0
          DO I=1,NSCAN
            IF(IFXSKY(I))THEN
              Y1=B(1)+B(2)*REAL(I)
              Y2=I0*10**(-3.33*((ABS(X(I)-X0)/RE)**0.25-1.0))
              WRITE(*,*)I,Y1,Y1+Y2,Y2
              NMEAN=NMEAN+1
              MEAN=MEAN+Y2
            END IF
          END DO
          MEAN=MEAN/REAL(NMEAN)
          SIGMA=0.
          DO I=1,NSCAN
            IF(IFXSKY(I))THEN
              SIGMA=SIGMA+(Y2-MEAN)*(Y2-MEAN)
            END IF
          END DO
          SIGMA=SQRT(SIGMA/(NMEAN-1))
C
          WRITE(*,100)'Mean no. of counts (galaxy): '
          WRITE(*,*)MEAN
          WRITE(*,100)'Std. deviation.... (galaxy): '
          WRITE(*,*)SIGMA
C..............................................................................
          CALL BUTTON(8,'[S]tatistics',0)
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.9)THEN
          IF(CMOUSE.EQ.'m')THEN
            CALL BUTTON(9,'Keyboard',0)
            CMOUSE='k'
          ELSE
            CALL BUTTON(9,'Mouse',0)
            CMOUSE='m'
          END IF
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.10)THEN
          CALL BUTTON(10,'Wei[g]hts',5)
          WRITE(CDUMMY,*) PESOSKY
          WRITE(*,100)'Weight for sky regions... '
          PESOSKY=READF(CDUMMY)
          WRITE(CDUMMY,*)PESOGAL
          WRITE(*,100)'Weight for galaxy regions '
          PESOGAL=READF(CDUMMY)
          CALL BUTTON(10,'Wei[g]hts',0)
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.11)THEN
          CALL BUTTON(11,'Sa[v]e...',5)
          WRITE(*,100)'Output file name'
          OUTFILE=OUTFILEX(30,'@',NSCAN,1,STWV,DISP,1,.FALSE.)
          DO I=1,NSCAN
            WRITE(30) I0*10**(-3.33*((ABS(X(I)-X0)/RE)**0.25-1.0))
          END DO
          CLOSE(30)
          CALL BUTTON(11,'Sa[v]e...',0)
C------------------------------------------------------------------------------
        ELSE
          WRITE(*,100)'Cursor at: '
          WRITE(*,*)XC,YC
        END IF
C------------------------------------------------------------------------------
        GOTO 50
C------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
        END
C
C******************************************************************************
C Funcion a minimizar: es la suma de un polinomio de primer grado y un 
C perfil de de Vaucouleurs
C B(1)=XX(1), I0=XX(2), RE=XX(3)
C
        REAL FUNCTION FUNKVAUC_0(XX)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        REAL XX(3)
C
        INTEGER I
        REAL S(NCMAX)
        REAL X0
        REAL FF
        REAL PESOSKY,PESOGAL
        LOGICAL IFFUNK1(NCMAX)                     !region a minimizar de cielo
        LOGICAL IFFUNK2(NCMAX)                   !region a minimizar de galaxia
C
        COMMON/BLKFUNK0/NSCAN
        COMMON/BLKFUNK1/X0
        COMMON/BLKFUNK2/S
        COMMON/BLKFUNK3/IFFUNK1,IFFUNK2
        COMMON/BLKFUNK4/PESOSKY,PESOGAL
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        FUNKVAUC_0=0.
C
        DO I=1,NSCAN
          IF(IFFUNK1(I))THEN
            FF=LOG(XX(2))-7.67*((ABS(REAL(I)-X0)/XX(3))**0.25-1.)
            FUNKVAUC_0=FUNKVAUC_0+PESOSKY*
     >       (LOG(S(I)-XX(1))-FF)*(LOG(S(I)-XX(1))-FF)
          END IF
          IF(IFFUNK2(I))THEN
            FF=LOG(XX(2))-7.67*((ABS(REAL(I)-X0)/XX(3))**0.25-1.)
            FUNKVAUC_0=FUNKVAUC_0+PESOGAL*
     >       (LOG(S(I)-XX(1))-FF)*(LOG(S(I)-XX(1))-FF)
          END IF
        END DO
C
        END
C
C******************************************************************************
C Funcion a minimizar: es la suma de un polinomio de primer grado y un 
C perfil de de Vaucouleurs
C B(1)=XX(1), B(2)=XX(2), I0=XX(3), RE=XX(4)
C
        REAL FUNCTION FUNKVAUC_1(XX)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        REAL XX(4)
C
        INTEGER I
        REAL S(NCMAX)
        REAL X0
        REAL FF
        REAL PESOSKY,PESOGAL
        LOGICAL IFFUNK1(NCMAX)                     !region a minimizar de cielo
        LOGICAL IFFUNK2(NCMAX)                   !region a minimizar de galaxia
C
        COMMON/BLKFUNK0/NSCAN
        COMMON/BLKFUNK1/X0
        COMMON/BLKFUNK2/S
        COMMON/BLKFUNK3/IFFUNK1,IFFUNK2
        COMMON/BLKFUNK4/PESOSKY,PESOGAL
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        FUNKVAUC_1=0.
C
        DO I=1,NSCAN
          IF(IFFUNK1(I))THEN
            FF=LOG(XX(3))-7.67*((ABS(REAL(I)-X0)/XX(4))**0.25-1.)
            FUNKVAUC_1=FUNKVAUC_1+PESOSKY*
     >       (LOG(S(I)-XX(1)-XX(2)*REAL(I))-FF)*
     >       (LOG(S(I)-XX(1)-XX(2)*REAL(I))-FF)
          END IF
          IF(IFFUNK2(I))THEN
            FF=LOG(XX(3))-7.67*((ABS(REAL(I)-X0)/XX(4))**0.25-1.)
            FUNKVAUC_1=FUNKVAUC_1+PESOGAL*
     >       (LOG(S(I)-XX(1)-XX(2)*REAL(I))-FF)*
     >       (LOG(S(I)-XX(1)-XX(2)*REAL(I))-FF)
          END IF
        END DO
C
        END
C
C******************************************************************************
C Funcion a minimizar: es la suma de un polinomio de primer grado y un 
C perfil de Sersic
C B(1)=XX(1), I0=XX(2), RE=XX(3), BN=XX(4), N=XX(5)
C
        REAL FUNCTION FUNKSERC_0(XX)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        REAL XX(5)
C
        INTEGER I
        REAL S(NCMAX)
        REAL X0
        REAL FF
        REAL PESOSKY,PESOGAL
        LOGICAL IFFUNK1(NCMAX)                     !region a minimizar de cielo
        LOGICAL IFFUNK2(NCMAX)                   !region a minimizar de galaxia
C
        COMMON/BLKFUNK0/NSCAN
        COMMON/BLKFUNK1/X0
        COMMON/BLKFUNK2/S
        COMMON/BLKFUNK3/IFFUNK1,IFFUNK2
        COMMON/BLKFUNK4/PESOSKY,PESOGAL
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        FUNKSERC_0=0.
C
        DO I=1,NSCAN
          IF(IFFUNK1(I))THEN
            FF=LOG(XX(2))-XX(4)*((ABS(REAL(I)-X0)/XX(3))**(1./XX(5))-1.)
            FUNKSERC_0=FUNKSERC_0+PESOSKY*
     >       (LOG(S(I)-XX(1))-FF)*(LOG(S(I)-XX(1))-FF)
          END IF
          IF(IFFUNK2(I))THEN
            FF=LOG(XX(2))-XX(4)*((ABS(REAL(I)-X0)/XX(3))**(1./XX(5))-1.)
            FUNKSERC_0=FUNKSERC_0+PESOGAL*
     >       (LOG(S(I)-XX(1))-FF)*(LOG(S(I)-XX(1))-FF)
          END IF
        END DO
C
        END
C
C******************************************************************************
C Funcion a minimizar: es la suma de un polinomio de primer grado y un 
C perfil de Sersic
C B(1)=XX(1), B(2)=XX(2), I0=XX(3), RE=XX(4), BN=XX(5), N=XX(6)
C
        REAL FUNCTION FUNKSERC_1(XX)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        REAL XX(6)
C
        INTEGER I
        REAL S(NCMAX)
        REAL X0
        REAL FF
        REAL PESOSKY,PESOGAL
        LOGICAL IFFUNK1(NCMAX)                     !region a minimizar de cielo
        LOGICAL IFFUNK2(NCMAX)                   !region a minimizar de galaxia
C
        COMMON/BLKFUNK0/NSCAN
        COMMON/BLKFUNK1/X0
        COMMON/BLKFUNK2/S
        COMMON/BLKFUNK3/IFFUNK1,IFFUNK2
        COMMON/BLKFUNK4/PESOSKY,PESOGAL
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        FUNKSERC_1=0.
C
        DO I=1,NSCAN
          IF(IFFUNK1(I))THEN
            FF=LOG(XX(3))-XX(5)*((ABS(REAL(I)-X0)/XX(4))**(1./XX(5))-1.)
            FUNKSERC_1=FUNKSERC_1+PESOSKY*
     >       (LOG(S(I)-XX(1)-XX(2)*REAL(I))-FF)*
     >       (LOG(S(I)-XX(1)-XX(2)*REAL(I))-FF)
          END IF
          IF(IFFUNK2(I))THEN
            FF=LOG(XX(3))-XX(5)*((ABS(REAL(I)-X0)/XX(4))**(1./XX(5))-1.)
            FUNKSERC_1=FUNKSERC_1+PESOGAL*
     >       (LOG(S(I)-XX(1)-XX(2)*REAL(I))-FF)*
     >       (LOG(S(I)-XX(1)-XX(2)*REAL(I))-FF)
          END IF
        END DO
C
        END
