!------------------------------------------------------------------------------
! Version 20-November-2026                                       file: testwc.f
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
! Program: testwc
! Classification: wavelengths
! Description: Test the wavelength calibration by measuring the positions of 
! sky lines.
!
!Comment
!
! comprueba la calibracion en l.d.o. de una imagen
!
        PROGRAM TESTWC
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER READILIM
        REAL READF
!
        INTEGER NLMAX
        PARAMETER(NLMAX=100)                           !numero maximo de lineas
!
        INTEGER I,J
        INTEGER J1,J2,IC
        INTEGER NB,NBLOCAL
        INTEGER NLINES,NL                        !numero de lineas seleccionado
        INTEGER NF                      !numero de puntos a ajustar a gaussiana
        INTEGER NC1FIT(NLMAX),NC2FIT(NLMAX)
        INTEGER NSKYLINEST,NSKYLINE
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        REAL A(NCMAX,NSMAX),ERR(NCMAX,NSMAX)
        REAL SP(NCMAX)
        REAL XC,YC
        REAL XF(NCMAX),YF(NCMAX)                  !puntos a ajustar a gaussiana
        REAL X0,SIGMA,AMP,Y0                 !parametros de ajuste a gaussianas
        REAL EEX0,EESIGMA,EEAMP,EEY0       !rms from DOWNHILL en los parametros
        REAL YRMSTOL
        REAL WV0
        REAL WV0_EXP(NLMAX)                    !l.d.o. esperada para cada linea
        REAL DWV0_MEAS(NSMAX,NLMAX)   !error l.d.o. medida en cada linea y scan
        REAL MINDWV,WV(NLINMAX)
        CHARACTER*1 CERR,CH,CSURE,CFITOK,CPLOT
        CHARACTER*8 CLABEL(NLINMAX)
        CHARACTER*50 CDUMMY
        CHARACTER*255 INFILE,ERRFILE
        LOGICAL IFSCAN_PLOT(NSMAX)
        LOGICAL LPSKYLINES,LBEXIST
        LOGICAL LCOLOR(MAX_ID_RED)
!
        COMMON/BLKDATA0A/NSCAN,NCHAN
        COMMON/BLKDATA0B/STWV,DISP
        COMMON/BLKDATA1/A
        COMMON/BLKDATA2/SP
        COMMON/BLKDATA3/DWV0_MEAS
        COMMON/BLKDATA4/IFSCAN_PLOT
        COMMON/BLKDATA5/NLINES
        COMMON/BLKFITG1/NF
        COMMON/BLKFITG2/XF,YF
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
!------------------------------------------------------------------------------
        OUTFILEX=OUTFILEX
!
        THISPROGRAM='testwc'
        CALL WELCOME('20-November-2026')
!------------------------------------------------------------------------------
        NLINES=0
        LPSKYLINES=.FALSE.
!------------------------------------------------------------------------------
        WRITE(*,100,ADVANCE='NO')'Work with error images (y/n) '
        CERR(1:1)=READC('n','yn')
!
        CALL RPGBEGIN(NTERM,IDN,LCOLOR)
!------------------------------------------------------------------------------
        WRITE(*,100,ADVANCE='NO')'Input file name'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        IF((STWV.EQ.0).OR.(DISP.EQ.0.))THEN
          WRITE(*,101)'FATAL ERROR: input file must be '//&
     &     'wavelength calibrated.'
          CLOSE(20)
          STOP
        END IF
        DO I=1,NSCAN
          READ(20) (A(J,I),J=1,NCHAN)
        END DO
        CLOSE(20)
        IF(CERR.EQ.'y')THEN
          CALL GUESSEF(INFILE,ERRFILE)
          ERRFILE=INFILEX(21,ERRFILE,NSCAN,NCHAN,STWV,DISP,21,.TRUE.) !...match
          DO I=1,NSCAN
            READ(21) (ERR(J,I),J=1,NCHAN)
          END DO
          CLOSE(21)
        END IF
!------------------------------------------------------------------------------
        WRITE(*,101)'* Define scan region to obtain averaged spectrum:'
        CALL DEFINESCAN
!------------------------------------------------------------------------------
        CALL BUTTON(1,'[z]oom',0)
        CALL BUTTON(2,'[w]hole',0)
        CALL BUTTON(3,'sky [l]ines',0)
        CALL BUTTON(4,'[s]cans',0)
        CALL BUTTON(7,'[f]it line',0)
        CALL BUTTON(8,'[g]o',0)
        CALL BUTTON(8,'[g]o',3)
        CALL BUTTON(12,'E[x]it',0)
!------------------------------------------------------------------------------
        CALL PLOTSP(1,NCHAN)
!------------------------------------------------------------------------------
10      CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
        CALL IFBUTTON(XC,YC,NB)
!
        CALL CHLOWER(CH)
        NBLOCAL=INDEX('zwls  fg   x',CH)
        IF((NBLOCAL.NE.0).AND.(CH.NE.' '))THEN
          CALL BUTTQEX(NBLOCAL,LBEXIST)
          IF(LBEXIST) NB=NBLOCAL
        END IF
!..............................................................................
11      IF(NB.EQ.0)THEN
!..............................................................................
        ELSEIF(NB.EQ.1)THEN
          CALL BUTTON(1,'[z]oom',5)
          WRITE(*,100,ADVANCE='NO')'Point #1: press mouse button...'
          IF(LCOLOR(1)) CALL PGSCI(5)
          CALL RPGBAND(6,0,0.,0.,XC,YC,CH)
          IF(LCOLOR(1)) CALL PGSCI(1)
          WRITE(*,101)' OK!'
          J1=NINT(XC)
          IF(J1.LT.1) J1=1
          IF(J1.GT.NCHAN) J1=NCHAN
          WRITE(*,100,ADVANCE='NO')'Point #2: press mouse button...'
          IF(LCOLOR(1)) CALL PGSCI(5)
          CALL RPGBAND(4,0,REAL(J1),0.,XC,YC,CH)
          IF(LCOLOR(1)) CALL PGSCI(1)
          WRITE(*,101)' OK!'
          J2=NINT(XC)
          IF(J2.LT.1) J2=1
          IF(J2.GT.NCHAN) J2=NCHAN
          IF(J2.LT.J1)THEN
            IC=J1
            J1=J2
            J2=IC
          END IF
          CALL PLOTSP(J1,J2)
          IF(LPSKYLINES) CALL PLOTSKYLINES
          CALL BUTTON(1,'[z]oom',0)
!..............................................................................
        ELSEIF(NB.EQ.2)THEN
          CALL BUTTON(2,'[w]hole',5)
          CALL PLOTSP(1,NCHAN)
          IF(LPSKYLINES) CALL PLOTSKYLINES
          IF(NLINES.GT.0)THEN
            DO NL=1,NLINES
              NF=0
              DO J=NC1FIT(NL),NC2FIT(NL)
                NF=NF+1
                XF(NF)=REAL(J)
                YF(NF)=SP(J)
              END DO
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                IF(LCOLOR(ITERM)) CALL PGSCI(2)
                CALL PGBIN(NF,XF,YF,.TRUE.)
                IF(LCOLOR(ITERM)) CALL PGSCI(1)
              END DO
              CALL GAUSCFIT(X0,SIGMA,AMP,Y0,EEX0,EESIGMA,EEAMP,EEY0,&
     &         YRMSTOL)
              NF=0
              DO J=NC1FIT(NL),NC2FIT(NL)
                NF=NF+1
                XF(NF)=REAL(J)
                YF(NF)=Y0+AMP*EXP(-(REAL(J)-X0)*(REAL(J)-X0)/&
     &           (2.*SIGMA*SIGMA))
              END DO
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                IF(LCOLOR(ITERM)) CALL PGSCI(3)
                CALL PGLINE(NF,XF,YF)
                IF(LCOLOR(ITERM)) CALL PGSCI(1)
              END DO
            END DO
          END IF
          CALL BUTTON(2,'[w]hole',0)
!..............................................................................
        ELSEIF(NB.EQ.3)THEN
          IF(LPSKYLINES)THEN
            CALL BUTTON(3,'sky [l]ines',0)
            LPSKYLINES=.FALSE.
          ELSE
            CALL BUTTON(3,'sky [l]ines',1)
            LPSKYLINES=.TRUE.
            CALL PLOTSKYLINES
          END IF
!..............................................................................
        ELSEIF(NB.EQ.4)THEN
          CALL BUTTON(4,'[s]cans',5)
          WRITE(*,101)'* Define scan region to obtain averaged'//&
     &     ' spectrum:'
          CALL DEFINESCAN
          CALL BUTTON(4,'[s]cans',0)
          NB=2
          GOTO 11
!..............................................................................
        ELSEIF(NB.EQ.7)THEN
          CALL BUTTON(7,'[f]it line',5)
          WRITE(*,*)
          WRITE(*,101)'* Enter region to be fitted:'
          WRITE(*,100,ADVANCE='NO')'Point #1: press mouse button...'
          IF(LCOLOR(1)) CALL PGSCI(5)
          CALL RPGBAND(6,0,0.,0.,XC,YC,CH)
          IF(LCOLOR(1)) CALL PGSCI(1)
          WRITE(*,101)' OK!'
          J1=NINT(XC)
          IF(J1.LT.1) J1=1
          IF(J1.GT.NCHAN) J1=NCHAN
          WRITE(*,100,ADVANCE='NO')'Point #2: press mouse button...'
          IF(LCOLOR(1)) CALL PGSCI(5)
          CALL RPGBAND(6,0,REAL(J1),0.,XC,YC,CH)
          IF(LCOLOR(1)) CALL PGSCI(1)
          WRITE(*,101)' OK!'
          J2=NINT(XC)
          IF(J2.LT.1) J2=1
          IF(J2.GT.NCHAN) J2=NCHAN
          IF(J2.LT.J1)THEN
            IC=J1
            J1=J2
            J2=IC
          END IF
          NC1FIT(NLINES+1)=J1
          NC2FIT(NLINES+1)=J2
!
          NF=0
          DO J=NC1FIT(NLINES+1),NC2FIT(NLINES+1)
            NF=NF+1
            XF(NF)=REAL(J)
            YF(NF)=SP(J)
          END DO
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(2)
            CALL PGBIN(NF,XF,YF,.TRUE.)
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END DO
!
          CALL GAUSCFIT(X0,SIGMA,AMP,Y0,EEX0,EESIGMA,EEAMP,EEY0,&
     &     YRMSTOL)
          NF=0
          DO J=NC1FIT(NLINES+1),NC2FIT(NLINES+1)
            NF=NF+1
            XF(NF)=REAL(J)
            YF(NF)=Y0+AMP*EXP(-(REAL(J)-X0)*(REAL(J)-X0)/&
     &       (2.*SIGMA*SIGMA))
          END DO
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(3)
            CALL PGLINE(NF,XF,YF)
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END DO
!
          WRITE(*,100,ADVANCE='NO')'Is this fit OK (y/n) '
          CFITOK(1:1)=READC('y','yn')
          IF(CFITOK.EQ.'y')THEN
            NLINES=NLINES+1
            IF(NLINES.EQ.1) CALL BUTTON(8,'[g]o',0)
            WRITE(*,100,ADVANCE='NO')'Measured central wavelength: '
            WV0=(X0-1.)*DISP+STWV
            WRITE(*,*) WV0
            CALL SELLINES(2,NSKYLINEST,WV,CLABEL)           !typical sky lines
            MINDWV=ABS(WV0-WV(1))
            WV0_EXP(NLINES)=WV(1)
            DO NSKYLINE=1,NSKYLINEST
              IF(ABS(WV0-WV(NSKYLINE)).LT.MINDWV)THEN
                MINDWV=ABS(WV0-WV(NSKYLINE))
                WV0_EXP(NLINES)=WV(NSKYLINE)
              END IF
            END DO
            WRITE(*,100,ADVANCE='NO')&
     &         'Enter expected central wavelength '
            WRITE(CDUMMY,*)WV0_EXP(NLINES)
            WV0_EXP(NLINES)=READF(CDUMMY)
          END IF
!
          CALL BUTTON(7,'[f]it line',0)
!..............................................................................
        ELSEIF(NB.EQ.8)THEN
          CALL BUTTON(8,'[g]o',5)
!
          WRITE(*,110)'>>> Current number of defined lines: ',NLINES
          WRITE(*,101)'Fitting lines (please wait)...'
          DO I=1,NSCAN
            WRITE(*,'(A,I4.4)',ADVANCE='NO')'\b\b\b\b',I
            IF(IFSCAN_PLOT(I))THEN
              DO NL=1,NLINES
                NF=0
                DO J=NC1FIT(NL),NC2FIT(NL)
                  NF=NF+1
                  XF(NF)=REAL(J)
                  YF(NF)=A(J,I)
                END DO
                CALL GAUSCFIT(X0,SIGMA,AMP,Y0,EEX0,EESIGMA,EEAMP,EEY0,&
     &           YRMSTOL)
                DWV0_MEAS(I,NL)=(X0-1.)*DISP+STWV-WV0_EXP(NL)
              END DO
            END IF
          END DO
          WRITE(*,*)
          WRITE(*,101)'  OK!'
!
          CPLOT='y'
          DO WHILE(CPLOT.EQ.'y')
            WRITE(*,100,ADVANCE='NO')&
     &         'Line number to be plotted (-1=ALL,0=EXIT) '
            NL=READILIM('0',-1,NLINES)
            IF(NL.EQ.0)THEN
              CPLOT='n'
            ELSE
              CALL RPGERASW(0.,1.,0.,0.80)
              CALL PLOTLINE(NL)
            END IF
          END DO
!
          CALL BUTTON(8,'[g]o',0)
          NB=2
          GOTO 11
!..............................................................................
        ELSEIF(NB.EQ.12)THEN
          CALL BUTTON(12,'E[x]it',5)
          WRITE(*,100,ADVANCE='NO')'Do you really want to Exit (y/n) '
          CSURE(1:1)=READC('n','yn')
          IF(CSURE.EQ.'y') GOTO 90
          CALL BUTTON(12,'E[x]it',0)
!..............................................................................
        END IF
        GOTO 10
!------------------------------------------------------------------------------
90      CALL PGEND
!
        STOP
100     FORMAT(A)
101     FORMAT(A)
110     FORMAT(A,I6)
        END
!
!******************************************************************************
! Dibuja el espectro promedio entre NS1 y NS2, usando como limites del plot
! (en canales) NC1 y NC2
        SUBROUTINE PLOTSP(NC1,NC2)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INTEGER NC1,NC2
!
        INTEGER I,J
        INTEGER NSADD
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        REAL A(NCMAX,NSMAX)
        REAL SP(NCMAX)
        REAL XP(NCMAX)
        REAL XMIN,XMAX,YMIN,YMAX,DY
        REAL XMINL,XMAXL
        LOGICAL IFSCAN_PLOT(NSMAX)
        LOGICAL LCOLOR(MAX_ID_RED)
!
        COMMON/BLKDATA0A/NSCAN,NCHAN
        COMMON/BLKDATA0B/STWV,DISP
        COMMON/BLKDATA1/A
        COMMON/BLKDATA2/SP
        COMMON/BLKDATA4/IFSCAN_PLOT
        COMMON/BLKPLIMITS/XMIN,XMAX,YMIN,YMAX
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
!------------------------------------------------------------------------------
        DO J=1,NCHAN
          XP(J)=REAL(J)
          SP(J)=0.
        END DO
        NSADD=0
        DO I=1,NSCAN
          IF(IFSCAN_PLOT(I))THEN
            NSADD=NSADD+1
            DO J=NC1,NC2
              SP(J)=SP(J)+A(J,I)
            END DO
          END IF
        END DO
!
        DO J=NC1,NC2
          SP(J)=SP(J)/REAL(NSADD)
        END DO
!
        XMIN=REAL(NC1)
        XMAX=REAL(NC2)
        YMIN=SP(NC1)
        YMAX=YMIN
        DO J=NC1+1,NC2
          IF(YMIN.GT.SP(J)) YMIN=SP(J)
          IF(YMAX.LT.SP(J)) YMAX=SP(J)
        END DO
        DY=YMAX-YMIN
        YMIN=YMIN-DY/30.
        YMAX=YMAX+DY/30.
!
        XMINL=(XMIN-1.)*DISP+STWV
        XMAXL=(XMAX-1.)*DISP+STWV
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(ITERM.EQ.1)THEN
            CALL RPGERASW(0.,1.,0.,0.80)
            CALL RPGENV(XMIN,XMAX,YMIN,YMAX,0,-2)
          ELSE
            CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,-2)
          END IF
          CALL PGSWIN(XMINL,XMAXL,YMIN,YMAX)
          CALL PGBOX('CMTS',0.,0,'BCNST',0.,0)
          CALL PGSWIN(XMIN,XMAX,YMIN,YMAX)
          CALL PGBOX('BNTS',0.,0,char(32),0.,0)
          CALL PGIDEN_RED
          CALL PGBIN(NCHAN,XP,SP,.TRUE.)
          CALL PGLABEL('channel','No. counts',CHAR(32))
        END DO
!
        END
!
!******************************************************************************
! Dibuja la diferencia entre la l.d.o. esperada y la l.d.o. calculada en los
! ajustes de la linea NL
        SUBROUTINE PLOTLINE(NL)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INTEGER NL
!
        INTEGER NLMAX
        PARAMETER(NLMAX=100)                           !numero maximo de lineas
!
        INTEGER NPLOT,NLINES
        INTEGER NL0,NL1,NL2
        INTEGER I,L,NCOLOR,K
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        REAL XMIN,XMAX,YMIN,YMAX,DY
        REAL XPLOT(NSMAX),YPLOT(NSMAX)
        REAL DWV0_MEAS(NSMAX,NLMAX)
        DOUBLE PRECISION MEAN(NLMAX),SIGMA(NLMAX)
        CHARACTER*50 CDUMMY
        LOGICAL IFSCAN_PLOT(NSMAX)
        LOGICAL LCOLOR(MAX_ID_RED)
!
        COMMON/BLKDATA0A/NSCAN,NCHAN
        COMMON/BLKDATA3/DWV0_MEAS
        COMMON/BLKDATA4/IFSCAN_PLOT
        COMMON/BLKDATA5/NLINES
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
!------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
!
        NCOLOR=0
!
        IF(NL.EQ.-1)THEN
          NL1=1
          NL2=NLINES
        ELSE
          NL1=NL
          NL2=NL
        END IF
!
        XMIN=1.
        XMAX=REAL(NSCAN)
        YMIN=+1.E20
        YMAX=-1.E20
        DO NL0=NL1,NL2
          NPLOT=0
          DO I=1,NSCAN
            IF(IFSCAN_PLOT(I))THEN
              NPLOT=NPLOT+1
              YPLOT(NPLOT)=DWV0_MEAS(I,NL0)
            END IF
          END DO
          DO I=1,NPLOT
            IF(YMIN.GT.YPLOT(I)) YMIN=YPLOT(I)
            IF(YMAX.LT.YPLOT(I)) YMAX=YPLOT(I)
          END DO
        END DO
        DY=YMAX-YMIN
        YMIN=YMIN-DY/30.
        YMAX=YMAX+DY/30.
!
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(ITERM.EQ.1)THEN
            CALL RPGERASW(0.,1.,0.,0.80)
            CALL RPGENV(XMIN,XMAX,YMIN,YMAX,0,0)
          ELSE
            CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,0)
          END IF
          IF(NL.EQ.-1)THEN
            CALL PGLABEL('scan','walength difference',&
     &       'all lines')
          ELSE
            WRITE(CDUMMY,*)NL
            CALL RMBLANK(CDUMMY,CDUMMY,L)
            CALL PGLABEL('scan','walength difference',&
     &       'Line #'//CDUMMY(1:L))
          END IF
        END DO
!
        DO NL0=NL1,NL2
!..............................................................................
! dibujamos
          NCOLOR=NCOLOR+1
          IF(NCOLOR.GT.12) NCOLOR=1
          NPLOT=0
          DO I=1,NSCAN
            IF(IFSCAN_PLOT(I))THEN
              NPLOT=NPLOT+1
              XPLOT(NPLOT)=REAL(I)
              YPLOT(NPLOT)=DWV0_MEAS(I,NL0)
            END IF
          END DO
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(NCOLOR)
            CALL PGPOINT(NPLOT,XPLOT,YPLOT,17)
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END DO
!..............................................................................
! estadistica
          MEAN(NL0)=0.D0
          DO K=1,NPLOT
            MEAN(NL0)=MEAN(NL0)+DBLE(YPLOT(K))
          END DO
          MEAN(NL0)=MEAN(NL0)/DBLE(NPLOT)
          SIGMA(NL0)=0.D0
          DO K=1,NPLOT
            SIGMA(NL0)=SIGMA(NL0)+&
     &       (DBLE(YPLOT(K))-MEAN(NL0))*(DBLE(YPLOT(K))-MEAN(NL0))
          END DO
          SIGMA(NL0)=DSQRT(SIGMA(NL0)/DBLE(NPLOT-1))
          WRITE(*,'(A,I2.2)')'* Sky line #',NL0
          WRITE(*,100,ADVANCE='NO')'> Mean.....: '
          WRITE(*,*) REAL(MEAN(NL0))
          WRITE(*,100,ADVANCE='NO')'> Sigma....: '
          WRITE(*,*) REAL(SIGMA(NL0))
!..............................................................................
        END DO
!
100     FORMAT(A)
        END
!
!******************************************************************************
! Dibuja las posiciones de las lineas tipicas del cielo nocturno
        SUBROUTINE PLOTSKYLINES
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
!
        INTEGER NLINEST,NLINE
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        REAL XMIN,XMAX,YMIN,YMAX,DY
        REAL WV(NLINMAX)
        REAL X0
        CHARACTER*8 CLABEL(NLINMAX)
        LOGICAL LCOLOR(MAX_ID_RED)
!
        COMMON/BLKDATA0B/STWV,DISP
        COMMON/BLKPLIMITS/XMIN,XMAX,YMIN,YMAX
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
!------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
!
        DY=YMAX-YMIN
        CALL SELLINES(2,NLINEST,WV,CLABEL)   !lineas tipicas del cielo nocturno 
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM)) CALL PGSCI(5)
          DO NLINE=1,NLINEST
            X0=(WV(NLINE)-STWV)/DISP+1.
            IF((X0.GE.XMIN).AND.(X0.LE.XMAX))THEN
              CALL PGSLS(2)
              CALL PGMOVE(X0,YMAX)
              CALL PGDRAW(X0,YMIN)
              CALL PGSLS(1)
              CALL PGSCH(0.8)
              CALL PGPTEXT(X0,YMAX+DY/100.,0.,.5,CLABEL(NLINE))
              CALL PGSCH(1.0)
            END IF
          END DO
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
        END DO
!
        END
!
!******************************************************************************
! Definimos scans a utilizar en la variable IFSCAN
        SUBROUTINE DEFINESCAN
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
!
        INTEGER I,NS0,NS1,NS2
        LOGICAL IFSCAN_PLOT(NSMAX)
!
        COMMON/BLKDATA0A/NSCAN,NCHAN
        COMMON/BLKDATA4/IFSCAN_PLOT
!------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
!
        DO I=1,NSCAN
          IFSCAN_PLOT(I)=.FALSE.
        END DO
10      WRITE(*,100,ADVANCE='NO')'Scan region (0,0=EXIT) '
        CALL READ2I('0,0',NS1,NS2)
        IF((NS1.EQ.0).AND.(NS2.EQ.0)) GOTO 20
        IF((NS1.LT.1).OR.(NS2.GT.NSCAN).OR.(NS1.GT.NS2))THEN
          WRITE(*,101)'ERROR: numbers out of range. Try again.'
          GOTO 10
        END IF
        DO I=NS1,NS2
          IFSCAN_PLOT(I)=.TRUE.
        END DO
        GOTO 10
20      NS0=0
        DO I=1,NSCAN
          IF(IFSCAN_PLOT(I)) NS0=NS0+1
        END DO
        IF(NS0.EQ.0)THEN
          WRITE(*,101)'ERROR: number of defined scans = 0. Try again.'
          GOTO 10
        END IF
!
100     FORMAT(A)
101     FORMAT(A)
        END
