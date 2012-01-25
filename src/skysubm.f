C------------------------------------------------------------------------------
C Version 07-September-2007                                    file: skysubm.f
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
C Program: skysubm
C Classification: sky subtraction
C Description: Calculates and subtracts a sky image by fitting polynomials to 
C each channel.
C
Comment
C
C Calcula el cielo en una imagen ajustando polinomios a cada canal. El
C programa resta el cielo a la imagen inicial y permite salvar tambien
C el cielo calculado. Trabaja con las imagenes de errores si se solicita.
C
C En la version modificada en Lick he cambiado POLFIT por POLFITSIG, para
C poder asi hacer ajuste eliminando puntos malos (cosmic rays,...)
C
        PROGRAM SKYSUBM
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER TRUELEN
        INTEGER READI
        INTEGER READILIM
        REAL READF
C
        INTEGER NMAXSIMUL
        PARAMETER (NMAXSIMUL=1000)     !numero maximo simulaciones para errores
        REAL PI2
        PARAMETER(PI2=6.283185307)
C
        REAL RANRED
C
        INTEGER I,J,M,L
        INTEGER I1,I2
        INTEGER NI,NI0,NC1,NC2
        INTEGER NPTS,NTERMS
        INTEGER NSKYREG,L1(20),L2(20),L0
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        INTEGER NSIMUL,NSIMULT,NSEED
        REAL G(NCMAX,NSMAX),SKY(NCMAX,NSMAX)
        REAL ERR(NCMAX,NSMAX),ESKY(NCMAX,NSMAX),ESKYG(NCMAX,NSMAX)
        REAL ERRSKY,ERRSKYG,SKYSIGMA1(NSMAX),SKYSIGMA2(NSMAX)
        REAL SKYSIGMA1G(NSMAX),SKYSIGMA2G(NSMAX)
        REAL XMIN,XMAX,YMIN,YMAX,DX,DY
        REAL X(NSMAX),Y(NSMAX)
        REAL XT(NSMAX),YT(NSMAX)
        REAL A(20)
!       REAL CHISQR
        REAL XP,POL,SK(NSMAX),SKSIGMA1(NSMAX),SKSIGMA2(NSMAX)
        REAL SKYRAN(NSMAX,NMAXSIMUL)
        REAL XC,YC
        REAL R1,R2,GXX
        REAL MEAN,SIGMA
        REAL TSIGMA
        CHARACTER*1 CLIMITS,CPLOT,CSAVSK,CERR,CERR2,CERRG,CMOUSE,CH
        CHARACTER*50 CDUMMY
        CHARACTER*80 GLABEL
        CHARACTER*75 INFILE,SKYFILE,OUTFILE,ERRFILE
        LOGICAL FIRSTPLOT,LLOOP
        LOGICAL LCOLOR(MAX_ID_RED)
        LOGICAL IFSCAN(NSMAX)
        LOGICAL LFIT2(NSMAX)
C
        COMMON/BLKSKY/NSKYREG,L1,L2
        COMMON/BLKLIM/XMIN,XMAX
C------------------------------------------------------------------------------
        THISPROGRAM='skysubm'
        CALL WELCOME('07-September-2007')
C
        NSIMULT=100
        TSIGMA=3.0
        ERRSKY=0.0 !evita un warning de compilacion
C------------------------------------------------------------------------------
        WRITE(*,100) 'Work with error images (y/n) '
        CERR(1:1)=READC('y','yn')
        WRITE(*,100) 'Do you want to estimate errors from the'
        WRITE(*,100) ' r.m.s. in the sky fits (y/n) '
        CERRG(1:1)=READC('y','yn')
C
        FIRSTPLOT=.TRUE.
        NTERMS=1                              !por defecto polinomio de grado 0
        CALL PIDEGTER(NTERM,IDN,LCOLOR)
C
        WRITE(*,100) 'Input file name......'
        INFILE=INFILEX(13,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        IF(CERR.EQ.'y') CALL GUESSEF(INFILE,ERRFILE)
        IF(TRUELEN(OBJECT).GT.0)THEN
          INFILE=INFILE(1:TRUELEN(INFILE))//' ['//
     +     OBJECT(1:TRUELEN(OBJECT))//']'
        END IF
        DO I=1,NSCAN
          READ(13) (G(J,I),J=1,NCHAN)
        END DO
        CLOSE(13)
        IF(CERR.EQ.'y')THEN
          WRITE(*,100) 'Input error file name '
          ERRFILE=INFILEX(14,ERRFILE,NSCAN,NCHAN,STWV,DISP,21,.TRUE.) !...match
          DO I=1,NSCAN
            READ(14) (ERR(J,I),J=1,NCHAN)
          END DO
          CLOSE(14)
        END IF
C
C calculamos corte espacial promedio
10      WRITE(*,100) '* Enter first and last channel to compute'
        WRITE(*,101) ' averaged spatial profile:'
        WRITE(*,100) 'Channels '
        WRITE(CDUMMY,'(A2,I10)') '1,',NCHAN
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        CALL READ2I(CDUMMY(1:L),NC1,NC2)
        IF((NC1.LT.1).OR.(NC2.GT.NCHAN).OR.(NC1.GT.NC2))THEN
          WRITE(*,101) 'ERROR: invalid entry. Try again.'
          GOTO 10
        END IF
        DO I=1,NSCAN
          X(I)=REAL(I)
          Y(I)=0.
          DO J=NC1,NC2
            Y(I)=Y(I)+G(J,I)
          END DO
          Y(I)=Y(I)/REAL(NC2-NC1+1)
        END DO
C
C dibujamos corte espacial promedio
        XMIN=1.
        XMAX=REAL(NSCAN)
        DX=XMAX-XMIN
        XMIN=XMIN-DX/50.
        XMAX=XMAX+DX/50.
C
20      YMIN=Y(1)
        YMAX=YMIN
        DO I=2,NSCAN
          IF(Y(I).LT.YMIN)YMIN=Y(I)
          IF(Y(I).GT.YMAX)YMAX=Y(I)
        END DO
        DY=YMAX-YMIN
        YMIN=YMIN-DY/50.
        YMAX=YMAX+DY/50.
C
21      DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,-2)
          CALL PGBOX('BCTNSI',0.0,0,'BCTNSI',0.0,0)
          CALL PGIDEN_RED
          CALL PGLABEL('Scan','No. counts (averaged)',INFILE)
          CALL PGBIN(NSCAN,X,Y,.TRUE.)
        END DO
        WRITE(*,*)
        WRITE(*,101) '(1) Change x/y-limits'
        WRITE(*,101) '(2) Compute new spatial profile'
        WRITE(*,101) '(0) continue'
        WRITE(*,100) 'Option (0/1/2) '
        CLIMITS(1:1)=READC('0','012')
C
        IF(CLIMITS.EQ.'1')THEN
          WRITE(*,100) 'XMIN, XMAX (0,0=automatic) '
          CALL READ2I('0,0',I1,I2)
          IF((I1.EQ.0).AND.(I2.EQ.0))THEN
            XMIN=1.
            XMAX=REAL(NSCAN)
            DX=XMAX-XMIN
            XMIN=XMIN-DX/50.
            XMAX=XMAX+DX/50.
          ELSE
            XMIN=REAL(I1)
            XMAX=REAL(I2)
          END IF
          WRITE(*,100) 'YMIN, YMAX (0,0=automatic) '
          CALL READ2I('0,0',I1,I2)
          IF((I1.EQ.0).AND.(I2.EQ.0)) GOTO 20
          YMIN=REAL(I1)
          YMAX=REAL(I2)
          GOTO 21
        ELSEIF(CLIMITS.EQ.'2')THEN
          GOTO 10
        END IF
C------------------------------------------------------------------------------
        WRITE(*,100) 'Select sky regions with [m]ouse or [k]eyboard '//
     +   '(m/k) '
        CMOUSE(1:1)=READC('m','mk')
C------------------------------------------------------------------------------
        IF(CMOUSE.EQ.'k')THEN               !seleccionamos regiones con teclado
          WRITE(*,100) 'No. of sky regions to be employed'
          NSKYREG=READI('@')
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(3)
          END DO
          DO I=1,NSKYREG
            WRITE(*,*)
28          WRITE(*,'(A,I2)') 'Region #',I
            WRITE(*,100) 'left limit, right limit'
            CALL READ2I('@',L1(I),L2(I))
            IF((L1(I).LT.1).OR.(L2(I).GT.NSCAN).OR.(L2(I).LT.L1(I)))THEN
              WRITE(*,101) 'ERROR: numbers out of range. Try again.'
              GOTO 28
            END IF
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              CALL PGMOVE(REAL(L1(I)),YMIN)
              CALL PGDRAW(REAL(L1(I)),YMAX)
              CALL PGMOVE(REAL(L2(I)),YMIN)
              CALL PGDRAW(REAL(L2(I)),YMAX)
              CALL PGRECT(REAL(L1(I)),REAL(L2(I)),
     +         YMIN,YMIN+(YMAX-YMIN)/80.)
              CALL PGRECT(REAL(L1(I)),REAL(L2(I)),
     +         YMAX,YMAX-(YMAX-YMIN)/80.)
            END DO
          END DO
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END DO
        ELSE                                  !seleccionamos regiones con raton
          LLOOP=.TRUE.
          NSKYREG=0
          DO WHILE(LLOOP)
            WRITE(*,101)'(Press <q>, <x> or right mouse button to EXIT)'
            WRITE(*,'(A,I2)') '>>> Next region is #',NSKYREG+1
            WRITE(*,100) 'Press mouse button...'
            IF(LCOLOR(1)) CALL PGSCI(2)
            CALL PGBAND(6,0,0.,0.,XC,YC,CH)
            IF((CH.EQ.'q').OR.(CH.EQ.'Q').OR.
     +       (CH.EQ.'x').OR.(CH.EQ.'X'))THEN
              LLOOP=.FALSE.
            END IF
            IF(LLOOP)THEN
              L1(NSKYREG+1)=NINT(XC)
              IF(L1(NSKYREG+1).LT.1) L1(NSKYREG+1)=1
              IF(L1(NSKYREG+1).GT.NSCAN) L1(NSKYREG+1)=NSCAN
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                IF(LCOLOR(ITERM)) CALL PGSCI(3)
                CALL PGMOVE(REAL(L1(NSKYREG+1)),YMIN)
                CALL PGDRAW(REAL(L1(NSKYREG+1)),YMAX)
              END DO
              WRITE(*,110) '   point #1, scan ',L1(NSKYREG+1)
              WRITE(*,100) 'Press mouse button...'
              IF(LCOLOR(1)) CALL PGSCI(2)
              CALL PGBAND(6,0,0.,0.,XC,YC,CH)
              IF((CH.EQ.'q').OR.(CH.EQ.'Q').OR.
     +         (CH.EQ.'x').OR.(CH.EQ.'X'))THEN
                LLOOP=.FALSE.
              END IF
              IF(LLOOP)THEN
                L2(NSKYREG+1)=NINT(XC)
                IF(L2(NSKYREG+1).LT.1) L2(NSKYREG+1)=1
                IF(L2(NSKYREG+1).GT.NSCAN) L2(NSKYREG+1)=NSCAN
                DO ITERM=NTERM,1,-1
                  CALL PGSLCT(IDN(ITERM))
                  IF(LCOLOR(ITERM)) CALL PGSCI(3)
                  CALL PGMOVE(REAL(L2(NSKYREG+1)),YMIN)
                  CALL PGDRAW(REAL(L2(NSKYREG+1)),YMAX)
                END DO
                WRITE(*,110) '   point #2, scan ',L2(NSKYREG+1)
                IF(L2(NSKYREG+1).LT.L1(NSKYREG+1))THEN
                  L0=L1(NSKYREG+1)
                  L1(NSKYREG+1)=L2(NSKYREG+1)
                  L2(NSKYREG+1)=L0
                END IF
                NSKYREG=NSKYREG+1
                CALL PGRECT(REAL(L1(NSKYREG)),REAL(L2(NSKYREG)),
     +           YMIN,YMIN+(YMAX-YMIN)/80.)
                CALL PGRECT(REAL(L1(NSKYREG)),REAL(L2(NSKYREG)),
     +           YMAX,YMAX-(YMAX-YMIN)/80.)
              END IF
            END IF
            WRITE(*,*)
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              IF(LCOLOR(ITERM)) CALL PGSCI(1)
            END DO
          END DO
        END IF
C declaramos una variable logica para evitar el problema del solape de las
C regiones a utilizar como cielo
        DO J=1,NSCAN
          IFSCAN(J)=.FALSE.
        END DO
        DO I=1,NSKYREG
          DO J=L1(I),L2(I)
            IFSCAN(J)=.TRUE.
          END DO
        END DO
C------------------------------------------------------------------------------
        WRITE(*,*)
        WRITE(*,100) 'Plot individual fits (y/n) '
        CPLOT(1:1)=READC('y','yn')
        IF(CPLOT.EQ.'n')THEN
          CALL PGEND
          GOTO 130
        END IF
C
        IF(CERR.EQ.'y')THEN
          NSEED=-1                      !inicializaremos RANRED la primera vez
          WRITE(*,100) 'No. of simulations to compute errors '
          WRITE(CDUMMY,*) NSIMULT
          NSIMULT=READILIM(CDUMMY,10,NMAXSIMUL)
        END IF
C------------------------------------------------------------------------------
        WRITE(*,100) 'Polynomial degree '
        NTERMS=READILIM('@',0,19)
        NTERMS=NTERMS+1
        WRITE(*,100) 'Times sigma to remove points '
        TSIGMA=READF('3.0')
C
117     CONTINUE
        WRITE(*,*)
        WRITE(*,101) 'ENTER CHANNEL NUMBER TO BE FITTED'
        WRITE(*,101) ' N=Channel number (1.LE.N.LE.NCHAN)'
        WRITE(*,101) ' 0=exit'
        IF(.NOT.FIRSTPLOT)THEN
          WRITE(*,101) '-1=change pol. degree'
          WRITE(*,101) '-2=change sky regions'
          WRITE(*,101) '-3=change limits of previous plot'
          IF(CERR.EQ.'y') WRITE(*,101) '-4=change no. of simulations'
        END IF
        WRITE(*,100) 'Channel'
        NI0=READI('@')
        IF(NI0.GT.NCHAN)THEN
          WRITE(*,101) 'ERROR: Invalid entry. Try again.'
          GOTO 117
        END IF
        IF(CERR.EQ.'y')THEN
          IF(NI0.LT.-4)THEN
            WRITE(*,101) 'ERROR: Invalid entry. Try again.'
            GOTO 117
          END IF
        ELSE
          IF(NI0.LT.-3)THEN
            WRITE(*,101) 'ERROR: Invalid entry. Try again.'
            GOTO 117
          END IF
        END IF
C
        IF(FIRSTPLOT)THEN
          IF(NI0.LT.0)THEN
            WRITE(*,101) 'ERROR: Invalid entry. Try again.'
            GOTO 117
          END IF
        ELSE
          IF(NI0.EQ.-3)THEN
            WRITE(*,100) 'XMIN, XMAX (0,0=automatic) '
            CALL READ2I('0,0',I1,I2)
            IF((I1.EQ.0).AND.(I2.EQ.0))THEN
              XMIN=1.
              XMAX=REAL(NSCAN)
              DX=XMAX-XMIN
              XMIN=XMIN-DX/50.
              XMAX=XMAX+DX/50.
            ELSE
              XMIN=REAL(I1)
              XMAX=REAL(I2)
            END IF
            WRITE(*,100) 'YMIN, YMAX (0,0=automatic) '
            CALL READ2I('0,0',I1,I2)
            IF((I1.EQ.0).AND.(I2.EQ.0))THEN
              YMIN=YT(1)
              YMAX=YMIN
              DO I=2,NSCAN
                IF(YT(I).LT.YMIN)YMIN=YT(I)
                IF(YT(I).GT.YMAX)YMAX=YT(I)
              END DO
              DY=YMAX-YMIN
              YMIN=YMIN-DY/50.
              YMAX=YMAX+DY/50.
            ELSE
              YMIN=REAL(I1)
              YMAX=REAL(I2)
            END IF
            GOTO 30
          END IF
        END IF
        NI=NI0
C
        IF(FIRSTPLOT) FIRSTPLOT=.FALSE.
C
        IF(NI.EQ.0)GOTO 120
        IF(NI.EQ.-1)THEN
          WRITE(*,100) 'Polynomial degree '
          NTERMS=READILIM('@',0,19)
          NTERMS=NTERMS+1
          GOTO 117
        END IF
C
        IF(NI.EQ.-2)THEN
          CALL CHANGESKY
          DO J=1,NSCAN
            IFSCAN(J)=.FALSE.
          END DO
          DO I=1,NSKYREG
            DO J=L1(I),L2(I)
              IFSCAN(J)=.TRUE.
            END DO
          END DO
          GOTO 117
        END IF
C
        IF(NI.EQ.-4)THEN
          NSEED=-1                                      !inicializaremos RANRED
          WRITE(*,100) 'No. of simulations to compute errors '
          WRITE(CDUMMY,*) NSIMULT
          NSIMULT=READILIM(CDUMMY,10,NMAXSIMUL)
          GOTO 117
        END IF
C
30      CONTINUE
        M=0
        DO J=1,NSCAN
          IF(IFSCAN(J))THEN
            M=M+1
            X(M)=REAL(J)
            Y(M)=G(NI,J)
            END IF
        END DO
        NPTS=M
C
        DO I=1,NSCAN
          XT(I)=REAL(I)
          YT(I)=G(NI,I)
        END DO

ccc        CALL POLFIT(X,Y,Y,NPTS,NTERMS,0,A,CHISQR)
        CALL POLFITSIG(NPTS,X,Y,TSIGMA,NTERMS-1,A,LFIT2)

        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,-2)
          CALL PGBOX('BCTNSI',0.0,0,'BCTNSI',0.0,0)
          CALL PGIDEN_RED
          WRITE(GLABEL,'(A,I4)') 'Channel #',NI
          CALL PGLABEL('Scan','No. counts',GLABEL(1:TRUELEN(GLABEL))
     +     //'   File: '//INFILE)
          CALL PGBIN(NSCAN,XT,YT,.TRUE.)
        END DO
C
        DO I=1,NSCAN
          XP=XT(I)
          POL=A(NTERMS)
          DO J=NTERMS-1,1,-1
             POL=POL*XP+A(J)
          END DO
          SK(I)=POL
        END DO
C
        IF(CERRG.EQ.'y')THEN !si hemos decidido estimar errores con el r.m.s.
          ERRSKYG=0. !almacenamos aqui el sumatorio
          IF(NPTS.GT.1)THEN
            DO J=1,NSCAN
              IF(IFSCAN(J))THEN
                ERRSKYG=ERRSKYG+(G(NI,J)-SK(J))*(G(NI,J)-SK(J))
              END IF
            END DO
            ERRSKYG=SQRT(ERRSKYG/REAL(NPTS-1))
            DO J=1,NSCAN
              SKYSIGMA1G(J)=SK(J)+ERRSKYG
              SKYSIGMA2G(J)=SK(J)-ERRSKYG
            END DO
          END IF
        END IF
C
        IF(CERR.EQ.'y')THEN
          DO NSIMUL=1,NSIMULT
            M=0
            ERRSKY=0.
            DO J=1,NSCAN
              IF(IFSCAN(J))THEN
                M=M+1
                X(M)=REAL(J)
                Y(M)=G(NI,J)
                IF(ERR(NI,J).GT.0.0)THEN
                  ERRSKY=ERRSKY+1./(ERR(NI,J)*ERR(NI,J))
                END IF
                R1=RANRED(NSEED)
                R2=RANRED(NSEED)
                GXX=1.414213562*ERR(NI,J)*SQRT(-1.*LOG(1.-R1))*
     >           COS(PI2*R2)
                Y(M)=Y(M)+GXX
                END IF
            END DO
            NPTS=M
            IF(ERRSKY.GT.0.0)THEN
              ERRSKY=1./SQRT(ERRSKY)
            ELSE
              ERRSKY=0.0
            END IF
ccc            CALL POLFIT(X,Y,Y,NPTS,NTERMS,0,A,CHISQR)
            CALL POLFITSIG(NPTS,X,Y,TSIGMA,NTERMS-1,A,LFIT2)
            DO I=1,NSCAN
              XP=XT(I)
              POL=A(NTERMS)
              DO J=NTERMS-1,1,-1
                 POL=POL*XP+A(J)
              END DO
              SKYRAN(I,NSIMUL)=POL
            END DO
          END DO
          DO I=1,NSCAN
            MEAN=0.
            DO NSIMUL=1,NSIMULT
              MEAN=MEAN+SKYRAN(I,NSIMUL)
            END DO
            MEAN=MEAN/REAL(NSIMULT)
            SIGMA=0.
            DO NSIMUL=1,NSIMULT
              SIGMA=SIGMA+(MEAN-SKYRAN(I,NSIMUL))*
     >         (MEAN-SKYRAN(I,NSIMUL))
            END DO
            SIGMA=SQRT(SIGMA/REAL(NSIMULT-1))
            SKSIGMA1(I)=SK(I)+SIGMA
            SKSIGMA2(I)=SK(I)-SIGMA
            SKYSIGMA1(I)=SK(I)+ERRSKY
            SKYSIGMA2(I)=SK(I)-ERRSKY
          END DO
        END IF
C
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM))CALL PGSCI(2)
          CALL PGLINE(NSCAN,XT,SK)
          IF(CERR.EQ.'y')THEN
            CALL PGSLS(2)
            IF(LCOLOR(ITERM))CALL PGSCI(5)
            CALL PGLINE(NSCAN,XT,SKSIGMA1)
            CALL PGLINE(NSCAN,XT,SKSIGMA2)
            CALL PGSLS(4)
            IF(LCOLOR(ITERM))CALL PGSCI(4)
            CALL PGLINE(NSCAN,XT,SKYSIGMA1)
            CALL PGLINE(NSCAN,XT,SKYSIGMA2)
          END IF
          IF(CERRG.EQ.'y')THEN
            IF(LCOLOR(ITERM))CALL PGSCI(7)
            CALL PGLINE(NSCAN,XT,SKYSIGMA1G)
            CALL PGLINE(NSCAN,XT,SKYSIGMA2G)
          END IF
          CALL PGSLS(1)
          IF(LCOLOR(ITERM))CALL PGSCI(3)
          DO I=1,NSKYREG
            CALL PGMOVE(REAL(L1(I)),YMIN)
            CALL PGDRAW(REAL(L1(I)),YMAX)
            CALL PGMOVE(REAL(L2(I)),YMIN)
            CALL PGDRAW(REAL(L2(I)),YMAX)
            CALL PGRECT(REAL(L1(I)),REAL(L2(I)),
     +       YMIN,YMIN+(YMAX-YMIN)/80.)
            CALL PGRECT(REAL(L1(I)),REAL(L2(I)),
     +       YMAX,YMAX-(YMAX-YMIN)/80.)
          END DO
          IF(LCOLOR(ITERM))CALL PGSCI(1)
        END DO
        IF(CERR.EQ.'y')THEN
          WRITE(*,101) '-> Dashed line (light blue): realistic errors'
          WRITE(*,101) '-> Dotted line (dark blue): approximated errors'
        END IF
        IF(CERRG.EQ.'y')THEN
          WRITE(*,101) '-> Dotted line (yellow): error from r.m.s.'
        END IF
        GOTO 117
120     CONTINUE
        CALL PGEND
C------------------------------------------------------------------------------
130     WRITE(*,*)
        WRITE(*,101) 'Final Sky Region Settings:'
        DO I=1,NSKYREG
          WRITE(*,'(A,I2.2,5X,I4,2X,I4)') 'Region #',I,L1(I),L2(I)
        END DO
        IF(CPLOT.NE.'n')WRITE(*,'(A,I2)') 'Polynomial degree: ',NTERMS-1
        WRITE(*,*)
C
        IF(CERR.EQ.'y')THEN
          WRITE(*,101) '(1) Errors using approximated formulae'
          WRITE(*,101) '(2) More realistic errors using simulations'
          WRITE(*,101) 'NOTE: (1)=(2) for pol. degree=0'
          WRITE(*,100) 'Option (1/2) '
          CERR2(1:1)=READC('@','12')
          IF(CERR2.EQ.'2')THEN
            NSEED=-1                                    !inicializaremos RANRED
            WRITE(*,100) 'No. of simulations to compute errors '
            WRITE(CDUMMY,*) NSIMULT
            NSIMULT=READILIM(CDUMMY,10,NMAXSIMUL)
          END IF
        END IF
        WRITE(*,*)
C
        WRITE(*,100) 'Polynomial degree '
        WRITE(CDUMMY,*) NTERMS-1
        NTERMS=READILIM(CDUMMY,0,19)
        NTERMS=NTERMS+1
        WRITE(CDUMMY,*) TSIGMA
        WRITE(*,100) 'Times sigma to remove points '
        TSIGMA=READF(CDUMMY)
        WRITE(*,*)
        WRITE(*,100) '* Fitting sky...'
C------------------------------------------------------------------------------
        DO NI=1,NCHAN
          M=0
          DO J=1,NSCAN
            IF(IFSCAN(J))THEN
              M=M+1
              X(M)=REAL(J)
              Y(M)=G(NI,J)
              END IF
          END DO
          NPTS=M
ccc          CALL POLFIT(X,Y,Y,NPTS,NTERMS,0,A,CHISQR)
          CALL POLFITSIG(NPTS,X,Y,TSIGMA,NTERMS-1,A,LFIT2)
          DO I=1,NSCAN
            XP=REAL(I)
            POL=A(NTERMS)
            DO J=NTERMS-1,1,-1
               POL=POL*XP+A(J)
            END DO
            SKY(NI,I)=POL
          END DO
          IF(CERRG.EQ.'y')THEN !si hemos decidido estimar errores con el r.m.s.
            ESKYG(NI,1)=0. !almacenamos aqui el sumatorio
            IF(NPTS.GT.1)THEN
              DO J=1,NSCAN
                IF(IFSCAN(J))THEN
                  ESKYG(NI,1)=ESKYG(NI,1)+
     +             (G(NI,J)-SKY(NI,J))*(G(NI,J)-SKY(NI,J))
                END IF
              END DO
              ESKYG(NI,1)=SQRT(ESKYG(NI,1)/REAL(NPTS-1))
            END IF
            DO J=2,NSCAN !metemos en toda la columna el mismo valor del error
              ESKYG(NI,J)=ESKYG(NI,1)
            END DO
          END IF
        END DO
        WRITE(*,101) '  ...OK'
C------------------------------------------------------------------------------
        IF(CERR.EQ.'y')THEN
          IF(CERR2.EQ.'2')THEN             !errores calculados con simulaciones
            WRITE(*,101) '* Computing errors:'
            DO NI=1,NCHAN
              WRITE(*,'(A,I4.4,$)') '\b\b\b\b',NI
              DO NSIMUL=1,NSIMULT
                M=0
                DO J=1,NSCAN
                  IF(IFSCAN(J))THEN
                    M=M+1
                    X(M)=REAL(J)
                    Y(M)=G(NI,J)
                    R1=RANRED(NSEED)
                    R2=RANRED(NSEED)
                    GXX=1.414213562*ERR(NI,J)*SQRT(-1.*LOG(1.-R1))*
     >               COS(PI2*R2)
                    Y(M)=Y(M)+GXX
                    END IF
                END DO
                NPTS=M
ccc                CALL POLFIT(X,Y,Y,NPTS,NTERMS,0,A,CHISQR)
                CALL POLFITSIG(NPTS,X,Y,TSIGMA,NTERMS-1,A,LFIT2)
                DO I=1,NSCAN
                  XP=REAL(I)
                  POL=A(NTERMS)
                  DO J=NTERMS-1,1,-1
                     POL=POL*XP+A(J)
                  END DO
                  SKYRAN(I,NSIMUL)=POL
                END DO
              END DO
              DO I=1,NSCAN
                MEAN=0.
                DO NSIMUL=1,NSIMULT
                  MEAN=MEAN+SKYRAN(I,NSIMUL)
                END DO
                MEAN=MEAN/REAL(NSIMULT)
                SIGMA=0.
                DO NSIMUL=1,NSIMULT
                  SIGMA=SIGMA+(MEAN-SKYRAN(I,NSIMUL))*
     >             (MEAN-SKYRAN(I,NSIMUL))
                END DO
                SIGMA=SIGMA/REAL(NSIMULT-1)                   !OJO: es sigma**2
                ERR(NI,I)=ERR(NI,I)*ERR(NI,I)+SIGMA          !error al cuadrado
                ESKY(NI,I)=SQRT(SIGMA)                  !error solo en el cielo
              END DO
            END DO
            WRITE(*,*)
          ELSE                          !errores calculados de forma aproximada
            WRITE(*,100) '* Computing errors...'
            DO NI=1,NCHAN
              ERRSKY=0.
              DO J=1,NSCAN
                IF(IFSCAN(J))THEN
                  IF(ERR(NI,J).GT.0.0)THEN
                    ERRSKY=ERRSKY+1./(ERR(NI,J)*ERR(NI,J))
                  END IF
                  END IF
              END DO
              IF(ERRSKY.GT.0.0)THEN
                ERRSKY=1./SQRT(ERRSKY)
              ELSE
                ERRSKY=0.0
              END IF
              DO J=1,NSCAN
                ERR(NI,J)=ERR(NI,J)*ERR(NI,J)+ERRSKY*ERRSKY  !error al cuadrado
                ESKY(NI,I)=ERRSKY                       !error solo en el cielo
              END DO
            END DO
            WRITE(*,101) '  ...OK'
          END IF
        END IF
C------------------------------------------------------------------------------
        WRITE(*,100) 'Save sky fit into file (y/n) '
        CSAVSK(1:1)=READC('n','yn')
        IF(CSAVSK.EQ.'y')THEN
          WRITE(*,100) 'Output sky fit file name'
          SKYFILE=OUTFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
          DO I=1,NSCAN
            WRITE(20) (SKY(J,I),J=1,NCHAN)
          END DO
          CLOSE(20)
          IF(CERR.EQ.'y')THEN
            WRITE(*,100) 'Output error file name '
            CALL GUESSEF(SKYFILE,ERRFILE)
            OUTFILE=OUTFILEX(21,ERRFILE,NSCAN,NCHAN,STWV,DISP,1,.TRUE.)
            DO I=1,NSCAN
              WRITE(21) (ESKY(J,I),J=1,NCHAN)
            END DO
            CLOSE(21)
          END IF
        END IF
C
        DO I=1,NSCAN
          DO J=1,NCHAN
            G(J,I)=G(J,I)-SKY(J,I)
          END DO
        END DO
        WRITE(*,100) 'Output file name'
        OUTFILE=OUTFILEX(15,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        DO I=1,NSCAN
          WRITE(15) (G(J,I),J=1,NCHAN)
        END DO
        CLOSE(15)
        IF(CERR.EQ.'y')THEN
          DO I=1,NSCAN
            DO J=1,NCHAN
              ERR(J,I)=SQRT(ERR(J,I))
            END DO
          END DO
          WRITE(*,100) 'Output error file name '
          CALL GUESSEF(OUTFILE,ERRFILE)
          OUTFILE=OUTFILEX(16,ERRFILE,NSCAN,NCHAN,STWV,DISP,1,.TRUE.)
          DO I=1,NSCAN
            WRITE(16) (ERR(J,I),J=1,NCHAN)
          END DO
          CLOSE(16)
        END IF
        IF(CERRG.EQ.'y')THEN
          WRITE(*,100) 'Output file name for error estimates '
          WRITE(*,100) 'computed from r.m.s. in fits '
          OUTFILE=OUTFILEX(17,'@',NSCAN,NCHAN,STWV,DISP,1,.TRUE.)
          DO I=1,NSCAN
            WRITE(17) (ESKYG(J,I),J=1,NCHAN)
          END DO
          CLOSE(17)
        END IF
C
        STOP
100     FORMAT(A,$)
101     FORMAT(A)
110     FORMAT(A,I6)
C
        END
C
C******************************************************************************
C
        SUBROUTINE CHANGESKY
        IMPLICIT NONE
        INCLUDE 'futils.inc'
        INTEGER READILIM
C
        INTEGER I,NSKY
        INTEGER NSKYREG,L1(20),L2(20)
        CHARACTER*1 COPC
        COMMON/BLKSKY/NSKYREG,L1,L2
C
10      CONTINUE
        WRITE(*,*)
        WRITE(*,101) '----------------------------------'
        DO I=1,NSKYREG
          WRITE(*,'(A,I2.2,5X,I4,2X,I4)') 'Region #',I,L1(I),L2(I)
        END DO
        WRITE(*,101) '----------------------------------'
        WRITE(*,101) '(1) Remove region'
        WRITE(*,101) '(2) Add region'
        WRITE(*,101) '(3) Change region limits'
        WRITE(*,101) '(0) EXIT'
        WRITE(*,100) 'Option '
        COPC(1:1)=READC('0','0123')
        IF(COPC.EQ.'0')RETURN
C
        IF(COPC.EQ.'1')THEN
          WRITE(*,100) 'Region number to be removed (0=none) '
          NSKY=READILIM('@',1,NSKYREG)
          DO I=NSKY,NSKYREG-1
            L1(I)=L1(I+1)
            L2(I)=L2(I+1)
          END DO
          NSKYREG=NSKYREG-1
        END IF
C
        IF(COPC.EQ.'2')THEN
          NSKYREG=NSKYREG+1
          IF(NSKYREG.GT.20)THEN
            NSKYREG=20
            WRITE(*,100) 'Number of sky region = maximum number = 20'
            GOTO 10
          END IF
          WRITE(*,100) 'left limit, right limit'
          CALL READ2I('@',L1(NSKYREG),L2(NSKYREG))
        END IF
C
        IF(COPC.EQ.'3')THEN
          WRITE(*,100) 'Region number to be changed '
          NSKY=READILIM('@',1,NSKYREG)
          WRITE(*,100) 'New limits: left limit, right limit'
          CALL READ2I('@',L1(NSKY),L2(NSKY))
        END IF
C
        GOTO 10
C
100     FORMAT(A,$)
101     FORMAT(A)
        END
