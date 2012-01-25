C------------------------------------------------------------------------------
C Version 26-November-2007                                     file: fit1dpol.f
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
C Program: fit1dpol
C Classification: arithmetic & manipulations
C Description: Fits 1-D polynomials.
C
Comment
C
C Ajusta un polinomio a cortes individuales/sumados de una imagen
C
        PROGRAM FIT1DPOL
C
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER TRUELEN
        INTEGER READILIM
        REAL READF
        REAL FPOLY
C
        INTEGER MAXNTERMS
        PARAMETER(MAXNTERMS=20)
C
        INTEGER I,J,L
        INTEGER K,KREMOVED
        INTEGER NS,NS1,NS2
        INTEGER NC,NC1,NC2
        INTEGER NF1,NF2
        INTEGER NCOEFF,NPTOS
        INTEGER DEG
        INTEGER NPT
        INTEGER NSMAX_LOCAL,NCMAX_LOCAL
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        INTEGER IMODE,I0
        REAL S(NCMAX,NSMAX)
        REAL X(NCMAX),Y(NCMAX)                            !mayor de NCMAX/NSMAX
        REAL XP(NCMAX),YP(NCMAX)                          !mayor de NCMAX/NSMAX
        REAL XF(NCMAX),YF(NCMAX)
        REAL XREM(NCMAX),YREM(NCMAX)
        REAL XXX,POL
        REAL XMAX,XMIN,YMAX,YMIN
        REAL CHISQR,A(MAXNTERMS)
        REAL TSIGMA
        REAL S_SPL(NCMAX),A_SPL(NCMAX),B_SPL(NCMAX),C_SPL(NCMAX)
        REAL YRMSTOL,PSEUDO_WEIGHT,PSEUDO_POWER
        DOUBLE PRECISION DSIGMA
        CHARACTER*1 COPC,CSAVE,COTH,CLUP
        CHARACTER*5 COPCALL
        CHARACTER*50 CDUMMY
        CHARACTER*75 INFILE,OUTFILE
        LOGICAL IFSCAN(NSMAX)
        LOGICAL IFCHAN(NCMAX)
        LOGICAL IFFIT(NCMAX)                             !mayor de NCMAX/NSMAX
        LOGICAL LFIRST
        LOGICAL IFREMOVED1(NCMAX),IFREMOVED2(NCMAX),LREPEAT
        LOGICAL LCOLOR(MAX_ID_RED)
        LOGICAL PSEUDO_LUP
C
        COMMON/BLKESC/XMIN,XMAX,YMIN,YMAX
C------------------------------------------------------------------------------
        THISPROGRAM='fit1dpol'
        CALL WELCOME('26-November-2007')
C
        NSMAX_LOCAL=NSMAX
        NCMAX_LOCAL=NCMAX
        IF(NSMAX_LOCAL.GT.NCMAX_LOCAL)STOP 'FATAL ERROR: NSMAX.GT.NCMAX'
C
        LFIRST=.TRUE.
        COPCALL='0'
        YRMSTOL=0.000001
        PSEUDO_WEIGHT=100.0
        PSEUDO_POWER=2.0
        PSEUDO_LUP=.TRUE.
C
        CALL PIDEGTER(NTERM,IDN,LCOLOR)
C
C leemos imagen
        WRITE(*,100)'Input file name'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        IF(TRUELEN(OBJECT).GT.0)THEN
          INFILE=INFILE(1:TRUELEN(INFILE))//' ['//
     +     OBJECT(1:TRUELEN(OBJECT))//']'
        END IF
        DO I=1,NSCAN
          READ(20) (S(J,I),J=1,NCHAN)
        END DO
        CLOSE(20)
C------------------------------------------------------------------------------
C establecemos opciones posibles
5       WRITE(*,*)
        IF(NCHAN.GT.1)THEN
           WRITE(*,101)'(1) Fit a single spectrum'
           L=TRUELEN(COPCALL)
           COPCALL=COPCALL(1:L)//'1'
           IF(NSCAN.GT.1)THEN
             WRITE(*,101)'(2) Fit averaged added spectra'
             L=TRUELEN(COPCALL)
             COPCALL=COPCALL(1:L)//'2'
           END IF
        END IF
        IF(NSCAN.GT.1)THEN
          WRITE(*,101)'(3) Fit a single spatial direction'
          L=TRUELEN(COPCALL)
          COPCALL=COPCALL(1:L)//'3'
          IF(NCHAN.GT.1)THEN
            WRITE(*,101)'(4) Fit averaged added spatial directions'
            L=TRUELEN(COPCALL)
            COPCALL=COPCALL(1:L)//'4'
          END IF
        END IF
C
        WRITE(*,101)'(0) EXIT'
        WRITE(*,*)
ccc10      WRITE(*,100)'Option'
        WRITE(*,100)'Option'
        L=TRUELEN(COPCALL)
        IF(LFIRST)THEN
          LFIRST=.FALSE.
          IF(L.EQ.2)THEN
            WRITE(*,100)' '
            COPC(1:1)=READC(COPCALL(2:2),COPCALL(1:L))
          ELSE
            COPC(1:1)=READC('@',COPCALL(1:L))
          END IF
        ELSE
          WRITE(*,100)' '
          COPC(1:1)=READC('0',COPCALL(1:L))
        END IF
        IF(COPC.EQ.'0')THEN
          CALL PGEND
          STOP
        END IF
C------------------------------------------------------------------------------
ccc15      IF(COPC.EQ.'1')THEN
        IF(COPC.EQ.'1')THEN
          NPT=NCHAN
          IF(NSCAN.GT.1)THEN
            WRITE(*,100)'Valid region: from 1 to '
            WRITE(CDUMMY,*)NSCAN
            CALL RMBLANK(CDUMMY,CDUMMY,L)
            WRITE(*,101)CDUMMY(1:L)
            WRITE(*,100)'Spectrum number '
            NS=READILIM('@',1,NSCAN)
          ELSE
            NS=1
          END IF
          DO J=1,NPT
            X(J)=REAL(J)
            Y(J)=S(J,NS)
          END DO
C------------------------------------------------------------------------------
        ELSEIF(COPC.EQ.'2')THEN
          NPT=NCHAN
          DO I=1,NSCAN
            IFSCAN(I)=.FALSE.
          END DO
          WRITE(*,100)'Valid region: from 1 to '
          WRITE(CDUMMY,*)NSCAN
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          WRITE(*,101)CDUMMY(1:L)
22        WRITE(*,100)'First and last scan (0,0=EXIT) '
          CALL READ2I('0,0',NS1,NS2)
          IF((NS1.EQ.0).AND.(NS2.EQ.0)) GOTO 24
          IF((NS1.LT.1).OR.(NS2.GT.NSCAN).OR.(NS1.GT.NS2))THEN
            WRITE(*,101)'ERROR: invalid entry. Try again.'
            GOTO 22
          END IF
          DO I=NS1,NS2
            IFSCAN(I)=.TRUE.
          END DO
          GOTO 22
24        DO J=1,NPT
            Y(J)=0.
          END DO
          K=0
          DO I=1,NSCAN
            IF(IFSCAN(I))THEN
              K=K+1
              DO J=1,NPT
                Y(J)=Y(J)+S(J,I)
              END DO
            END IF
          END DO
          IF(K.EQ.0)THEN
            WRITE(*,101)'ERROR: no scans added.'
            GOTO 5
          END IF
          DO J=1,NPT
            X(J)=REAL(J)
            Y(J)=Y(J)/REAL(K)
          END DO
C------------------------------------------------------------------------------
        ELSEIF(COPC.EQ.'3')THEN
          NPT=NSCAN
          IF(NCHAN.GT.1)THEN
            WRITE(*,100)'Valid region: from 1 to '
            WRITE(CDUMMY,*)NCHAN
            CALL RMBLANK(CDUMMY,CDUMMY,L)
            WRITE(*,101)CDUMMY(1:L)
            WRITE(*,100)'Channel number '
            NC=READILIM('@',1,NCHAN)
          ELSE
            NC=1
          END IF
          DO I=1,NPT
            X(I)=REAL(I)
            Y(I)=S(NC,I)
          END DO
C------------------------------------------------------------------------------
        ELSEIF(COPC.EQ.'4')THEN
          NPT=NSCAN
          DO J=1,NCHAN
            IFCHAN(J)=.FALSE.
          END DO
          WRITE(*,100)'Valid region: from 1 to '
          WRITE(CDUMMY,*)NCHAN
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          WRITE(*,101)CDUMMY(1:L)
28        WRITE(*,100)'First and last channel (0,0=EXIT) '
          CALL READ2I('0,0',NC1,NC2)
          IF((NC1.EQ.0).AND.(NC2.EQ.0)) GOTO 30
          IF((NC1.LT.1).OR.(NC2.GT.NCHAN).OR.(NC1.GT.NC2))THEN
            WRITE(*,101)'ERROR: invalid entry. Try again.'
            GOTO 28
          END IF
          DO I=NC1,NC2
            IFCHAN(I)=.TRUE.
          END DO
          GOTO 28
30        DO J=1,NPT
            Y(J)=0.
          END DO
          K=0
          DO J=1,NCHAN
            IF(IFCHAN(J))THEN
              K=K+1
              DO I=1,NPT
                Y(I)=Y(I)+S(J,I)
              END DO
            END IF
          END DO
          IF(K.EQ.0)THEN
            WRITE(*,101)'ERROR: no channels added.'
            GOTO 5
          END IF
          DO J=1,NPT
            X(J)=REAL(J)
            Y(J)=Y(J)/REAL(K)
          END DO
        END IF
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM))THEN
            CALL PLOTSP(NPT,X,Y,.TRUE.,7,1,ITERM)
            CALL PGSCI(1)
          ELSE
            CALL PLOTSP(NPT,X,Y,.TRUE.,1,1,ITERM)
          END IF
          IF((COPC.EQ.'1').OR.(COPC.EQ.'2'))THEN
            CALL PGLABEL('channel','No. counts','file: '//INFILE)
          ELSE
            CALL PGLABEL('scan','No. counts','file: '//INFILE)
          END IF
        END DO
C------------------------------------------------------------------------------
C determinamos puntos a ajustar
        DO J=1,NPT
          IFFIT(J)=.FALSE.
        END DO
        WRITE(*,*)
        WRITE(*,100)'Valid region: from 1 to '
        WRITE(CDUMMY,*)NPT
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,101)CDUMMY(1:L)
        WRITE(*,101)'Enter regions to be fitted:'
50      WRITE(*,100)'First and last point (0,0=EXIT) '
        CALL READ2I('0,0',NF1,NF2)
        IF((NF1.EQ.0).AND.(NF2.EQ.0)) GOTO 60
        IF((NF1.LT.1).OR.(NF2.GT.NPT).OR.(NF1.GT.NF2))THEN
          WRITE(*,101)'ERROR: invalid entry. Try again.'
          GOTO 50
        END IF
        DO J=NF1,NF2
          IFFIT(J)=.TRUE.
        END DO
        GOTO 50
60      K=0
        DO J=1,NPT
          IF(IFFIT(J))THEN
            K=K+1
            XF(K)=X(J)
            YF(K)=Y(J)
          END IF
        END DO
        NPTOS=K
        WRITE(*,*)
        IF(NPTOS.EQ.0)THEN
          WRITE(*,101)'ERROR: No. of points = 0'
          GOTO 5
        ELSE
          WRITE(*,110)'>>> Total no. of points for initial '//
     +     'fit: ',NPTOS
        END IF
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM)) CALL PGSCI(3)
          CALL PGPOINT(NPTOS,XF,YF,17)
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
        END DO
C------------------------------------------------------------------------------
70      WRITE(*,100)'Polynomial degree (-1=splines,-2=pseudocontinuum)'
        DEG=READILIM('@',-2,MAXNTERMS-1)
        IF(DEG.EQ.-2)THEN
          WRITE(CDUMMY,*) YRMSTOL
          WRITE(*,100)'YRMSTOL for DOWNHILL '
          YRMSTOL=READF(CDUMMY)
          WRITE(CDUMMY,*) PSEUDO_WEIGHT
          WRITE(*,100)'WEIGHT for PSEUDOFIT '
          PSEUDO_WEIGHT=READF(CDUMMY)
          WRITE(CDUMMY,*) PSEUDO_POWER
          WRITE(*,100)'POWER for PSEUDOFIT '
          PSEUDO_POWER=READF(CDUMMY)
          IF(PSEUDO_LUP)THEN
            CLUP='y'
          ELSE
            CLUP='n'
          END IF
          WRITE(*,100)'Fit upper side (y/n) '
          CLUP(1:1)=READC(CLUP,'yn')
          PSEUDO_LUP=(CLUP.EQ.'y')
          WRITE(*,100)'Polynomial degree'
          NCOEFF=READILIM('@',0,MAXNTERMS-1)
          NCOEFF=NCOEFF+1
          CALL PSEUDOFIT(XF,YF,NPTOS,NCOEFF,YRMSTOL,
     +     PSEUDO_WEIGHT,PSEUDO_POWER,PSEUDO_LUP,A)
          DO J=1,NPT
            XP(J)=REAL(J)
            YP(J)=FPOLY(NCOEFF-1,A,XP(J))
          END DO
          WRITE(*,101)'Coefficients from fit: (y=a0+a1*x+a2*x*x+...)'
          DO I=1,NCOEFF
            IF(I.LT.11)THEN
              WRITE(*,'(A,I1,A,$)')'> a(',I-1,') : '
            ELSE
              WRITE(*,'(A,I2,A,$)')'> a(',I-1,'): '
            END IF
            WRITE(*,*)A(I)
          END DO
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM))THEN
              CALL PGSCI(5)
              CALL PGLINE(NPT,XP,YP)
              CALL PGSCI(1)
            ELSE
              CALL PGLINE(NPT,XP,YP)
            END IF
          END DO
          GOTO 80
        ELSEIF(DEG.EQ.-1)THEN
          WRITE(*,100) 'IMODE'
          IMODE=READILIM('@',1,4)
          CALL CUBSPL(XF,YF,NPTOS,IMODE,S_SPL,A_SPL,B_SPL,C_SPL)
          DO J=1,NPT
            XP(J)=REAL(J)
            CALL CUBSPLX(XF,YF,A_SPL,B_SPL,C_SPL,NPTOS,I0,XP(J),YP(J))
          END DO
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM))THEN
              CALL PGSCI(5)
              CALL PGLINE(NPT,XP,YP)
              CALL PGSCI(1)
            ELSE
              CALL PGLINE(NPT,XP,YP)
            END IF
          END DO
          GOTO 80
        END IF
        NCOEFF=DEG+1
C
71      WRITE(*,100)'Times sigma to remove points '//
     +   '(0.0 = do not remove) '
        TSIGMA=READF('3.0')
        IF(TSIGMA.LT.0.0)THEN
          WRITE(*,101)'ERROR: this value must be > 0!'
          GOTO 71
        ELSEIF(TSIGMA.EQ.0.0)THEN
        ELSE
          DO I=1,NPT
            IFREMOVED1(I)=.FALSE.
          END DO
        END IF
C------------------------------------------------------------------------------
C ajustamos polinomio
73      CALL POLFIT(XF,YF,YF,NPTOS,NCOEFF,0,A,CHISQR)
C si no iteramos, hemos terminado
        KREMOVED=0
        IF(TSIGMA.EQ.0.0) GOTO 75
C calculamos sigma solo con los puntos que hemos ajustado
        DSIGMA=0.D0
        IF(NPTOS.GT.1)THEN
          DO J=1,NPTOS
            XXX=XF(J)
            POL=A(NCOEFF)
            DO K=NCOEFF-1,1,-1
              POL=POL*XXX+A(K)
            END DO
            DSIGMA=DSIGMA+DBLE(YF(J)-POL)*DBLE(YF(J)-POL)
          END DO
          DSIGMA=SQRT(DSIGMA/DBLE(NPTOS-1))
        END IF
        WRITE(*,100)'>>> Sigma around the fit: '
        WRITE(*,*) REAL(DSIGMA)
C eliminamos los puntos que se desvian
        DO J=1,NPT
          IF(IFFIT(J))THEN
            XXX=X(J)
            POL=A(NCOEFF)
            DO K=NCOEFF-1,1,-1
              POL=POL*XXX+A(K)
            END DO
            IF(ABS(Y(J)-POL).GT.REAL(DSIGMA)*TSIGMA)THEN
              KREMOVED=KREMOVED+1
              XREM(KREMOVED)=X(J)
              YREM(KREMOVED)=Y(J)
              IFREMOVED2(J)=.TRUE.
            ELSE
              IFREMOVED2(J)=.FALSE.
            END IF
          END IF
        END DO
        WRITE(*,110)'>>> No. of points removed: ',KREMOVED
C comprobamos si hemos cambiado algo
        LREPEAT=.FALSE.
        DO I=1,NPT
          IF(IFREMOVED1(I).NEQV.IFREMOVED2(I)) LREPEAT=.TRUE.
        END DO
        IF(.NOT.LREPEAT) GOTO 75
        DO I=1,NPT
          IFREMOVED1(I)=IFREMOVED2(I)
        END DO
        K=0
        DO J=1,NPT
          IF(IFFIT(J).AND.(.NOT.IFREMOVED1(J)))THEN
            K=K+1
            XF(K)=X(J)
            YF(K)=Y(J)
          END IF
        END DO
        NPTOS=K
        WRITE(*,*)
        IF(NPTOS.EQ.0)THEN
          WRITE(*,101)'ERROR: No. of points = 0'
          GOTO 5
        ELSE
          WRITE(*,110)'>>> Total no. of points for next '//
     +     'fit: ',NPTOS
        END IF
        GOTO 73
C------------------------------------------------------------------------------
75      DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM))THEN
            CALL PLOTSP(NPT,X,Y,.FALSE.,5,1,ITERM)
            CALL PGSCI(1)
          ELSE
            CALL PLOTSP(NPT,X,Y,.FALSE.,1,1,ITERM)
          END IF
          IF((COPC.EQ.'1').OR.(COPC.EQ.'2'))THEN
            CALL PGLABEL('channel','No. counts','file: '//INFILE)
          ELSE
            CALL PGLABEL('scan','No. counts','file: '//INFILE)
          END IF
          IF(LCOLOR(ITERM)) CALL PGSCI(3)
          CALL PGPOINT(NPTOS,XF,YF,17)
          IF(KREMOVED.GT.0)THEN         !mostramos puntos eliminados del ajuste
            IF(LCOLOR(ITERM)) CALL PGSCI(2)
            CALL PGPOINT(KREMOVED,XREM,YREM,5)
          END IF
        END DO
C mostramos el ajuste
        CALL POLFIT(XF,YF,YF,NPTOS,NCOEFF,0,A,CHISQR)
        DO J=1,NPT
          XP(J)=REAL(J)
          YP(J)=A(NCOEFF)
          DO K=NCOEFF-1,1,-1
            YP(J)=YP(J)*XP(J)+A(K)
          END DO
        END DO
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM))THEN
            CALL PGSCI(5)
            CALL PGLINE(NPT,XP,YP)
            CALL PGSCI(1)
          ELSE
            CALL PGLINE(NPT,XP,YP)
          END IF
        END DO
C mostramos los coeficientes
        WRITE(*,*)
        WRITE(*,101)'>>> Fit result: y= a0 + a1 x + a2 * x^2 + ...'
        DO I=1,NCOEFF
          IF(I.LT.11)THEN
            WRITE(*,'(A,I1.1,A,$)')'    a[',I-1,']= '
          ELSE
            WRITE(*,'(A,I2.2,A,$)')'   a[',I-1,']= '
          END IF
          WRITE(*,*)A(I)
        END DO
        WRITE(*,*)
C------------------------------------------------------------------------------
80      WRITE(*,100)'Other fit (y/n) '
        COTH(1:1)=READC('n','yn')
        IF(COTH.EQ.'y') GOTO 70
C------------------------------------------------------------------------------
ccc90      WRITE(*,100)'Save fit into output file (y/n) '
        WRITE(*,100)'Save fit into output file (y/n) '
        CSAVE(1:1)=READC('y','yn')
        IF(CSAVE.EQ.'y')THEN
          WRITE(*,100)'Output file name'
          IF((COPC.EQ.'1').OR.(COPC.EQ.'2'))THEN
            OUTFILE=OUTFILEX(30,'@',1,NCHAN,STWV,DISP,1,.FALSE.)
            WRITE(30) (YP(J),J=1,NCHAN)
          ELSE
            OUTFILE=OUTFILEX(30,'@',NSCAN,1,0.,0.,1,.FALSE.)
            DO I=1,NSCAN
              WRITE(30) YP(I)
            END DO
          END IF
          CLOSE(30)
        END IF
        GOTO 5
C------------------------------------------------------------------------------
C
100     FORMAT(A,$)
101     FORMAT(A)
110     FORMAT(A,I6)
        END
C
C*************************************************************************
C
        SUBROUTINE PLOTSP(N,X,Y,LLIMIT,NCOLOR,NSYMB,ITERM)
        IMPLICIT NONE
        INCLUDE 'futils.inc'
        REAL READF
C
        INTEGER N,NCOLOR,NSYMB
        REAL X(N),Y(N)
        LOGICAL LLIMIT
        INTEGER ITERM
C
        INTEGER I
        REAL XMIN,XMAX,YMIN,YMAX
        REAL DX,DY
        CHARACTER*1 CCHANGE
        CHARACTER*50 CDUMMY
C
        COMMON/BLKESC/XMIN,XMAX,YMIN,YMAX
C-------------------------------------------------------------------------
        IF(N.LE.1)THEN
          WRITE(*,101)'ERROR: number of points.LE.1 (PLOTSP)'
          RETURN
        END IF
        IF(.NOT.LLIMIT) GOTO 5
        XMIN=X(1)
        XMAX=X(1)
        YMIN=Y(1)
        YMAX=Y(1)
        DO I=2,N
          IF(XMIN.GT.X(I)) XMIN=X(I)
          IF(XMAX.LT.X(I)) XMAX=X(I)
          IF(YMIN.GT.Y(I)) YMIN=Y(I)
          IF(YMAX.LT.Y(I)) YMAX=Y(I)
        END DO
        DX=XMAX-XMIN
        DY=YMAX-YMIN
        XMIN=XMIN-DX*.05
        XMAX=XMAX+DX*.05
        YMIN=YMIN-DY*.05
        YMAX=YMAX+DY*.05
5       CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,0)
        CALL PGIDEN_RED
ccc10      CALL PGSCI(NCOLOR)
        CALL PGSCI(NCOLOR)
        CALL PGPOINT(N,X,Y,NSYMB)
        CALL PGSCI(1)
        IF(.NOT.LLIMIT) RETURN
C
        IF(ITERM.EQ.1)THEN
          WRITE(*,100)'Change plot-limits (y/n) '
          CCHANGE(1:1)=READC('n','yn')
          IF(CCHANGE.EQ.'y')THEN
            WRITE(CDUMMY,*) XMIN
            WRITE(*,100)'Xmin '
            XMIN=READF(CDUMMY)
            WRITE(CDUMMY,*) XMAX
            WRITE(*,100)'Xmax '
            XMAX=READF(CDUMMY)
            WRITE(CDUMMY,*) YMIN
            WRITE(*,100)'Ymin '
            YMIN=READF(CDUMMY)
            WRITE(CDUMMY,*) YMAX
            WRITE(*,100)'Ymax '
            YMAX=READF(CDUMMY)
            GOTO 5 
          END IF
        END IF
C
100     FORMAT(A,$)
101     FORMAT(A)
        END
