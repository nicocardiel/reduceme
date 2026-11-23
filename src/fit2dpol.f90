!------------------------------------------------------------------------------
! Version 24-November-2026                                   file: fit2dpol.f90
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
! Program: fit2dpol
! Classification: arithmetic & manipulations
! Description: Fits 2-D polynomials.
!
!Comment
! Images shared by the main program and its procedures (they were in
! COMMON /BLK2/ and /BLK3/)
        MODULE FIT2DPOL_DATA_M
        IMPLICIT NONE
        REAL, ALLOCATABLE :: SB(:,:),SBB(:,:)
        END MODULE FIT2DPOL_DATA_M
!
! Ajusta un polinomio bidimensional a una imagen, por minimos
! cuadrados (previamente, se realiza un binning en la imagen).
!
        PROGRAM FIT2DPOL
        USE RED_M
        USE FSPEC_M
        USE FUTILS_M
        USE REDLIB_M
        USE FIT2DPOL_DATA_M
!
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
!
        INTEGER MAXNCOEFF
        PARAMETER(MAXNCOEFF=100)
!
        INTEGER NBINX,NBINY
        INTEGER NS1B,NS2B,NC1B,NC2B,NSCANB,NCHANB
        INTEGER NX,NY
        INTEGER NS1,NS2,NC1,NC2
        INTEGER I,J,K,L,M,P,Q
        INTEGER II,JJ
        INTEGER GX,GY,NCOEFF
        INTEGER ORDER(MAXNCOEFF),IOK,IPAR
        INTEGER NSMAX_LOCAL,NCMAX_LOCAL
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        REAL SUM
        REAL, ALLOCATABLE :: S(:,:)
        REAL X(NCMAX),Y(NSMAX)
        REAL XESC,YESC
        REAL XX,YY
        REAL A(MAXNCOEFF,MAXNCOEFF),B(MAXNCOEFF)
        REAL SCALEROW(MAXNCOEFF),XSOL(MAXNCOEFF)
        REAL FFACTOR
        REAL YMIN,YMAX,OFF
        REAL DEV,MEAN
        CHARACTER*1 CSAVE,CWHOLE,COTHER,CALL
        CHARACTER*255 INFILE,OUTFILE
        CHARACTER*80 GLABEL
        LOGICAL REPLOT
        LOGICAL LCOLOR(MAX_ID_RED)
!
        COMMON/BLK1/NX,NY
        COMMON/BLKSAVE/YMIN,YMAX,OFF
!-------------------------------------------------------------------------
        THISPROGRAM='fit2dpol'
        CALL WELCOME('24-November-2026')
!
        NSMAX_LOCAL=NSMAX
        NCMAX_LOCAL=NCMAX
        IF(NSMAX_LOCAL.GT.NCMAX_LOCAL)STOP 'FATAL ERROR: NSMAX.GT.NCMAX'
!
        REPLOT=.FALSE.
!
        CALL PIDEGTER(NTERM,IDN,LCOLOR)
!
! leemos imagen
        WRITE(*,100,ADVANCE='NO')'Input file name'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        ! allocate with the actual image dimensions (static arrays
        ! were limited by NCMAX,NSMAX); initialised as before; SB and SBB
        ! store the binned image in SB(1:NX,1:NY) and the whole image later
        ALLOCATE(S(NCHAN,NSCAN),SB(NCHAN,NSCAN),SBB(NCHAN,NSCAN))
        S=0.0
        SB=0.0
        SBB=0.0
        IF(TRUELEN(OBJECT).GT.0)THEN
          INFILE=INFILE(1:TRUELEN(INFILE))//' ['//&
     &     OBJECT(1:TRUELEN(OBJECT))//']'
        END IF
        DO I=1,NSCAN
          READ(20) (S(J,I),J=1,NCHAN)
        END DO
        CLOSE(20)
!
! pregunta region que va a ajustarse
10      WRITE(*,100,ADVANCE='NO')&
     & 'Are you fitting the whole image (y/n) '
        CWHOLE(1:1)=READC('y','yn')
        IF(CWHOLE.EQ.'y')THEN
          NS1B=1
          NS2B=NSCAN
          NSCANB=NSCAN
          NC1B=1
          NC2B=NCHAN
          NCHANB=NCHAN
        ELSE
          WRITE(*,101)'Define the region to be fitted (must be '//&
     &     ' a rectangle): '
          WRITE(*,100,ADVANCE='NO')'First and last channel'
          CALL READ2I('@',NC1B,NC2B)
          IF((NC1B.LT.1).OR.(NC2B.GT.NCHAN).OR.(NC1B.GT.NC2B))THEN
            WRITE(*,101)'ERROR: invalid numbers.'
            GOTO 10
          END IF
          NCHANB=NC2B-NC1B+1
          WRITE(*,100,ADVANCE='NO')'First and last scan   '
          CALL READ2I('@',NS1B,NS2B)
          IF((NS1B.LT.1).OR.(NS2B.GT.NSCAN).OR.(NS1B.GT.NS2B))THEN
            WRITE(*,101)'ERROR: invalid numbers.'
            GOTO 10
          END IF
          NSCANB=NS2B-NS1B+1
        END IF
! definimos el binning en la direccion X e Y
        WRITE(*,101)'INTRODUCE BINNING: '
!cc11      WRITE(*,100)'Binning in X-direction '
        WRITE(*,100,ADVANCE='NO')'Binning in X-direction '
        NBINX=READILIM('@',1,NCHANB)
!cc12      WRITE(*,100)'Binniny in Y-direction '
        WRITE(*,100,ADVANCE='NO')'Binniny in Y-direction '
        NBINY=READILIM('@',1,NSCANB)
!
! Calculamos puntos centrales en el grid definido por NBINX, NBINY.
! Obtenemos un grid nuevo de NX * NY elementos que almacenamos en SB(,)
        NX=NCHANB/NBINX
        IF(MOD(NCHANB,NBINX).NE.0) NX=NX+1
        NY=NSCANB/NBINY
        IF(MOD(NSCANB,NBINY).NE.0) NY=NY+1
        WRITE(*,'(A,I5,1X,I5)')'Dimensions of binned image: ',NX,NY
!
! definimos coordenadas X,Y de la imagen con binning de tal forma que
! su recorrido este en el intervalo (1,2]
        DO I=1,NY
          NS1=(I-1)*NBINY+NS1B
          NS2=NS1+NBINY
          IF(NS2.GT.NS2B) NS2=NS2B
          Y(I)=REAL(NS2+NS1)/2.
        END DO
        YESC=Y(NY)
        DO I=1,NY
          Y(I)=Y(I)/YESC+1.
        END DO
        DO J=1,NX
          NC1=(J-1)*NBINX+NC1B
          NC2=NC1+NBINX
          IF(NC2.GT.NC2B) NC2=NC2B
          X(J)=REAL(NC2+NC1)/2.
        END DO
        XESC=X(NX)
        DO J=1,NX
          X(J)=X(J)/XESC+1.
        END DO
!
        MEAN=0.
20      DO I=1,NY
          DO J=1,NX
            K=0
            SUM=0.
            NS1=(I-1)*NBINY+NS1B
            NS2=NS1+NBINY
            IF(NS2.GT.NS2B) NS2=NS2B
            DO II=NS1,NS2
              NC1=(J-1)*NBINX+NC1B
              NC2=NC1+NBINX
              IF(NC2.GT.NC2B) NC2=NC2B
              DO JJ=NC1,NC2
                K=K+1
                SUM=SUM+S(JJ,II)
              END DO
            END DO
            SUM=SUM/REAL(K)
            SB(J,I)=SUM
            MEAN=MEAN+SB(J,I)
          END DO
        END DO
        MEAN=MEAN/REAL((NX+1)*(NY+1))
        DEV=0.
        DO I=1,NY
          DO J=1,NX
            DEV=DEV+ABS(SB(J,I)-MEAN)
          END DO
        END DO
        DEV=DEV/REAL((NX+1)*(NY+1))
        WRITE(*,100,ADVANCE='NO')'Mean value in binned image        : '
        WRITE(*,*)MEAN
        WRITE(*,100,ADVANCE='NO')'Averaged deviation around the mean: '
        WRITE(*,*)DEV
        IF(REPLOT) GOTO 34
!
! dibujamos superficie con binning
!cc30      WRITE(*,100)'All the plots in the same page (y/n) '
        WRITE(*,100,ADVANCE='NO')'All the plots in the same page (y/n) '
        CALL(1:1)=READC('y','yn')
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(CALL.EQ.'y')THEN
            CALL PGSUBP(2,3)
          ELSE
            CALL PGSUBP(1,1)
          END IF
        END DO
34      DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM)) CALL PGSCI(2)
          CALL PLOTSB(.FALSE.,.TRUE.,LCOLOR(ITERM))
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
          CALL PGLABEL(CHAR(32),CHAR(32),INFILE)
        END DO
!
! grado de los polinomios en X,Y
40      WRITE(*,100,ADVANCE='NO')'Polynomial degree in X-direction '
        GX=READILIM('@',0,19)
        WRITE(*,100,ADVANCE='NO')'Polynomial degree in Y-direction '
        GY=READILIM('@',0,19)
        NCOEFF=(GX+1)*(GY+1)
        IF(NCOEFF.GT.MAXNCOEFF)THEN
          WRITE(*,101)'ERROR: polynomial degrees too large. Try '//&
     &     'with lower figures.'
          GOTO 40
        END IF
!
        WRITE(*,'(A,I3,A,I3,A)')'---> Solving ',NCOEFF,&
     &   ' equations with ',NCOEFF,' unknowns...'
!
! definimos sistema de (GX+1)*(GY+1) ecuaciones con (GX+1)*(GY+1) incognitas
! II es el numero de ecuacion, y JJ el numero de incognita, por lo que
! el sistema de ecuaciones es A(II,JJ) * X = B(II)
        II=0
        DO P=0,GX
          DO Q=0,GY
            II=II+1
            JJ=0
            DO I=0,GX
              DO J=0,GY
                JJ=JJ+1
                A(II,JJ)=0.
                DO L=1,NX
                  DO M=1,NY
                    FFACTOR=1.
                    IF(I+P.NE.0) FFACTOR=FFACTOR*(X(L)**(I+P))
                    IF(J+Q.NE.0) FFACTOR=FFACTOR*(Y(M)**(J+Q))
                    A(II,JJ)=A(II,JJ)+FFACTOR
!cc                 A(II,JJ)=A(II,JJ)+(X(L)**(I+P))*(Y(M)**(J+Q))
                  END DO
                END DO
              END DO
            END DO
            B(II)=0.
            DO L=1,NX
              DO M=1,NY
                FFACTOR=SB(L,M)
                IF(P.NE.0) FFACTOR=FFACTOR*(X(L)**(P))
                IF(Q.NE.0) FFACTOR=FFACTOR*(Y(M)**(Q))
                B(II)=B(II)+FFACTOR
!cc             B(II)=B(II)+SB(L,M)*(X(L)**(P))*(Y(M)**(Q))
              END DO
            END DO
          END DO
        END DO
! resolvemos el sistema de ecuaciones
        WRITE(*,100,ADVANCE='NO')'---> LU Descomposition...'
        CALL LUDCMP(A,NCOEFF,MAXNCOEFF,ORDER,SCALEROW,IOK,IPAR)
        WRITE(*,101)' OK!'
        WRITE(*,100,ADVANCE='NO')&
     &     '---> Forward substitution and backsubstitution...'
        CALL LUSOLV(A,NCOEFF,MAXNCOEFF,ORDER,SCALEROW,B,XSOL)
        WRITE(*,101)' OK!'
!
! dibujamos la imagen ajustada
        DO M=1,NY
          DO L=1,NX
            SBB(L,M)=SB(L,M)
            SB(L,M)=0.
            K=0
            DO I=0,GX
              DO J=0,GY
                K=K+1
                FFACTOR=XSOL(K)
                IF(I.NE.0) FFACTOR=FFACTOR*(X(L)**(I))
                IF(J.NE.0) FFACTOR=FFACTOR*(Y(M)**(J))
                SB(L,M)=SB(L,M)+FFACTOR
!cc                SB(L,M)=SB(L,M)+XSOL(K)*(X(L)**(I))*(Y(M)**(J))
              END DO
            END DO
          END DO
        END DO
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM)) CALL PGSCI(3)
          CALL PLOTSB(.FALSE.,.FALSE.,LCOLOR(ITERM))
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
          CALL PGLABEL(CHAR(32),CHAR(32),'polynomial fit')
          IF(LCOLOR(ITERM))THEN
            CALL PLOTCUTS(.TRUE.,.TRUE.,ITERM)
          ELSE
            CALL PLOTCUTS(.FALSE.,.TRUE.,ITERM)
          END IF
        END DO
        DO M=1,NY
          DO L=1,NX
            SB(L,M)=SBB(L,M)-SB(L,M)
          END DO
        END DO
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM)) CALL PGSCI(4)
          CALL PLOTSB(.FALSE.,.TRUE.,LCOLOR(ITERM))
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
          CALL PGLABEL('residual image',CHAR(32),CHAR(32))
          IF(CALL.EQ.'y')THEN
            CALL PGENV(0.,1.,0.,5.,0,-2)
            CALL PGIDEN_RED
            WRITE(GLABEL,'(A)')'Initial file name:'
            CALL PGPTEXT(.5,1.6,0.,1.,GLABEL)
            CALL PGPTEXT(.6,1.6,0.,0.,INFILE)
            WRITE(GLABEL,'(A)')'Image size (NCHAN,NSCAN): '
            CALL PGPTEXT(.5,1.2,0.,1.,GLABEL)
            WRITE(GLABEL,*)NCHAN,NSCAN
            CALL PGPTEXT(.5,1.2,0.,0.,GLABEL)
            WRITE(GLABEL,'(A)')'First and last channel (binning):'
            CALL PGPTEXT(.5,1.,0.,1.,GLABEL)
            WRITE(GLABEL,*)NC1B,NC2B
            CALL PGPTEXT(.5,1.,0.,0.,GLABEL)
            WRITE(GLABEL,'(A)')'First and last scan (binning): '
            CALL PGPTEXT(.5,.8,0.,1.,GLABEL)
            WRITE(GLABEL,*)NS1B,NS2B
            CALL PGPTEXT(.5,.8,0.,0.,GLABEL)
            WRITE(GLABEL,'(A)')'Binning X-direction: '
            CALL PGPTEXT(.5,.6,0.,1.,GLABEL)
            WRITE(GLABEL,*)NBINX
            CALL PGPTEXT(.5,.6,0.,0.,GLABEL)
            WRITE(GLABEL,'(A)')'Binning Y-direction: '
            CALL PGPTEXT(.5,.4,0.,1.,GLABEL)
            WRITE(GLABEL,*)NBINY
            CALL PGPTEXT(.5,.4,0.,0.,GLABEL)
            WRITE(GLABEL,'(A)')'Pol. degree X-direction: '
            CALL PGPTEXT(.5,.2,0.,1.,GLABEL)
            WRITE(GLABEL,*)GX
            CALL PGPTEXT(.5,.2,0.,0.,GLABEL)
            WRITE(GLABEL,'(A)')'Pol. degree Y-direction: '
            CALL PGPTEXT(.5,.0,0.,1.,GLABEL)
            WRITE(GLABEL,*)GY
            CALL PGPTEXT(.5,.0,0.,0.,GLABEL)
          END IF
        END DO
! superficie ajustada pero sin binning
        DO M=1,NSCAN
          DO L=1,NCHAN
            SB(L,M)=0.
            K=0
            DO I=0,GX
              DO J=0,GY
                K=K+1
                XX=REAL(L)/XESC+1.
                YY=REAL(M)/YESC+1.
                FFACTOR=XSOL(K)
                IF(I.NE.0) FFACTOR=FFACTOR*(XX**(I))
                IF(J.NE.0) FFACTOR=FFACTOR*(YY**(J))
                SB(L,M)=SB(L,M)+FFACTOR
!cc             SB(L,M)=SB(L,M)+XSOL(K)*(XX**(I))*(YY**(J))
              END DO
            END DO
          END DO
        END DO
!
! otro ajuste o salvar/salir
        WRITE(*,100,ADVANCE='NO')'Other fit (y/n) '
        COTHER(1:1)=READC('y','yn')
        IF(COTHER.EQ.'y')THEN
          REPLOT=.TRUE.
          GOTO 20
        END IF
! salvando imagenes
        WRITE(*,100,ADVANCE='NO')'Save fitted image -not binned- (y/n) '
        CSAVE(1:1)=READC('n','yn')
        IF(CSAVE.EQ.'y')THEN
          WRITE(*,100,ADVANCE='NO')'Output file name'
          OUTFILE=OUTFILEX(30,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
            DO I=1,NSCAN
              WRITE(30) (SB(J,I),J=1,NCHAN)
            END DO
          CALL CLOSEFILEX(30)
        END IF
!
        WRITE(*,100,ADVANCE='NO')&
     &     'Save residual image -not binned- (y/n) '
        CSAVE(1:1)=READC('n','yn')
        IF(CSAVE.EQ.'y')THEN
          WRITE(*,100,ADVANCE='NO')'Output file name'
          OUTFILE=OUTFILEX(30,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
            DO I=1,NSCAN
              DO J=1,NCHAN
                S(J,I)=S(J,I)-SB(J,I)
              END DO
            WRITE(30) (S(J,I),J=1,NCHAN)
          END DO
          CALL CLOSEFILEX(30)
        END IF
!
        CALL PGEND
        STOP
!
100     FORMAT(A)
101     FORMAT(A)
        END
!
!*************************************************************************
!
! dibuja la superficie de puntos
        SUBROUTINE PLOTSB(OVERP,LLIMITS,LCOLOR)
        USE RED_M
        USE FUTILS_M
        USE REDLIB_M
        USE FIT2DPOL_DATA_M
        IMPLICIT NONE
        LOGICAL OVERP,LLIMITS
        LOGICAL LCOLOR
        INCLUDE 'redlib.inc'
!
        INTEGER I,J
        INTEGER NX,NY
        INTEGER LASTCOLOR
! XP,YP hay que dimensionarlos al mayor de NCMAX,NSMAX
        REAL XP(NCMAX),YP(NCMAX)
        REAL YMIN,YMAX,OFF
        REAL FACTOR
!
        COMMON/BLK1/NX,NY
        COMMON/BLKSAVE/YMIN,YMAX,OFF
!
!-------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
! FACTOR determina como aparece de ampliada la estructura en el grafico
! que simula una representacion tridimensional de la superficie
        FACTOR=4.
!
        IF(.NOT.OVERP)THEN
          IF(LLIMITS)THEN
            YMIN=SB(1,1)
            YMAX=YMIN
            DO I=1,NY
              DO J=1,NX
                IF(SB(J,I).LT.YMIN) YMIN=SB(J,I)
                IF(SB(J,I).GT.YMAX) YMAX=SB(J,I)
              END DO
            END DO
            OFF=YMAX-YMIN
            YMIN=YMIN-OFF*.05
          END IF
          IF(LCOLOR)THEN
            CALL PGQCI(LASTCOLOR)
            CALL PGSCI(1)
          END IF
          CALL PGENV(0.,REAL(NX)*4./3.+1,YMIN,YMAX+FACTOR*OFF,0,-2)
          CALL PGBOX(' ',0.0,0,'BTNS',0.0,0)
          CALL PGIDEN_RED
          IF(LCOLOR) CALL PGSCI(LASTCOLOR)
        END IF
!
        DO I=1,NY
          DO J=1,NX
            IF(NY.EQ.1)THEN
              XP(J)=REAL(J)
              YP(J)=SB(J,I)
            ELSE
              XP(J)=REAL(J)+REAL(I-1)/REAL(NY-1)*REAL(NX)/3.
              YP(J)=SB(J,I)+REAL(I-1)/REAL(NY-1)*FACTOR*OFF
            END IF
          END DO
          CALL PGLINE(NX,XP,YP)
        END DO
        IF(NY.GT.1)THEN
          DO J=1,NX
            DO I=1,NY
              XP(I)=REAL(J)+REAL(I-1)/REAL(NY-1)*REAL(NX)/3.
              YP(I)=SB(J,I)+REAL(I-1)/REAL(NY-1)*FACTOR*OFF
            END DO
            CALL PGLINE(NY,XP,YP)
          END DO
        END IF
        END
!
!*************************************************************************
!
! dibuja los ajustes en la direccion espacial y espectral
        SUBROUTINE PLOTCUTS(IFCOLOR,LPLOTS,ITERM)
        USE RED_M
        USE FUTILS_M
        USE REDLIB_M
        USE FIT2DPOL_DATA_M
        IMPLICIT NONE
        LOGICAL IFCOLOR,LPLOTS
        INTEGER ITERM
        INCLUDE 'redlib.inc'
!
        INTEGER I,J
        INTEGER NX,NY
! XP,YP,DEV hay que dimensionarlos al mayor de NCMAX,NSMAX
        REAL XP(NCMAX),YP(NCMAX)
        REAL YMIN,YMAX,DY
        REAL DEV(NCMAX),MDEV
!
        COMMON/BLK1/NX,NY
!
!-------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
        DY=0.0 !avoid compilation warning
!
        IF(LPLOTS)THEN
           YMIN=SB(1,1)
          YMAX=YMIN
          DO I=1,NY
            DO J=1,NX
              IF(SB(J,I).LT.YMIN) YMIN=SB(J,I)
              IF(SBB(J,I).LT.YMIN) YMIN=SBB(J,I)
              IF(SB(J,I).GT.YMAX) YMAX=SB(J,I)
              IF(SBB(J,I).GT.YMAX) YMAX=SBB(J,I)
            END DO
          END DO
          DY=YMAX-YMIN
          YMIN=YMIN-0.05*DY
          YMAX=YMAX+0.05*DY
          DY=YMAX-YMIN
        END IF
!
! cortes direccion espectral
        IF(LPLOTS)THEN
          CALL PGSCI(1)
          CALL PGENV(0.,REAL(NX)+1.,YMIN,YMIN+REAL(NY)*DY,0,0)
          CALL PGIDEN_RED
          CALL PGLABEL('X-direction',CHAR(32),'Binned & Fitted Image')
        END IF
        DO I=1,NY
          IF(LPLOTS)THEN
            DO J=1,NX
              XP(J)=REAL(J)
              YP(J)=SBB(J,I)+REAL(I-1)*DY
            END DO
            IF(IFCOLOR) CALL PGSCI(2)
            CALL PGSLS(4)
            CALL PGLINE(NX,XP,YP)
            CALL PGPOINT(NX,XP,YP,17)
            CALL PGSLS(1)
          END IF
          DEV(I)=0.
          DO J=1,NX
            IF(LPLOTS) YP(J)=SB(J,I)+REAL(I-1)*DY
            DEV(I)=DEV(I)+ABS(SB(J,I)-SBB(J,I))
          END DO
          IF(LPLOTS)THEN
            IF(IFCOLOR) CALL PGSCI(3)
            CALL PGLINE(NX,XP,YP)
          END IF
          DEV(I)=DEV(I)/REAL(NX)
        END DO
        MDEV=0.
        DO I=1,NY
          MDEV=MDEV+DEV(I)
        END DO
        MDEV=MDEV/REAL(NY)
        IF(ITERM.EQ.1)THEN
          WRITE(*,'(A,I4,A)',ADVANCE='NO')'Mean deviation over ',NY,&
     &     '  scans   : '
          WRITE(*,*)MDEV
        END IF
!
! cortes en la direccion espectral
        IF(LPLOTS)THEN
          CALL PGSCI(1)
          CALL PGENV(0.,REAL(NY)+1.,YMIN,YMIN+REAL(NX)*DY,0,0)
          CALL PGIDEN_RED
          CALL PGLABEL('Y-direction',CHAR(32),'Binned & Fitted Image')
        END IF
        DO J=1,NX
          IF(LPLOTS)THEN
            DO I=1,NY
              XP(I)=REAL(I)
              YP(I)=SBB(J,I)+REAL(J-1)*DY
            END DO
            IF(IFCOLOR) CALL PGSCI(2)
            CALL PGSLS(4)
            CALL PGLINE(NY,XP,YP)
            CALL PGPOINT(NY,XP,YP,17)
            CALL PGSLS(1)
          END IF
          DEV(J)=0.
          DO I=1,NY
            IF(LPLOTS) YP(I)=SB(J,I)+REAL(J-1)*DY
            DEV(J)=DEV(J)+ABS(SB(J,I)-SBB(J,I))
          END DO
          IF(LPLOTS)THEN
            IF(IFCOLOR) CALL PGSCI(3)
            CALL PGLINE(NY,XP,YP)
          END IF
          DEV(J)=DEV(J)/REAL(NY)
        END DO
        MDEV=0.
        DO J=1,NX
          MDEV=MDEV+DEV(J)
        END DO
        MDEV=MDEV/REAL(NX)
        IF(ITERM.EQ.1)THEN
          WRITE(*,'(A,I4,A)',ADVANCE='NO')'Mean deviation over ',NX,&
     &     '  channels: '
          WRITE(*,*)MDEV
        END IF
!
        END
