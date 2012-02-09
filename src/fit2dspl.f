C------------------------------------------------------------------------------
C Version 19-June-1998                                         file: fit2dspl.f
C @ ncl & fjg
C------------------------------------------------------------------------------
Comment
C
C Program: fit2dspl
C Classification: arithmetic & manipulations
C Description: Fits 2-D splines and local polynomials.
C
Comment
C
C Ajusta una superficie cualquiera mediante superficies polinomicas mas
C pequen~as. Tambien permite ajustar mediante "bicubic splines". El programa 
C ajusta por splines monodimensionales el corte espacial promedio y el
C corte en longitud de onda promedio. Estos ajustes son refinados por la
C tecnica usada en splfit para buscar las mejores posiciones de los knots.
C Las posiciones anteriores seran usadas (invariantemente) para calcular la
C superficie bidimensional. 
C
C Se realiza un ajuste inicial bidimensional a superficies polinomicas de bajo 
C grado alrededor de cada Knot. Se utilizan estos ajustes para determinar el
C valor de la imagen en cada Knot, asi como generar una imagen (FINAL) con los
C valores locales de las superficies ajustadas. 
C
C NOTA: Si la imagen a ajustar tiene valores con diferencias muy pequenhas,
C es conveniente restar a la imagen el valor medio antes del ajuste (luego el
C valor medio vuelve a sumarse al final), y tambien posiblemente multiplicar
C la imagen por un factor de escala, para conseguir un mejor rango en los
C calculos.
C
        PROGRAM FIT2DSPL
C
        IMPLICIT NONE
C
C Include NSMAX,NCMAX,NSCAN,NCHAN
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER TRUELEN
        INTEGER READI
        INTEGER READILIM
        REAL READF
C
        INTEGER NMAXKNOT                    !numero maximo de Knots en cada eje
        PARAMETER(NMAXKNOT=30)                 !modificar tambien en subrutinas
C
        INTEGER I,J,K,L
        INTEGER NPLOT
        INTEGER KK1,KK2,NPTKNOT
        INTEGER NS1CX,NS2CX,NC1CY,NC2CY,NSCAN_ADDED,NCHAN_ADDED
        INTEGER NKNOTX,NKNOTY,NF
        INTEGER KNOTX(NMAXKNOT),KNOTY(NMAXKNOT)
        INTEGER KREP
        INTEGER NBINX,NBINY,NBINX_OLD,NBINY_OLD
        INTEGER GX,GY,BINMODE,GX_OLD,GY_OLD,BINMODE_OLD
        INTEGER NITER,NKJUMP,PESOPOWER
        INTEGER NSMAX_LOCAL,NCMAX_LOCAL
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        REAL S(NCMAX,NSMAX),SB(NCMAX,NSMAX),SBFIT(NCMAX,NSMAX)
        REAL FINAL(NCMAX,NSMAX),FINAL_PESO(NCMAX,NSMAX)
        REAL CORTEX(NCMAX),CORTEY(NSMAX)
        REAL AA(NCMAX,NSMAX),BB(NCMAX,NSMAX),CC(NCMAX,NSMAX)
        REAL X(NCMAX),Y(NCMAX),XXX(NCMAX),YYY(NCMAX)
        REAL X0,Y0,Z0
        REAL XMIN,XMAX,YMIN,YMAX,DX,DY
        REAL YRMSTOL,SIGMAF
        REAL FKNOTX(NMAXKNOT),FKNOTY(NMAXKNOT)
        REAL FCTE1,FCTE2
        REAL FG,BG,TR(6)
        REAL XMIN_IMA,XMAX_IMA,YMIN_IMA,YMAX_IMA
        REAL POLYSURF,FX,FY,BINSIGMA,BINSIGMA_OLD
        REAL SIGMAREJ
        DOUBLE PRECISION DMEAN,DSIGMA
        CHARACTER*1 CMASK,CSAVE,CREP,CSUB,CYLIMITS,CNEXT,CRESTART
        CHARACTER*1 CPOK
        CHARACTER*50 CDUMMY,GLABEL
        CHARACTER*75 INFILE,MASKFILE,OUTFILE
        LOGICAL LFITX(NCMAX),LFITY(NSMAX)
        LOGICAL LMASK(NCMAX,NSMAX)
        LOGICAL IFSCAN(NSMAX),IFCHAN(NCMAX)
        LOGICAL LREP
        LOGICAL LCOLOR(MAX_ID_RED)
C
        COMMON/BLKPIDE1/NSCAN,NCHAN
        COMMON/BLKPIDE2/CMASK
        COMMON/BLKYLIMITS/YMIN,YMAX
        COMMON/BLKPOLYNOM1/S,SB
        COMMON/BLKPOLYNOM2/NBINX,NBINY,BINMODE
        COMMON/BLKPOLYNOM3/LFITX,LFITY
        COMMON/BLKPOLYNOM3BIS/LMASK
        COMMON/BLKPOLYNOM4/TR,BINSIGMA
        COMMON/BLKPOLYNOM5/SIGMAREJ,FX,FY,FKNOTX,FKNOTY
        COMMON/BLKPOLYNOM6/NITER,PESOPOWER
        COMMON/BLKPOLYNOM7/FINAL,FINAL_PESO
        COMMON/BLKPOLYNOM8/KNOTX,KNOTY,NKNOTX,NKNOTY
        COMMON/BLKPOLYNOM9/CNEXT
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
        THISPROGRAM='fit2dspl'
        CALL WELCOME('6-December-1996')
C
        NSMAX_LOCAL=NSMAX
        NCMAX_LOCAL=NCMAX
        IF(NSMAX_LOCAL.GT.NCMAX_LOCAL)STOP 'FATAL ERROR: NSMAX.GT.NCMAX'
C
        TR(1)=0.
        TR(2)=1.
        TR(3)=0.
        TR(4)=0.
        TR(5)=0.
        TR(6)=1.
C
        BINSIGMA=3.0
        YRMSTOL=1.0E-6
        PESOPOWER=2
C Evitamos warnings de compilacion
        KREP=0
        NKJUMP=0
C Abrimos salida grafica
        CALL PIDEGTER(NTERM,IDN,LCOLOR)
C------------------------------------------------------------------------------
C Leemos imagen
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
C Permitimos usar una máscara. En este caso, para facilitar la introducción de
C este cambio en el código, no permitimos ni binning ni eliminar filas o
C columnas
        WRITE(*,100)'Are you using a mask file (y/n) '
        CMASK(1:1)=READC('n','yn')
        IF(CMASK.EQ.'y')THEN
          WRITE(*,101)'> NOTE: Pixels with signal != 0.0 will be masked'
          WRITE(*,100)'Input mask file name'
          MASKFILE=INFILEX(25,'@',NSCAN,NCHAN,STWV,DISP,21,.FALSE.)
          DO I=1,NSCAN
            READ(25) (SB(J,I),J=1,NCHAN)  !usamos temporalmente el array SB
          END DO
          CLOSE(25)
          DO I=1,NSCAN
            DO J=1,NCHAN
              LMASK(J,I)=(SB(J,I).EQ.0.0)
            END DO
          END DO
        ELSE
          DO I=1,NSCAN
            DO J=1,NCHAN
              LMASK(J,I)=.TRUE. !ningún píxel estará enmascarado
            END DO
          END DO
        END IF
C------------------------------------------------------------------------------
        FCTE1=0.0
        FCTE2=1.0
        WRITE(*,100)'Are you using (data-cte)*factor (y/n) '
        CSUB(1:1)=READC('n','yn')
        IF(CSUB.EQ.'y')THEN
          WRITE(*,100)'Cte...'
          WRITE(CDUMMY,*) FCTE1
          FCTE1=READF(CDUMMY)
          WRITE(*,100)'Factor'
          WRITE(CDUMMY,*) FCTE2
          FCTE2=READF(CDUMMY)
          DO I=1,NSCAN
            DO J=1,NCHAN
              S(J,I)=(S(J,I)-FCTE1)*FCTE2
            END DO
          END DO
        END IF
C------------------------------------------------------------------------------
C Ajustamos corte en direccion X (longitud de onda)
10      WRITE(*,*)
        WRITE(*,101)'>>> KNOTS in the X-direction:'
        IF(CMASK.EQ.'y')THEN
          DO I=1,NSCAN
            IFSCAN(I)=.TRUE.
          END DO
          NSCAN_ADDED=NSCAN
        ELSE
          DO I=1,NSCAN
            IFSCAN(I)=.FALSE.
          END DO
          WRITE(*,101)'Enter scan region to obtain averaged spectrum:'
          WRITE(CDUMMY,'(A,I6)')'1,',NSCAN
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          WRITE(*,101)'Valid region is: '//CDUMMY(1:L)
11        WRITE(*,100)'1st and last scan (0,0=EXIT) '
          CALL READ2I('0,0',NS1CX,NS2CX)
          IF((NS1CX.EQ.0).AND.(NS2CX.EQ.0)) GOTO 12
          IF((NS1CX.LT.1).OR.(NS2CX.GT.NSCAN).OR.(NS1CX.GT.NS2CX))THEN
            WRITE(*,101)'ERROR: numbers out of range. Try again.'
            GOTO 11
          END IF
          DO I=NS1CX,NS2CX
            IFSCAN(I)=.TRUE.
          END DO
          GOTO 11
C
12        NSCAN_ADDED=0
          DO I=1,NSCAN
            IF(IFSCAN(I)) NSCAN_ADDED=NSCAN_ADDED+1
          END DO
          IF(NSCAN_ADDED.EQ.0)THEN
            WRITE(*,101)'ERROR: no. of scans added = 0. Try again.'
            GOTO 11
          END IF
        END IF
C
        DO J=1,NCHAN
          CORTEX(J)=0.
          X(J)=REAL(J)
        END DO
        DO I=1,NSCAN
          IF(IFSCAN(I))THEN
            DO J=1,NCHAN
              CORTEX(J)=CORTEX(J)+S(J,I)
            END DO
          END IF
        END DO
        DO J=1,NCHAN
          CORTEX(J)=CORTEX(J)/REAL(NSCAN_ADDED)
        END DO
        XMIN=1.
        XMAX=REAL(NCHAN)
        CALL FINDMM(NCHAN,CORTEX,YMIN,YMAX)
        DX=XMAX-XMIN
        DY=YMAX-YMIN
        XMIN=XMIN-DX/50.
        XMAX=XMAX+DX/50.
        YMIN=YMIN-DY/50.
        YMAX=YMAX+DY/50.
13        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,0)
          CALL PGIDEN_RED
          CALL PGLABEL('channel','Averaged no. of counts','File: '//
     +     INFILE)
          IF(LCOLOR(ITERM)) CALL PGSCI(2)
          CALL PGBIN(NCHAN,X,CORTEX,.TRUE.)
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
        END DO
        WRITE(*,100)'Change Y-limits (y/n) '
        CYLIMITS(1:1)=READC('n','yn')
        IF(CYLIMITS.EQ.'y')THEN
          WRITE(CDUMMY,*) YMIN
          WRITE(*,100)'New Ymin value '
          YMIN=READF(CDUMMY)
          WRITE(CDUMMY,*) YMAX
          WRITE(*,100)'New Ymax value '
          YMAX=READF(CDUMMY)
          GOTO 13
        END IF
C
        CALL PIDEKNOT(NCHAN,LFITX,CORTEX,NKNOTX,KNOTX,YRMSTOL)
        NF=0
        DO J=1,NCHAN
          IF(LFITX(J))THEN
            NF=NF+1
            X(NF)=REAL(J)
            Y(NF)=CORTEX(J)
          END IF
        END DO
        DO K=1,NKNOTX
          FKNOTX(K)=REAL(KNOTX(K))
        END DO
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM)) CALL PGSCI(3)
        END DO
        CALL SPLFIT(NF,X,Y,NKNOTX,FKNOTX,YRMSTOL,NCHAN,XXX,YYY,
     +   1.,REAL(NCHAN),SIGMAF,.TRUE.)
C redondeamos las posiciones de los Knots a numeros enteros
        DO K=1,NKNOTX
          KNOTX(K)=NINT(FKNOTX(K))
        END DO
C eliminamos posibles knots repetidos
15      LREP=.FALSE.
        DO K=2,NKNOTX
          IF(KNOTX(K).EQ.KNOTX(K-1))THEN
            LREP=.TRUE.
            KREP=K
          END IF
        END DO
        IF(LREP)THEN
          WRITE(*,101)'Removing knot...'
          IF(KREP.LT.NKNOTX)THEN
            DO K=KREP,NKNOTX-1
              KNOTX(K)=KNOTX(K+1)
            END DO
          END IF
          NKNOTX=NKNOTX-1
          GOTO 15
        END IF
C comprobamos que hay datos para ajustar a ambos lados de cada knot 
C (salvo en los bordes, donde solo exigimos a un lado)
        WRITE(*,*)
        DO K=1,NKNOTX
          KK1=K-1
          KK2=K+1
          IF(KK1.LT.1)THEN
            KK1=KNOTX(1)
          ELSE
            KK1=KNOTX(KK1)+1
          END IF
          IF(KK2.GT.NKNOTX)THEN
            KK2=KNOTX(NKNOTX)
          ELSE
            KK2=KNOTX(KK2)-1
          END IF
          NPTKNOT=0
          DO J=KK1,KK2
            IF(LFITX(J)) NPTKNOT=NPTKNOT+1
          END DO
          WRITE(*,'(A,I2,A,I6)')'Knot #',K,'  located at: ',
     +     KNOTX(K)
          IF(NPTKNOT.LT.1)THEN
            WRITE(*,101)'ERROR: insuficient number of points.'
            WRITE(*,100)'(press <CR> to continue...)'
            READ(*,*)
            GOTO 10
          END IF
        END DO
C
        WRITE(*,*)
        WRITE(*,100)'Repeat fit (y/n) '
        CREP(1:1)=READC('n','yn')
        IF(CREP.EQ.'y') GOTO 10
        DO K=1,NKNOTX
          FKNOTX(K)=REAL(KNOTX(K))
        END DO
C------------------------------------------------------------------------------
        WRITE(*,100)'Save last fit (wavelength cross section) (y/n) '
        CSAVE(1:1)=READC('n','yn')
        IF(CSAVE.EQ.'y')THEN
          WRITE(*,100)'Output file name'
          OUTFILE=OUTFILEX(30,'@',1,NCHAN,STWV,DISP,1,.FALSE.)
          IF(CSUB.EQ.'y')THEN
            DO J=1,NCHAN
              YYY(J)=YYY(J)/FCTE2+FCTE1
            END DO
          END IF
          WRITE(30) (YYY(J),J=1,NCHAN)
          CLOSE(30)
        END IF
C------------------------------------------------------------------------------
C Ajustamos corte en direccion Y (direccion espacial)
20      WRITE(*,*)
        WRITE(*,101)'>>> KNOTS in the Y-direction:'
        IF(CMASK.EQ.'y')THEN
          DO J=1,NCHAN
            IFCHAN(J)=.TRUE.
          END DO
          NCHAN_ADDED=NCHAN
        ELSE
          DO J=1,NCHAN
            IFCHAN(J)=.FALSE.
          END DO
          WRITE(*,101)'Enter channel region to obtain averaged '//
     +     'spatial form:'
          WRITE(CDUMMY,'(A,I6)')'1,',NCHAN
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          WRITE(*,101)'Valid region is: '//CDUMMY(1:L)
21          WRITE(*,100)'1st and last channel (0,0=EXIT) '
          CALL READ2I('0,0',NC1CY,NC2CY)
          IF((NC1CY.EQ.0).AND.(NC2CY.EQ.0)) GOTO 22
          IF((NC1CY.LT.1).OR.(NC2CY.GT.NCHAN).OR.(NC1CY.GT.NC2CY))THEN
            WRITE(*,101)'ERROR: numbers out of range. Try again.'
            GOTO 21
          END IF
          DO J=NC1CY,NC2CY
            IFCHAN(J)=.TRUE.
          END DO
          GOTO 21
C
22        NCHAN_ADDED=0
          DO J=1,NCHAN
            IF(IFCHAN(J)) NCHAN_ADDED=NCHAN_ADDED+1
          END DO
          IF(NCHAN_ADDED.EQ.0)THEN
            WRITE(*,101)'ERROR: no. of channels added = 0. Try again.'
            GOTO 21
          END IF
        END IF
C
        DO I=1,NSCAN
          CORTEY(I)=0.
          Y(I)=REAL(I)
        END DO
        DO J=1,NCHAN
          IF(IFCHAN(J))THEN
            DO I=1,NSCAN
              CORTEY(I)=CORTEY(I)+S(J,I)
            END DO
          END IF
        END DO
        DO I=1,NSCAN
          CORTEY(I)=CORTEY(I)/REAL(NCHAN_ADDED)
        END DO
        XMIN=1.
        XMAX=REAL(NSCAN)
        CALL FINDMM(NSCAN,CORTEY,YMIN,YMAX)
        DX=XMAX-XMIN
        DY=YMAX-YMIN
        XMIN=XMIN-DX/50.
        XMAX=XMAX+DX/50.
        YMIN=YMIN-DY/50.
        YMAX=YMAX+DY/50.
23        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,0)
          CALL PGIDEN_RED
          CALL PGLABEL('scan','Averaged no. of counts','File: '//INFILE)
          IF(LCOLOR(ITERM)) CALL PGSCI(2)
          CALL PGBIN(NSCAN,Y,CORTEY,.TRUE.)
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
        END DO
        WRITE(*,100)'Change Y-limits (y/n) '
        CYLIMITS(1:1)=READC('n','yn')
        IF(CYLIMITS.EQ.'y')THEN
          WRITE(CDUMMY,*) YMIN
          WRITE(*,100)'New Ymin value '
          YMIN=READF(CDUMMY)
          WRITE(CDUMMY,*) YMAX
          WRITE(*,100)'New Ymax value '
          YMAX=READF(CDUMMY)
          GOTO 23
        END IF
C
        CALL PIDEKNOT(NSCAN,LFITY,CORTEY,NKNOTY,KNOTY,YRMSTOL)
        NF=0
        DO I=1,NSCAN
          IF(LFITY(I))THEN
            NF=NF+1
            X(NF)=REAL(I)
            Y(NF)=CORTEY(I)
          END IF
        END DO
        DO K=1,NKNOTY
          FKNOTY(K)=REAL(KNOTY(K))
        END DO
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM)) CALL PGSCI(3)
        END DO
        CALL SPLFIT(NF,X,Y,NKNOTY,FKNOTY,YRMSTOL,NSCAN,XXX,YYY,
     +   1.,REAL(NSCAN),SIGMAF,.TRUE.)
C redondeamos las posiciones de los Knots a numeros enteros
        DO K=1,NKNOTY
          KNOTY(K)=NINT(FKNOTY(K))
        END DO
C eliminamos posibles knots repetidos
25      LREP=.FALSE.
        DO K=2,NKNOTY
          IF(KNOTY(K).EQ.KNOTY(K-1))THEN
            LREP=.TRUE.
            KREP=K
          END IF
        END DO
        IF(LREP)THEN
          WRITE(*,101)'Removing knot...'
          IF(KREP.LT.NKNOTY)THEN
            DO K=KREP,NKNOTY-1
              KNOTY(K)=KNOTY(K+1)
            END DO
          END IF
          NKNOTY=NKNOTY-1
          GOTO 25
        END IF
C comprobamos que hay datos para ajustar a ambos lados de cada knot 
C (salvo en los bordes, donde solo exigimos a un lado)
        WRITE(*,*)
        DO K=1,NKNOTY
          KK1=K-1
          KK2=K+1
          IF(KK1.LT.1)THEN
            KK1=KNOTY(1)
          ELSE
            KK1=KNOTY(KK1)+1
          END IF
          IF(KK2.GT.NKNOTY)THEN
            KK2=KNOTY(NKNOTY)
          ELSE
            KK2=KNOTY(KK2)-1
          END IF
          NPTKNOT=0
          DO I=KK1,KK2
            IF(LFITY(I)) NPTKNOT=NPTKNOT+1
          END DO
          WRITE(*,'(A,I2,A,I6)')'Knot #',K,'  located at: ',
     +     KNOTY(K)
          IF(NPTKNOT.LT.1)THEN
            WRITE(*,101)'ERROR: insuficient number of points.'
            WRITE(*,100)'(press <CR> to continue...)'
            READ(*,*)
            GOTO 20
          END IF
        END DO
C
        WRITE(*,*)
        WRITE(*,100)'Repeat fit (y/n) '
        CREP(1:1)=READC('n','yn')
        IF(CREP.EQ.'y') GOTO 20
        DO K=1,NKNOTY
          FKNOTY(K)=REAL(KNOTY(K))
        END DO
C------------------------------------------------------------------------------
        WRITE(*,100)'Save last fit (spatial cross section) (y/n) '
        CSAVE(1:1)=READC('n','yn')
        IF(CSAVE.EQ.'y')THEN
          WRITE(*,100)'Output file name'
          OUTFILE=OUTFILEX(30,'@',NSCAN,1,STWV,DISP,1,.FALSE.)
          IF(CSUB.EQ.'y')THEN
            DO I=1,NSCAN
              YYY(I)=YYY(I)/FCTE2+FCTE1
            END DO
          END IF
          DO I=1,NSCAN
            WRITE(30) YYY(I)
          END DO
          CLOSE(30)
        END IF
C------------------------------------------------------------------------------
30        CONTINUE
C inicializamos la matriz FINAL
        DO I=1,NSCAN
          DO J=1,NCHAN
            FINAL(J,I)=0.
          END DO
        END DO
C inicializamos la matriz FINAL_PESO
        DO I=1,NSCAN
          DO J=1,NCHAN
            FINAL_PESO(J,I)=0.
          END DO
        END DO
C------------------------------------------------------------------------------
C dibujamos region que va a ser utilizada para realizar los ajustes
        DO I=1,NSCAN
          DO J=1,NCHAN
            SB(J,I)=S(J,I)
          END DO
        END DO
C media en toda la imagen que va a ser ajustada
        NPTKNOT=0
        DMEAN=0.D0
        DO I=1,NSCAN
          IF(LFITY(I))THEN
            DO J=1,NCHAN
              IF(LFITX(J))THEN
                IF(LMASK(J,I))THEN
                  DMEAN=DMEAN+DBLE(SB(J,I))
                  NPTKNOT=NPTKNOT+1
                END IF
              END IF
            END DO
          END IF
        END DO
        DMEAN=DMEAN/DBLE(NPTKNOT)
C desviacion estandard en toda la imagen que va a ser ajustada
        DSIGMA=0.D0
        IF(NPTKNOT.GT.1)THEN
          DO I=1,NSCAN
            IF(LFITY(I))THEN
              DO J=1,NCHAN
                IF(LFITX(J))THEN
                  IF(LMASK(J,I))THEN
                    DSIGMA=DSIGMA+
     +               (DBLE(SB(J,I))-DMEAN)*(DBLE(SB(J,I))-DMEAN)
                  END IF
                END IF
              END DO
            END IF
          END DO
          DSIGMA=DSQRT(DSIGMA/DBLE(NPTKNOT-1))
        ELSE
          DSIGMA=DMEAN                    !para que PGPLOT continue funcionando
        END IF
C mostramos la imagen con cortes a 3 sigma alrededor de la media
        BG=REAL(DMEAN-3.D0*DSIGMA)
        FG=REAL(DMEAN+3.D0*DSIGMA)
        XMIN_IMA=0.5
        XMAX_IMA=REAL(NCHAN)+0.5
        YMIN_IMA=0.5
        YMAX_IMA=REAL(NSCAN)+0.5
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGPAGE
          CALL PGVPORT(.07,.95,.50,.94)
          CALL PGWINDOW(XMIN_IMA,XMAX_IMA,YMIN_IMA,YMAX_IMA)
          CALL PGSCH(0.7)
          CALL PGLABEL('channel','scan','File: '//INFILE)
          CALL PGGRAY(SB,NCMAX,NSMAX,1,NCHAN,1,NSCAN,FG,BG,TR)
          CALL PGBOX('BCTNSI',0.0,0,'BCTNSI',0.0,0)
          CALL PGSCH(1.0)
          CALL PGIDEN_RED
        END DO
C
        WRITE(*,100)'Press <CR> to see image regions to be fitted...'
        READ(*,*)
C eliminamos de la imagen los puntos que no van a ser ajustados
        DO I=1,NSCAN
          IF(.NOT.LFITY(I))THEN
            DO J=1,NCHAN
              SB(J,I)=BG-FG
            END DO
          END IF
        END DO
        DO J=1,NCHAN
          IF(.NOT.LFITX(J))THEN
            DO I=1,NSCAN
              SB(J,I)=BG-FG
            END DO
          END IF
        END DO
        DO I=1,NSCAN
          DO J=1,NCHAN
            IF(.NOT.LMASK(J,I))THEN
              SB(J,I)=BG-FG
            END IF
          END DO
        END DO
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGGRAY(SB,NCMAX,NSMAX,1,NCHAN,1,NSCAN,FG,BG,TR)
        END DO
C
        WRITE(*,100)'Press <CR> to see Knots...'
        READ(*,*)
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM)) CALL PGSCI(5)
          DO K=1,NKNOTX
            CALL PGMOVE(FKNOTX(K),YMIN_IMA)
            CALL PGDRAW(FKNOTX(K),YMAX_IMA)
          END DO
          IF(LCOLOR(ITERM)) CALL PGSCI(4)
          DO K=1,NKNOTY
            CALL PGMOVE(XMIN_IMA,FKNOTY(K))
            CALL PGDRAW(XMAX_IMA,FKNOTY(K))
          END DO
          IF(LCOLOR(ITERM)) CALL PGSCI(2)
          DO I=1,NKNOTY
            DO J=1,NKNOTX
              CALL PGPOINT(1,FKNOTX(J),FKNOTY(I),17)
            END DO
          END DO
          IF(LCOLOR(ITERM)) CALL PGSCI(3)
          CALL PGSCH(0.6)
          DO I=1,NKNOTY
            DO J=1,NKNOTX
              WRITE(CDUMMY,'(I10,A1,I10)')(I-1)*NKNOTX+J,'/',
     +         NKNOTX*NKNOTY
              CALL RMBLANK(CDUMMY,CDUMMY,L)
              CALL PGPTEXT(FKNOTX(J)+REAL(NCHAN)/150,FKNOTY(I),
     +         0.,0.,CDUMMY(1:L))
            END DO
          END DO
          CALL PGSCH(1.0)
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
        END DO
C------------------------------------------------------------------------------
C definimos parametros para ajustar superficie polinomia inicial de grado bajo
        CALL PIDEPOLY(NBINX,NBINY,BINMODE,BINSIGMA,GX,GY,.FALSE.)
C------------------------------------------------------------------------------
C mostramos el efecto del peso aplicado para unir los polinomios
        CPOK='n'
        DO WHILE(CPOK.EQ.'n')
          WRITE(*,100)'* Enter power for the weighting around knots:'
          WRITE(CDUMMY,*)PESOPOWER
          PESOPOWER=READILIM(CDUMMY,1,100)
C
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGVPORT(.05,.95,.05,.40)
            CALL PGWINDOW(0.,20.,0.,1.1)
            CALL PGBOX('BC',0.0,0,'BCTNS',0.0,0)
            IF(LCOLOR(ITERM))CALL PGSCI(2)
            CALL PGMOVE(0.,1.)
            CALL PGDRAW(3.,1.)
            NPLOT=200
            DO I=1,NPLOT
              X0=3.+REAL(I-1)*7./REAL(NPLOT-1)
              Y0=(10.-X0)/7.
              Y0=Y0**PESOPOWER
              IF(I.EQ.1)THEN
                CALL PGMOVE(X0,Y0)
              ELSE
                CALL PGDRAW(X0,Y0)
              END IF
            END DO
            DO I=1,NPLOT
              X0=10.+REAL(I-1)*5./REAL(NPLOT-1)
              Y0=(15.-X0)/5.
              Y0=Y0**PESOPOWER
              IF(I.EQ.1)THEN
                CALL PGMOVE(X0,Y0)
              ELSE
                CALL PGDRAW(X0,Y0)
              END IF
            END DO
            CALL PGMOVE(15.,1.)
            CALL PGDRAW(20.,1.)
            IF(LCOLOR(ITERM))CALL PGSCI(3)
            CALL PGMOVE(3.,0.)
            CALL PGDRAW(3.,1.1)
            CALL PGMOVE(10.0,0.)
            CALL PGDRAW(10.0,1.1)
            CALL PGMOVE(15.0,0.)
            CALL PGDRAW(15.0,1.1)
            IF(LCOLOR(ITERM))CALL PGSCI(1)
            CALL RMBLANK(CDUMMY,CDUMMY,L)
            CALL PGSCH(0.7)
            CALL PGMTEXT('T',1.0,0.5,0.5,'Power: '//CDUMMY(1:L))
            CALL PGSCH(1.0)
          END DO
C
          WRITE(*,100)'Is the power ok (y/n) '
          CPOK(1:1)=READC('y','yn')
        END DO
C------------------------------------------------------------------------------
C bucle: para cada knot, ajustamos una superficie polinomica de grado bajo
C usando la region delimitada por los knots que lo rodean (solo se utiliza los
C pixels comprendidos entre el Knot problema y las fracciones FX,FY de la 
C distancia a los knots vecinos).
        FX=-1.0
        FY=FX
        DO WHILE((FX.LE.0.0).OR.(FX.GT.1.0).OR.
     +   (FY.LE.0.0).OR.(FY.GT.1.0))
          WRITE(*,101)'* Enter fraction of the distance to '//
     +     'neighbouring knots:'
          WRITE(*,100)'Fraction in X-axis '
          FX=READF('0.8')
          WRITE(*,100)'Fraction in Y-axis '
          FY=READF('0.8')
          IF((FX.LE.0.0).OR.(FX.GT.1.0).OR.
     +     (FY.LE.0.0).OR.(FY.GT.1.0))THEN
            WRITE(*,101)'ERROR: fractions must be real numbers'//
     +       ' in the range (0,1]. Try again.'
          END IF
        END DO
C ajustamos las superficies polinomicas
        WRITE(*,100)'No. of iterations to reject points (0=NONE) '
        NITER=READI('1')
        IF(NITER.GT.0)THEN
          WRITE(*,100)'Times sigma to reject points '
          SIGMAREJ=READF('3.0')
        ELSE
          SIGMAREJ=0.0
        END IF
        CNEXT='@'
        DO I=1,NKNOTY
          Y(I)=FKNOTY(I)
          DO J=1,NKNOTX
            X(J)=FKNOTX(J)
            WRITE(CDUMMY,'(I10,A1,I10)')(I-1)*NKNOTX+J,'/',NKNOTX*NKNOTY
            CALL RMBLANK(CDUMMY,CDUMMY,L)
            GLABEL=CDUMMY(1:L)
            WRITE(*,101)'>>> Next Knot is #'//CDUMMY(1:L)
            SBFIT(J,I)=POLYSURF(J,I,GX,GY,GLABEL)
            IF(CNEXT.EQ.'n') CNEXT='@'
            IF((CNEXT.EQ.'j').AND.((I-1)*NKNOTX+J.EQ.NKJUMP)) CNEXT='@'
            DO WHILE((CNEXT.NE.'n').AND.(CNEXT.NE.'g').AND.
     +       (CNEXT.NE.'j'))
              WRITE(*,101)'>>> Last Knot is #'//CDUMMY(1:L)
              IF((I-1)*NKNOTX+J.EQ.NKNOTX*NKNOTY)THEN
                WRITE(*,100)'[r]epeat, [n]ext, [g]o, e[x]it '//
     +           '(r/n/g/x) '
                CNEXT(1:1)=READC('n','rRnNgGxX')
              ELSE
                WRITE(*,100)'[r]epeat, [n]ext, [j]ump, [g]o, e[x]it '//
     +           '(r/n/j/g/x) '
                CNEXT(1:1)=READC('n','rRnNjJgGxX')
              END IF
              CALL CHLOWER(CNEXT)
              IF(CNEXT.EQ.'x') GOTO 90
              IF(CNEXT.EQ.'j')THEN
                WRITE(*,100)'Knot number'
                NKJUMP=READILIM('@',(I-1)*NKNOTX+J+1,NKNOTX*NKNOTY)
              END IF
              IF(CNEXT.EQ.'r')THEN
                NBINX_OLD=NBINX
                NBINY_OLD=NBINY
                BINMODE_OLD=BINMODE
                BINSIGMA_OLD=BINSIGMA
                GX_OLD=GX
                GY_OLD=GY
                CALL PIDEPOLY(NBINX,NBINY,BINMODE,BINSIGMA,GX,GY,
     +           .TRUE.)
                SBFIT(J,I)=POLYSURF(J,I,GX,GY,GLABEL)
                NBINX=NBINX_OLD
                NBINY=NBINY_OLD
                BINMODE=BINMODE_OLD
                BINSIGMA=BINSIGMA_OLD
                GX=GX_OLD
                GY=GY_OLD
              END IF
            END DO
          END DO
        END DO
C------------------------------------------------------------------------------
        WRITE(*,*)
        WRITE(*,101)'* Next frame will be created by using the '//
     +   'local polynomial fits:'
        WRITE(*,100)'Save fitted image -not binned-......(y/n)'
        CSAVE(1:1)=READC('@','yn')
C
        IF(CSAVE.EQ.'n') GOTO 80
C
        DO I=1,NSCAN
          DO J=1,NCHAN
            FINAL(J,I)=FINAL(J,I)/FINAL_PESO(J,I)
          END DO
        END DO
C
        WRITE(*,100)'Output file name'
        OUTFILE=OUTFILEX(30,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        IF(CSUB.EQ.'y')THEN
          DO I=1,NSCAN
            DO J=1,NCHAN
              FINAL(J,I)=FINAL(J,I)/FCTE2+FCTE1
            END DO
          END DO
        END IF
        DO I=1,NSCAN
          WRITE(30) (FINAL(J,I),J=1,NCHAN)
        END DO
        CLOSE(30)
C------------------------------------------------------------------------------
        WRITE(*,100)'Save residual image -not binned-....(y/n)'
        CSAVE(1:1)=READC('@','yn')
        IF(CSAVE.EQ.'y')THEN
          WRITE(*,100)'Output file name'
          OUTFILE=OUTFILEX(30,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
          IF(CSUB.EQ.'y')THEN
            DO I=1,NSCAN
              DO J=1,NCHAN
                SB(J,I)=S(J,I)/FCTE2+FCTE1
              END DO
            END DO
          ELSE
            DO I=1,NSCAN
              DO J=1,NCHAN
                SB(J,I)=S(J,I)
              END DO
            END DO
          END IF
          DO I=1,NSCAN
            DO J=1,NCHAN
              SB(J,I)=SB(J,I)-FINAL(J,I)
            END DO
            WRITE(30) (SB(J,I),J=1,NCHAN)
          END DO
          CLOSE(30)
        END IF
C------------------------------------------------------------------------------
C Calculamos la superficie mediante splines
80      WRITE(*,101)'* Next frame will be created by using '//
     +   'bicubic splines:'
        WRITE(*,100)'Save fitted image -not binned-......(y/n)'
        CSAVE(1:1)=READC('@','yn')
        IF(CSAVE.EQ.'n') GOTO 90
C
!?      CALL BICUBSPL(X,Y,SBFIT,NKNOTX,NKNOTY,NCMAX,NSMAX,AA,BB,CC)
        CALL BICUBSPL(Y,SBFIT,NKNOTX,NKNOTY,NCMAX,NSMAX,AA,BB,CC)
        DO I=1,NSCAN
          Y0=REAL(I)
          DO J=1,NCHAN
            X0=REAL(J)
            CALL BICUBSPLX(X,Y,SBFIT,NKNOTX,NKNOTY,NCMAX,NSMAX,
     +                     AA,BB,CC,X0,Y0,Z0)
            SB(J,I)=Z0
          END DO
        END DO
        WRITE(*,*) !no borrar
        WRITE(*,100)'Output file name'
        OUTFILE=OUTFILEX(30,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        IF(CSUB.EQ.'y')THEN
          DO I=1,NSCAN
            DO J=1,NCHAN
              SB(J,I)=SB(J,I)/FCTE2+FCTE1
            END DO
          END DO
        END IF
        DO I=1,NSCAN
          WRITE(30) (SB(J,I),J=1,NCHAN)
        END DO
        CLOSE(30)
C------------------------------------------------------------------------------
        WRITE(*,100)'Save residual image -not binned-....(y/n)'
        CSAVE(1:1)=READC('@','yn')
        IF(CSAVE.EQ.'y')THEN
          WRITE(*,100)'Output file name'
          OUTFILE=OUTFILEX(30,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
          IF(CSUB.EQ.'y')THEN
            DO I=1,NSCAN
              DO J=1,NCHAN
                SBFIT(J,I)=S(J,I)/FCTE2+FCTE1
              END DO
            END DO
          ELSE
            DO I=1,NSCAN
              DO J=1,NCHAN
                SBFIT(J,I)=S(J,I)
              END DO
            END DO
          END IF
          DO I=1,NSCAN
            DO J=1,NCHAN
              SBFIT(J,I)=SBFIT(J,I)-SB(J,I)
            END DO
            WRITE(30) (SBFIT(J,I),J=1,NCHAN)
          END DO
          CLOSE(30)
        END IF
C------------------------------------------------------------------------------
90        WRITE(*,100)'Repeat whole fit (y/n) '
        CRESTART(1:1)=READC('@','yn')
        IF(CRESTART.EQ.'y') GOTO 30
C------------------------------------------------------------------------------
        CALL PGEND
100     FORMAT(A,$)
101     FORMAT(A)
        STOP
C
        END
C
C******************************************************************************
C
C Devuelve los Knots seleccionados sobre un corte promediado (espacial/ldo).
C
        SUBROUTINE PIDEKNOT(N,LFIT,CORTE,NKNOT,KNOT,YRMSTOL)
        IMPLICIT NONE
C
        INCLUDE 'redlib.inc'
        INCLUDE 'futils.inc'
        INTEGER READI
        REAL READF
C
        INTEGER NMAXKNOT
        PARAMETER(NMAXKNOT=30)
C
        INTEGER N
        LOGICAL LFIT(N)
        REAL CORTE(N)
        INTEGER NKNOT
        INTEGER KNOT(NMAXKNOT)
        REAL YRMSTOL
C
        INTEGER I,K,L
        INTEGER N1,N2,NFIT
        INTEGER NP,ND0
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        REAL XP(NCMAX),YP(NCMAX)          !dimensionado al mayor de NCMAX,NSMAX
        REAL FKNOT(NMAXKNOT)
        REAL XC,YC
        REAL YMIN,YMAX
        CHARACTER*1 CKMOD,CH
        CHARACTER*50 CDUMMY
        LOGICAL LCOLOR(MAX_ID_RED)
C
        COMMON/BLKYLIMITS/YMIN,YMAX
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        DO I=1,N
          LFIT(I)=.FALSE.
        END DO
        WRITE(CDUMMY,'(A,I10)')'1,',N
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,*)
        WRITE(*,101)'* Fitting splines:'
        WRITE(*,101)'Valid region is: '//CDUMMY(1:L)
10        WRITE(*,100)'1st and last position to be used (0,0=EXIT) '
        CALL READ2I('0,0',N1,N2)
        IF((N1.EQ.0).AND.(N2.EQ.0)) GOTO 20
        IF((N1.LT.1).OR.(N2.GT.N).OR.(N1.GT.N2))THEN
          WRITE(*,101)'ERROR: numbers out of range. Try again.'
          GOTO 10
        END IF
        DO I=N1,N2
          LFIT(I)=.TRUE.
        END DO
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM)) CALL PGSCI(4)
          CALL PGMOVE(REAL(N1),YMIN)
          CALL PGDRAW(REAL(N1),YMAX)
          CALL PGMOVE(REAL(N2),YMIN)
          CALL PGDRAW(REAL(N2),YMAX)
          CALL PGRECT(REAL(N1),REAL(N2),YMIN,YMIN+(YMAX-YMIN)*.01)
          CALL PGRECT(REAL(N1),REAL(N2),YMAX,YMAX-(YMAX-YMIN)*.01)
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
        END DO
        GOTO 10
C
20        NFIT=0
        DO I=1,N
          IF(LFIT(I)) NFIT=NFIT+1
        END DO
        IF(NFIT.EQ.0)THEN
          WRITE(*,101)'ERROR: no. of points selected = 0. Try again.'
          GOTO 10
        END IF
C
        NP=0
        DO I=1,N
          IF(LFIT(I))THEN
            NP=NP+1
            XP(NP)=REAL(I)
            YP(NP)=CORTE(I)
          END IF
        END DO
C calculamos el punto minimo y maximo
        DO I=N,1,-1
          IF(LFIT(I)) N1=I
        END DO
        DO I=1,N
          IF(LFIT(I)) N2=I
        END DO
C
        NKNOT=2
        KNOT(1)=N1
        KNOT(2)=N2
        WRITE(*,110)'>>> Current number of Knots: ',NKNOT
        DO K=1,NKNOT
          WRITE(*,'(A,I2,A,I6)')'Knot #',K,'  located at: ',KNOT(K)
        END DO
C
        CKMOD=' '
C
        WRITE(CDUMMY,'(I10,A1,I10)')N1,',',N2
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
46          WRITE(*,'(A,I2,A,$)')'New Knot #',ND0+1,
     +     '   Press mouse (q=X=EXIT)...'
          CALL PGBAND(6,0,0.,0.,XC,YC,CH)
          WRITE(*,110)'    position: ',NINT(XC)
          IF((CH.EQ.'q').OR.(CH.EQ.'Q').OR.(CH.EQ.'X').OR.
     +     (CH.EQ.'x'))THEN
            GOTO 47
          END IF
          KNOT(NKNOT+ND0+1)=NINT(XC)
          IF( (KNOT(NKNOT+ND0+1).LE.N1).OR.
     +     (KNOT(NKNOT+ND0+1).GE.N2) )THEN
            WRITE(*,101)'ERROR: knot position out of range. Try again.'
            GOTO 46
          END IF
          ND0=ND0+1
          IF(NKNOT+ND0.EQ.NMAXKNOT)THEN
            WRITE(*,110)'WARNING: Maximum number of Knots is:',NMAXKNOT
            WRITE(*,101)'No more Knots allowed.'
            GOTO 47
          END IF
          GOTO 46
47          CONTINUE
        ELSE
48          WRITE(*,100)'No. of new Knots between '//CDUMMY(1:L)//
     +       ' (excluding '//CDUMMY(1:L)//')'
          ND0=READI('@')
          IF(NKNOT+ND0.GT.NMAXKNOT)THEN
            WRITE(*,110)'ERROR: Maximum number of knots is:',NMAXKNOT
            WRITE(*,101)'       Try again.'
            GOTO 48
          END IF
          IF(CKMOD.EQ.'a')THEN
            DO K=NKNOT+1,NKNOT+ND0
              KNOT(K)=NINT(REAL(N1)+REAL(K-NKNOT)*
     +         REAL(N2-N1)/REAL(ND0+1))
            END DO
          ELSE
            DO K=NKNOT+1,NKNOT+ND0
49              WRITE(*,'(A,I2.2,$)')'New Knot #',K-NKNOT
              KNOT(K)=READI('@')
              IF((KNOT(K).LE.N1).OR.(KNOT(K).GE.N2))THEN
                WRITE(*,101)'ERROR: knot position out of range. '//
     +           'Try again.'
                GOTO 49
              END IF
              DO L=1,K-1
                IF(KNOT(L).EQ.KNOT(K))THEN
                  WRITE(*,101)'ERROR: this Knot already exist. '//
     +             'Try again.'
                  GOTO 49
                END IF
              END DO
            END DO
          END IF
        END IF
C
        NKNOT=NKNOT+ND0
        WRITE(*,110)'Total number of Knots: ',NKNOT
        DO K=1,NKNOT
          FKNOT(K)=REAL(KNOT(K))
        END DO
        CALL ORDENA1F(NKNOT,FKNOT)
        DO K=1,NKNOT
          KNOT(K)=NINT(FKNOT(K))
        END DO
C
        IF(CKMOD.EQ.'a')THEN
          DO K=1,NKNOT
            WRITE(CDUMMY,*)KNOT(K)
            CALL RMBLANK(CDUMMY,CDUMMY,L)
            WRITE(*,'(A,I2,A)')'Knot #',K,' located at '//CDUMMY(1:L)
          END DO
        END IF
C
        WRITE(CDUMMY,*) YRMSTOL
        WRITE(*,100)'YRMSTOL for DOWNHILL '
        YRMSTOL=READF(CDUMMY)
C
100     FORMAT(A,$)
101     FORMAT(A)
110     FORMAT(A,I6)
        END
C
C******************************************************************************
C Ajusta una superficie polinomica de bajo grado entre los canales (JJ1,JJ2) y
C los scans (II1,II2), usando solo las regiones que pueden ser ajustadas
C (definidas por LFITX,LFITY). Los grados de los polinomios en X e Y son,
C respectivamente, GX y GY. Se realiza un binning en la region a ajustar
C usando los valores NBINX,NBINY. Este binning se lleva a cabo sustituyendo
C el nuevo pixel por la media o la mediana de los valores iniciales
C (dependiendo de la variable BINMODE). Como retorno la funcion devuelve el 
C valor de la superficie ajustada en el Knot numero (JKNOT,IKNOT).
C
        REAL FUNCTION POLYSURF(JKNOT,IKNOT,GGX,GGY,GLABEL)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INTEGER JKNOT,IKNOT
        INTEGER GGX,GGY
        CHARACTER*(*) GLABEL
C
        INTEGER MAXNCOEFF
        PARAMETER(MAXNCOEFF=100)
        INTEGER NMAXBINX,NMAXBINY
        PARAMETER(NMAXBINX=101,NMAXBINY=101)
        INTEGER NMAXKNOT
        PARAMETER(NMAXKNOT=30)
C
        INTEGER K,I,J,L,M
        INTEGER II,JJ,NN
        INTEGER III,JJJ
        INTEGER JJ1,JJ2,II1,II2
        INTEGER JJJ1,JJJ2,III1,III2
        INTEGER NX,NY
        INTEGER NBINX,NBINY
        INTEGER GX,GY,NCOEFF,BINMODE
        INTEGER NCHANB,NSCANB,NC0,NC1,NC2,NS0,NS1,NS2
        INTEGER P,Q,ORDER(MAXNCOEFF),IOK,IPAR
        INTEGER ITER,NITER,PESOPOWER
        INTEGER KNOTX(NMAXKNOT),KNOTY(NMAXKNOT),NKNOTX,NKNOTY
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        REAL FX,FY,FJKNOT,FIKNOT
        REAL FKNOTX(NMAXKNOT),FKNOTY(NMAXKNOT)
        REAL S(NCMAX,NSMAX),SB(NCMAX,NSMAX)
        REAL SS(NCMAX,NSMAX),SBIN(NCMAX,NSMAX)
        REAL FINAL(NCMAX,NSMAX),FINAL_PESO(NCMAX,NSMAX)
        REAL X(NCMAX),Y(NSMAX)
        REAL XX(NCMAX),SPX(NCMAX),YY(NSMAX),SPY(NSMAX)
        REAL PIXEL(NMAXBINX*NMAXBINY)
        REAL BINSIGMA,FMEAN1,FMEAN2,FMEDIAN,FMEDIAN1
        REAL A(MAXNCOEFF,MAXNCOEFF),B(MAXNCOEFF)
        REAL SCALEROW(MAXNCOEFF),XSOL(MAXNCOEFF)
        REAL TR(6),BG,FG,FMEAN,FSIGMA
        REAL XMIN,XMAX,YMIN,YMAX,DX,DY,YMIN1,YMIN2,YMAX1,YMAX2
        REAL SIGMAREJ,FSIGREJ,FSIGREJ0,FSIGREJMAX
        REAL CCX1,CCX2,CCY1,CCY2,FFACTOR
        REAL FPESO,FINAL0,PESO0
        CHARACTER*1 CNEXT
        CHARACTER*255 LOCALGLABEL
        LOGICAL LFITX(NCMAX),LFITY(NSMAX)
        LOGICAL LMASK(NCMAX,NSMAX)
        LOGICAL LFITBINX(NCMAX),LFITBINY(NSMAX)
        LOGICAL IFPIXEL(NCMAX,NSMAX)
        LOGICAL LPLOT
        LOGICAL LCOLOR(MAX_ID_RED)
C
        COMMON/BLKPIDE1/NSCAN,NCHAN
        COMMON/BLKPOLYNOM1/S,SB
        COMMON/BLKPOLYNOM2/NBINX,NBINY,BINMODE
        COMMON/BLKPOLYNOM3/LFITX,LFITY
        COMMON/BLKPOLYNOM3BIS/LMASK
        COMMON/BLKPOLYNOM4/TR,BINSIGMA
        COMMON/BLKPOLYNOM5/SIGMAREJ,FX,FY,FKNOTX,FKNOTY
        COMMON/BLKPOLYNOM6/NITER,PESOPOWER
        COMMON/BLKPOLYNOM7/FINAL,FINAL_PESO
        COMMON/BLKPOLYNOM8/KNOTX,KNOTY,NKNOTX,NKNOTY
        COMMON/BLKPOLYNOM9/CNEXT
        COMMON/BLKSS1/NX,NY
        COMMON/BLKSS2/X,Y,SS
        COMMON/BLKSS4/IFPIXEL
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        LPLOT=((CNEXT.NE.'j').AND.(CNEXT.NE.'g'))
        ITER=0
        FJKNOT=REAL(KNOTX(JKNOT))
        FIKNOT=REAL(KNOTY(IKNOT))
        GX=GGX                                   !duplicamos por si modificamos
        GY=GGY                                   !duplicamos por si modificamos
        NCOEFF=(GX+1)*(GY+1)
        IF(NCOEFF.GT.MAXNCOEFF)
     +   STOP 'FATAL ERROR in POLYSURF: NCOEFF too large.'
C------------------------------------------------------------------------------
C limites para realizar el ajuste polinomico
        II1=IKNOT-1
        II2=IKNOT+1
        IF(II1.LT.1)THEN
          II1=1
        ELSE
          II1=KNOTY(II1+1)-NINT((FKNOTY(II1+1)-FKNOTY(II1))*FY)
        END IF
        IF(II2.GT.NKNOTY)THEN
          II2=NSCAN
        ELSE
          II2=KNOTY(II2-1)+NINT((FKNOTY(II2)-FKNOTY(II2-1))*FY)
        END IF
C
        JJ1=JKNOT-1
        JJ2=JKNOT+1
        IF(JJ1.LT.1)THEN
          JJ1=1
        ELSE
          JJ1=KNOTX(JJ1+1)-NINT((FKNOTX(JJ1+1)-FKNOTX(JJ1))*FX)
        END IF
        IF(JJ2.GT.NKNOTX)THEN
          JJ2=NCHAN
        ELSE
          JJ2=KNOTX(JJ2-1)+NINT((FKNOTX(JJ2)-FKNOTX(JJ2-1))*FX)
        END IF
C------------------------------------------------------------------------------
C nuevo taman~o de la region debido al binning
        NCHANB=JJ2-JJ1+1
        NX=NCHANB/NBINX
        IF(MOD(NCHANB,NBINX).NE.0) NX=NX+1
        NSCANB=II2-II1+1
        NY=NSCANB/NBINY
        IF(MOD(NSCANB,NBINY).NE.0) NY=NY+1
        WRITE(*,100)'knot position (channel,scan): '
        WRITE(*,*)KNOTX(JKNOT),KNOTY(IKNOT)
        WRITE(*,100)'channels: '
        WRITE(*,*)jj1,jj2
        WRITE(*,100)'scans   : '
        WRITE(*,*)ii1,ii2
        WRITE(*,100)'initial dimension with binning (channels,scans): '
        WRITE(*,*)nx,ny
C------------------------------------------------------------------------------
        IF(.NOT.LPLOT) GOTO 7
C dibujamos la region de la imagen original que vamos a utilizar
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGVPORT(0.00,1.00,0.04,0.46)
          CALL PGSCI(0)
          CALL PGWINDOW(0.,1.,0.,1.)
          CALL PGRECT(0.,1.,0.,1.)
          CALL PGSCI(1)
          CALL PGVPORT(0.00,0.50,0.00,0.05)
          CALL PGSCI(0)
          CALL PGWINDOW(0.,1.,0.,1.)
          CALL PGRECT(0.,1.,0.,1.)
          CALL PGSCI(1)
          CALL PGVPORT(0.05,0.30,0.32,0.45)
          CALL PGWINDOW(REAL(JJ1)-.5,REAL(JJ2)+.5,
     +     REAL(II1)-.5,REAL(II2)+.5)
          CALL PGBOX('BC',0.0,0,'BC',0.0,0)
          CALL PGSCH(0.7)
          CALL PGMTEXT('L',1.0,0.5,0.5,'original')
          CALL PGSCH(1.0)
        END DO
C------------------------------------------------------------------------------
C calculamos nuevos cortes ajustados a la region a utilizar
        FMEAN=0.
        NN=0
        DO I=II1,II2
          IF(LFITY(I))THEN
            DO J=JJ1,JJ2
              IF(LFITX(J))THEN
                IF(LMASK(J,I))THEN
                  NN=NN+1
                  FMEAN=FMEAN+SB(J,I)
                END IF
              END IF
            END DO
          END IF
        END DO
        FMEAN=FMEAN/REAL(NN)
        FSIGMA=0.
        IF(NN.GT.1)THEN
          DO I=II1,II2
            IF(LFITY(I))THEN
              DO J=JJ1,JJ2
                IF(LFITX(J))THEN
                  IF(LMASK(J,I))THEN
                    FSIGMA=FSIGMA+(SB(J,I)-FMEAN)*(SB(J,I)-FMEAN)
                  END IF
                END IF
              END DO
            END IF
          END DO
        END IF
        FSIGMA=SQRT(FSIGMA/REAL(NN-1))
        BG=FMEAN-3.0*FSIGMA
        FG=FMEAN+3.0*FSIGMA
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGGRAY(SB,NCMAX,NSMAX,JJ1,JJ2,II1,II2,FG,BG,TR)
        END DO
C------------------------------------------------------------------------------
C determinamos que nuevos pixels (con binning) contienen pixels originales que
C pueden ser utilizados para ajustar la superficie polinomica y, de paso,
C calculamos tambien las coordenadas promedio de los nuevos pixels con binning
C (si estamos usando máscara, no hay binning)
7       DO JJ=1,NX
          NC1=JJ1+(JJ-1)*NBINX
          NC2=NC1+NBINX-1
          IF(NC2.GT.JJ2) NC2=JJ2
          LFITBINX(JJ)=.FALSE.
          NN=0
          X(JJ)=0.
          DO NC0=NC1,NC2
            IF(LFITX(NC0))THEN
              NN=NN+1
              X(JJ)=X(JJ)+REAL(NC0)
            END IF
          END DO
          IF(NN.GT.0)THEN
            LFITBINX(JJ)=.TRUE.
            X(JJ)=X(JJ)/REAL(NN)
          END IF
        END DO
C
        DO II=1,NY
          NS1=II1+(II-1)*NBINY
          NS2=NS1+NBINY-1
          IF(NS2.GT.II2) NS2=II2
          LFITBINY(II)=.FALSE.
          NN=0
          Y(II)=0.
          DO NS0=NS1,NS2
            IF(LFITY(NS0))THEN
              NN=NN+1
              Y(II)=Y(II)+REAL(NS0)
            END IF
          END DO
          IF(NN.GT.0)THEN
            LFITBINY(II)=.TRUE.
            Y(II)=Y(II)/REAL(NN)
          END IF
        END DO
C------------------------------------------------------------------------------
C normalizamos el recorrido de las variables X e Y al intervalo [-1,1]
        CCX1=2/(REAL(JJ2)-REAL(JJ1))
        CCX2=(REAL(JJ1)+REAL(JJ2))/(REAL(JJ2)-REAL(JJ1))
        DO JJ=1,NX
          X(JJ)=CCX1*X(JJ)-CCX2
        END DO
        CCY1=2/(REAL(II2)-REAL(II1))
        CCY2=(REAL(II1)+REAL(II2))/(REAL(II2)-REAL(II1))
        DO II=1,NY
          Y(II)=CCY1*Y(II)-CCY2
        END DO
C------------------------------------------------------------------------------
C calculamos el valor medio/mediana en cada pixel con binning que puede 
C ser utilizado
        DO II=1,NY
          DO JJ=1,NX
            SBIN(JJ,II)=BG-FG
          END DO
        END DO
C
        DO II=1,NY
          IF(LFITBINY(II))THEN
            NS1=II1+(II-1)*NBINY
            NS2=NS1+NBINY-1
            IF(NS2.GT.II2) NS2=II2
            DO JJ=1,NX
              IF(LFITBINX(JJ))THEN
                NC1=JJ1+(JJ-1)*NBINX
                NC2=NC1+NBINX-1
                IF(NC2.GT.JJ2) NC2=JJ2
                K=0
                DO NS0=NS1,NS2
                  IF(LFITY(NS0))THEN
                    DO NC0=NC1,NC2
                      IF(LFITX(NC0))THEN
                        K=K+1
                        PIXEL(K)=S(NC0,NS0)
                      END IF
                    END DO
                  END IF
                END DO
                IF(BINMODE.EQ.1)THEN                 !mean
                  SBIN(JJ,II)=FMEAN1(K,PIXEL)
                ELSEIF(BINMODE.EQ.2)THEN             !mean rejecting with sigma
                  SBIN(JJ,II)=FMEAN2(K,PIXEL,BINSIGMA)
                ELSEIF(BINMODE.EQ.3)THEN             !median
                  FMEDIAN=FMEDIAN1(K,PIXEL)
                  SBIN(JJ,II)=FMEDIAN
                ELSE
                  WRITE(*,101)'FATAL ERROR in subroutine POLYSURF. '//
     +             'Invalid BINMODE value.'
                  STOP
                END IF
              END IF
            END DO
          END IF
        END DO
C------------------------------------------------------------------------------
C dibujamos la imagen a utilizar tras el binning
        IF(LPLOT)THEN
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGVPORT(0.05,0.30,0.19,0.32)
            CALL PGWINDOW(1.0-.5,REAL(NX)+.5,1.0-.5,REAL(NY)+.5)
            CALL PGBOX('BC',0.0,0,'BC',0.0,0)
            CALL PGSCH(0.7)
            CALL PGMTEXT('L',1.0,0.5,0.5,'original binned')
            CALL PGSCH(1.0)
            CALL PGGRAY(SBIN,NCMAX,NSMAX,1,NX,1,NY,FG,BG,TR)
          END DO
        END IF
C------------------------------------------------------------------------------
C eliminamos las columnas de pixels con binning que no van a ajustarse
        JJJ=0
        DO JJ=1,NX
          IF(LFITBINX(JJ))THEN
            JJJ=JJJ+1
            X(JJJ)=X(JJ)
            DO II=1,NY
              SS(JJJ,II)=SBIN(JJ,II)
            END DO
          END IF
        END DO
        NX=JJJ
C------------------------------------------------------------------------------
C eliminamos las filas de pixels con binning que no van a ajustarse
        III=0
        DO II=1,NY
          IF(LFITBINY(II))THEN
            III=III+1
            Y(III)=Y(II)
            DO JJ=1,NX
              SS(JJ,III)=SS(JJ,II)
            END DO
          END IF
        END DO
        NY=III
C------------------------------------------------------------------------------
C almacenamos una copia de la imagen que va a ser ajustada
        DO II=1,NY
          DO JJ=1,NX
            SBIN(JJ,II)=SS(JJ,II)
          END DO
        END DO
C------------------------------------------------------------------------------
C nueva dimension para el ajuste
        WRITE(*,100)'final dimension with binning (channels,scans)..: '
        WRITE(*,*)nx,ny
        DO WHILE(NX.LT.GX+1)
          WRITE(*,101)'WARNING: image size with binning is too small.'
          WRITE(*,'(A,I2)')'Last polynomial degree in X-axis: ',GX
          GX=GX-1
          WRITE(*,'(A,I2)')'New  polynomial degree in X-axis: ',GX
        END DO
        DO WHILE(NY.LT.GY+1)
          WRITE(*,101)'WARNING: image size with binning is too small.'
          WRITE(*,'(A,I2)')'Last polynomial degree in Y-axis: ',GY
          GY=GY-1
          WRITE(*,'(A,I2)')'New  polynomial degree in Y-axis: ',GY
        END DO
C------------------------------------------------------------------------------
C inicialmente todos los puntos disponibles seran usados en el ajuste
        DO M=1,NY
          DO L=1,NX
!           IFPIXEL(L,M)=.TRUE.
            !introducimos aquí el uso de la máscara. OJO: estamos suponiendo 
            !que no hay binning
            IFPIXEL(L,M)=(LMASK(L+JJ1-1,M+II1-1))
          END DO
        END DO
C------------------------------------------------------------------------------
C dibujamos la superficie antes de ajustarla
        IF(LPLOT)THEN
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(4)
            CALL PLOTSS(.FALSE.,.TRUE.,1,LCOLOR(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END DO
        END IF
10      CONTINUE
        IF(LPLOT)THEN
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGSCH(0.7)
            CALL PGLABEL(' ','original binned',' ')
            CALL PGSCH(1.0)
          END DO
        END IF
C------------------------------------------------------------------------------
C definimos un sistema de (GX+1)*(GY+1) ecuaciones con (GX+1)*(GY+1) incognitas
C II es el numero de ecuacion, y JJ el numero de incognita, por lo que
C el sistema de ecuaciones es A(II,JJ) * X = B(II)
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
                    IF(IFPIXEL(L,M))THEN !incluimos la máscara
                      FFACTOR=1.
                      IF(I+P.NE.0) FFACTOR=FFACTOR*(X(L)**(I+P))
                      IF(J+Q.NE.0) FFACTOR=FFACTOR*(Y(M)**(J+Q))
                      A(II,JJ)=A(II,JJ)+FFACTOR
CCC hemos cambiado la linea de abajo por las cuatro de arriba para
CCC evitar la indeterminacion 0.**0, que debe dar 1 para que funcione bien.
ccc                   A(II,JJ)=A(II,JJ)+(X(L)**(I+P))*(Y(M)**(J+Q))
                    END IF
                  END DO
                END DO
              END DO
            END DO
            B(II)=0.
            DO L=1,NX
              DO M=1,NY
                IF(IFPIXEL(L,M))THEN
                  FFACTOR=SS(L,M)
                  IF(P.NE.0) FFACTOR=FFACTOR*(X(L)**(P))
                  IF(Q.NE.0) FFACTOR=FFACTOR*(Y(M)**(Q))
                  B(II)=B(II)+FFACTOR
CCC hemos cambiado la linea de abajo por las cuatro de arriba para
CCC evitar la indeterminacion 0.**0, que debe dar 1 para que funcione bien.
ccc               B(II)=B(II)+SS(L,M)*(X(L)**(P))*(Y(M)**(Q))
                END IF
              END DO
            END DO
          END DO
        END DO
C------------------------------------------------------------------------------
C resolvemos el sistema de ecuaciones
        WRITE(*,'(A,I3,A,I3,A)')'---> Solving ',NCOEFF,
     +   ' equations with ',NCOEFF,' unknowns...'
        WRITE(*,100)'---> LU Descomposition...'
        CALL LUDCMP(A,NCOEFF,MAXNCOEFF,ORDER,SCALEROW,IOK,IPAR)
        WRITE(*,101)' OK!'
        WRITE(*,100)'---> Forward substitution and back substitution...'
        CALL LUSOLV(A,NCOEFF,MAXNCOEFF,ORDER,SCALEROW,B,XSOL)
        WRITE(*,101)' OK!'
C------------------------------------------------------------------------------
C calculamos la superficie ajustada
        IF(LPLOT)THEN
          DO M=1,NY
            DO L=1,NX
              SS(L,M)=0.
              K=0
              DO I=0,GX
                DO J=0,GY
                  K=K+1
                  FFACTOR=XSOL(K)
                  IF(I.NE.0) FFACTOR=FFACTOR*(X(L)**(I))
                  IF(J.NE.0) FFACTOR=FFACTOR*(Y(M)**(J))
                  SS(L,M)=SS(L,M)+FFACTOR
CCC hemos cambiado la linea de abajo por las cuatro de arriba para
CCC evitar la indeterminacion 0.**0, que debe dar 1 para que funcione bien.
ccc                  SS(L,M)=SS(L,M)+XSOL(K)*(X(L)**(I))*(Y(M)**(J))
                END DO
              END DO
            END DO
          END DO
C------------------------------------------------------------------------------
C dibujamos la superficie ajustada
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(5)
            CALL PLOTSS(.FALSE.,.FALSE.,2,LCOLOR(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
            CALL PGSCH(0.7)
            CALL PGLABEL(' ','fitted binned',' ')
            CALL PGSCH(1.0)
          END DO
C------------------------------------------------------------------------------
C dibujamos la imagen suavizada sin binning
          DO II=II1,II2
            DO JJ=JJ1,JJ2
              SS(JJ,II)=0.
              K=0
              DO I=0,GX
                DO J=0,GY
                  K=K+1
                  FFACTOR=XSOL(K)
                  IF(I.NE.0) FFACTOR=FFACTOR*((CCX1*REAL(JJ)-CCX2)**(I))
                  IF(J.NE.0) FFACTOR=FFACTOR*((CCY1*REAL(II)-CCY2)**(J))
                  SS(JJ,II)=SS(JJ,II)+FFACTOR
CCC hemos cambiado las tres lineas de abajo por las cuatro de arriba para
CCC evitar la indeterminacion 0.**0, que debe dar 1 para que funcione bien.
ccc                  SS(JJ,II)=SS(JJ,II)+XSOL(K)*
ccc     +             ((CCX1*REAL(JJ)-CCX2)**(I))*
ccc     +             ((CCY1*REAL(II)-CCY2)**(J))
                END DO
              END DO
            END DO
          END DO
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGVPORT(0.05,0.30,0.06,0.19)
            CALL PGWINDOW(REAL(JJ1)-0.5,REAL(JJ2)+0.5,
     +       REAL(II1)-0.5,REAL(II2)+0.5)
            CALL PGBOX('BC',0.0,0,'BC',0.0,0)
            CALL PGSCH(0.7)
            CALL PGMTEXT('L',1.0,0.5,0.5,'fitted')
            CALL PGSCH(1.0)
            IF(LCOLOR(ITERM)) CALL PGSCI(3)
            LOCALGLABEL(1:6)='Knot #'
            LOCALGLABEL(7:)=GLABEL
            CALL PGMTEXT('B',1.5,0.5,0.5,LOCALGLABEL)
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
            CALL PGGRAY(SS,NCMAX,NSMAX,JJ1,JJ2,II1,II2,FG,BG,TR)
          END DO
        END IF
C------------------------------------------------------------------------------
C calculamos el valor en la imagen FINAL usando los pesos adecuados en funcion
C de la distancia a los Knots que van a influir en cada punto
        III1=IKNOT-1
        III2=IKNOT+1
        IF(III1.LT.1)THEN
          III1=1
        ELSE
          III1=KNOTY(III1)
        END IF
        IF(III2.GT.NKNOTY)THEN
          III2=NSCAN
        ELSE
          III2=KNOTY(III2)
        END IF
C
        JJJ1=JKNOT-1
        JJJ2=JKNOT+1
        IF(JJJ1.LT.1)THEN
          JJJ1=1
        ELSE
          JJJ1=KNOTX(JJJ1)
        END IF
        IF(JJJ2.GT.NKNOTX)THEN
          JJJ2=NCHAN
        ELSE
          JJJ2=KNOTX(JJJ2)
        END IF
C
        DO II=III1,III2
          DO JJ=JJJ1,JJJ2
            FINAL0=0.
            K=0
            DO I=0,GX
              DO J=0,GY
                K=K+1
                FFACTOR=XSOL(K)
                IF(I.NE.0) FFACTOR=FFACTOR*((CCX1*REAL(JJ)-CCX2)**(I))
                IF(J.NE.0) FFACTOR=FFACTOR*((CCY1*REAL(II)-CCY2)**(J))
                FINAL0=FINAL0+FFACTOR
CCC hemos cambiado las tres lineas de abajo por las cuatro de arriba para
CCC evitar la indeterminacion 0.**0, que debe dar 1 para que funcione bien.
ccc                FINAL0=FINAL0+XSOL(K)*
ccc     +           ((CCX1*REAL(JJ)-CCX2)**(I))*
ccc     +           ((CCY1*REAL(II)-CCY2)**(J))
              END DO
            END DO
            PESO0=FPESO(JKNOT,IKNOT,JJ,II)      !usaremos los pesos al cuadrado
            FINAL(JJ,II)=FINAL(JJ,II)+FINAL0*(PESO0**PESOPOWER)
            FINAL_PESO(JJ,II)=FINAL_PESO(JJ,II)+(PESO0**PESOPOWER)
          END DO
        END DO
C------------------------------------------------------------------------------
C dibujamos corte promedio en X e Y de la imagen original y del ajuste (sin 
C binning), usando solo las regiones permitidas en el ajuste
        IF(.NOT.LPLOT) GOTO 20
        NN=0
        DO II=II1,II2
          IF(LFITY(II)) NN=NN+1
        END DO
        DO JJ=JJ1,JJ2
          J=JJ-JJ1+1
          SPX(J)=0.
          DO II=II1,II2
            IF(LFITY(II)) SPX(J)=SPX(J)+S(JJ,II)
          END DO
          SPX(J)=SPX(J)/REAL(NN)
          XX(J)=REAL(JJ)
        END DO
        CALL FINDMM(JJ2-JJ1+1,SPX,YMIN1,YMAX1)
C
        NN=0
        DO JJ=JJ1,JJ2
          IF(LFITX(JJ)) NN=NN+1
        END DO
        DO II=II1,II2
          I=II-II1+1
          SPY(I)=0.
          DO JJ=JJ1,JJ2
            IF(LFITX(JJ)) SPY(I)=SPY(I)+S(JJ,II)
          END DO
          SPY(I)=SPY(I)/REAL(NN)
          YY(I)=REAL(II)
        END DO
        CALL FINDMM(II2-II1+1,SPY,YMIN2,YMAX2)
C
        IF(YMIN1.LT.YMIN2)THEN
          YMIN=YMIN1
        ELSE
          YMIN=YMIN2
        END IF
        IF(YMAX1.GT.YMAX2)THEN
          YMAX=YMAX1
        ELSE
          YMAX=YMAX2
        END IF
        DY=YMAX-YMIN
        YMIN=YMIN-DY/30.
        YMAX=YMAX+DY/30.
C
        XMIN=REAL(JJ1)
        XMAX=REAL(JJ2)
        DX=XMAX-XMIN
        XMIN=XMIN-DX/30.
        XMAX=XMAX+DX/30.
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGVPORT(0.70,0.95,0.29,0.44)
          CALL PGWINDOW(XMIN,XMAX,YMIN,YMAX)
          CALL PGSCH(0.7)
          CALL PGMTEXT('T',0.5,0.5,0.5,'channel')
          CALL PGBOX('BCNTS',0.0,0,'BCNTS',0.0,0)
          CALL PGSCH(1.0)
          IF(LCOLOR(ITERM)) CALL PGSCI(4)
          CALL PGBIN(JJ2-JJ1+1,XX,SPX,.TRUE.)
          CALL PGBBUF
          IF(LCOLOR(ITERM)) CALL PGSCI(2)
          DO JJ=JJ1,JJ2
            J=JJ-JJ1+1
            IF(.NOT.LFITX(JJ)) CALL PGPOINT(1,XX(J),SPX(J),5)
          END DO
          CALL PGEBUF
          DO JJ=JJ1,JJ2
            J=JJ-JJ1+1
            SPX(J)=0.
            DO II=II1,II2
              SPX(J)=SPX(J)+SS(JJ,II)
            END DO
            SPX(J)=SPX(J)/REAL(II2-II1+1)
            XX(J)=REAL(JJ)
          END DO
          IF(LCOLOR(ITERM)) CALL PGSCI(5)
          CALL PGBIN(JJ2-JJ1+1,XX,SPX,.TRUE.)
          IF(LCOLOR(ITERM)) CALL PGSCI(3)
          CALL PGSLS(2)
          CALL PGMOVE(FJKNOT,YMIN)
          CALL PGDRAW(FJKNOT,YMAX)
          CALL PGSLS(1)
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
        END DO
C
        XMIN=REAL(II1)
        XMAX=REAL(II2)
        DX=XMAX-XMIN
        XMIN=XMIN-DX/30.
        XMAX=XMAX+DX/30.
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGVPORT(0.70,0.95,0.07,0.22)
          CALL PGWINDOW(XMIN,XMAX,YMIN,YMAX)
          CALL PGSCH(0.7)
          CALL PGMTEXT('T',0.5,0.5,0.5,'scan')
          CALL PGBOX('BCNTS',0.0,0,'BCNTS',0.0,0)
          CALL PGSCH(1.0)
          IF(LCOLOR(ITERM)) CALL PGSCI(4)
          CALL PGBIN(II2-II1+1,YY,SPY,.TRUE.)
          CALL PGBBUF
          IF(LCOLOR(ITERM)) CALL PGSCI(2)
          DO II=II1,II2
            I=II-II1+1
            IF(.NOT.LFITY(II)) CALL PGPOINT(1,YY(I),SPY(I),5)
          END DO
          CALL PGEBUF
          DO II=II1,II2
            I=II-II1+1
            SPY(I)=0.
            DO JJ=JJ1,JJ2
              SPY(I)=SPY(I)+SS(JJ,II)
            END DO
            SPY(I)=SPY(I)/REAL(JJ2-JJ1+1)
            YY(I)=REAL(II)
          END DO
          IF(LCOLOR(ITERM)) CALL PGSCI(5)
          CALL PGBIN(II2-II1+1,YY,SPY,.TRUE.)
          IF(LCOLOR(ITERM)) CALL PGSCI(3)
          CALL PGSLS(2)
          CALL PGMOVE(FIKNOT,YMIN)
          CALL PGDRAW(FIKNOT,YMAX)
          CALL PGSLS(1)
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
        END DO
C------------------------------------------------------------------------------
C calculamos el valor de la superficie en el punto buscado
20      POLYSURF=0.
        K=0
        DO I=0,GX
          DO J=0,GY
            K=K+1
            FFACTOR=XSOL(K)
            IF(I.NE.0) FFACTOR=FFACTOR*((CCX1*FJKNOT-CCX2)**(I))
            IF(J.NE.0) FFACTOR=FFACTOR*((CCY1*FIKNOT-CCY2)**(J))
            POLYSURF=POLYSURF+FFACTOR
CCC hemos cambiado las tres lineas de abajo por las cuatro de arriba para
CCC evitar la indeterminacion 0.**0, que debe dar 1 para que funcione bien.
ccc            POLYSURF=POLYSURF+XSOL(K)*
ccc     +       ((CCX1*FJKNOT-CCX2)**(I))*
ccc     +       ((CCY1*FIKNOT-CCY2)**(J))
          END DO
        END DO
C------------------------------------------------------------------------------
C iteramos NITER veces
        IF(ITER.EQ.NITER) RETURN
C calculamos, otra vez, la superficie ajustada (con binning)
C NOTA: la superficie con binning con la que hay que comparar en SBIN
        DO M=1,NY
          DO L=1,NX
            SS(L,M)=0.
            K=0
            DO I=0,GX
              DO J=0,GY
                K=K+1
                FFACTOR=XSOL(K)
                IF(I.NE.0) FFACTOR=FFACTOR*(X(L)**(I))
                IF(J.NE.0) FFACTOR=FFACTOR*(Y(M)**(J))
                SS(L,M)=SS(L,M)+FFACTOR
CCC hemos cambiado la linea de abajo por las cuatro de arriba para
CCC evitar la indeterminacion 0.**0, que debe dar 1 para que funcione bien.
ccc                SS(L,M)=SS(L,M)+XSOL(K)*(X(L)**(I))*(Y(M)**(J))
              END DO
            END DO
          END DO
        END DO
C calculamos la desviacion alrededor de la superficie
        FSIGREJ=0.
        FSIGREJMAX=0.
        DO M=1,NY
          DO L=1,NX
            FSIGREJ0=ABS(SS(L,M)-SBIN(L,M))
            FSIGREJ=FSIGREJ+FSIGREJ0*FSIGREJ0
            IF(FSIGREJ0.GT.FSIGREJMAX) FSIGREJMAX=FSIGREJ0
          END DO
        END DO
        FSIGREJ=SQRT(FSIGREJ/REAL(NY*NX-1))
C mostramos SIGMA maximo
        WRITE(*,100)'Maximum deviation in sigma units: '
        WRITE(*,*)FSIGREJMAX/FSIGREJ
C eliminamos del ajuste puntos alejados
        DO M=1,NY
          DO L=1,NX
            IFPIXEL(L,M)=(ABS(SS(L,M)-SBIN(L,M)).LE.FSIGREJ*SIGMAREJ)
          END DO
        END DO
        !incluimos los puntos enmascarados
        DO M=1,NY
          DO L=1,NX
            IF(.NOT.LMASK(L+JJ1-1,M+II1-1)) IFPIXEL(L,M)=.FALSE.
          END DO
        END DO
C comprobamos que el numero de puntos a ajustar es suficiente para el grado
C de los polinomios
        NN=0
        DO M=1,NY
          DO L=1,NX
            IF(IFPIXEL(L,M)) NN=NN+1
          END DO
        END DO
        IF((GX+1)*(GY+1).GT.NN)THEN
          WRITE(*,101)'WARNING: insufficient no. of points'
          RETURN
        END IF
C
        WRITE(*,110)'No. of points available: ',NX*NY
        WRITE(*,110)'No. of points rejected : ',NX*NY-NN
        IF(NX*NY-NN.EQ.0) RETURN
C preparamos "hueco" para nuevos dibujos
        IF(LPLOT)THEN
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGVPORT(0.32,1.00,0.04,0.46)
            CALL PGSCI(0)
            CALL PGWINDOW(0.,1.,0.,1.)
            CALL PGRECT(0.,1.,0.,1.)
            CALL PGSCI(1)
          END DO
        END IF
C iteramos
        DO II=1,NY
          DO JJ=1,NX
            SS(JJ,II)=SBIN(JJ,II)
          END DO
        END DO
        IF(LPLOT)THEN
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(4)
            CALL PLOTSS(.FALSE.,.FALSE.,1,LCOLOR(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END DO
        END IF
        ITER=ITER+1
        GOTO 10
C------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
110     FORMAT(A,I6)
        END
C
C******************************************************************************
C Calcula el peso para el pixel JJ,II cuando la superficie ajustada es la
C calculada alrededor del knot numero (JK,IK). Los pesos que calcula
C esta funcion son lineales en cada eje, valiendo 1 para el pixel sobre el knot C y 0 cuando llegamos a los knots vecinos. 
        REAL FUNCTION FPESO(JK,IK,JJ,II)
        IMPLICIT NONE
        INTEGER JK,IK
        INTEGER JJ,II
C
        INTEGER NMAXKNOT
        PARAMETER(NMAXKNOT=30)
C
        INTEGER KNOTX(NMAXKNOT),KNOTY(NMAXKNOT),NKNOTX,NKNOTY
        REAL PESO_X,PESO_Y
C
        COMMON/BLKPOLYNOM8/KNOTX,KNOTY,NKNOTX,NKNOTY
C------------------------------------------------------------------------------
        IF(JJ.LE.KNOTX(JK))THEN                      !eje X, hacia la izquierda
          IF(JK.EQ.1)THEN
            PESO_X=1.
          ELSE
            PESO_X=REAL(JJ-KNOTX(JK-1))/REAL(KNOTX(JK)-KNOTX(JK-1))
          END IF
        ELSE                                           !eje X, hacia la derecha
          IF(JK.EQ.NKNOTX)THEN
            PESO_X=1.
          ELSE
            PESO_X=REAL(KNOTX(JK+1)-JJ)/REAL(KNOTX(JK+1)-KNOTX(JK))
          END IF
        END IF
C
        IF(II.LE.KNOTY(IK))THEN                             !eje Y, hacia abajo
          IF(IK.EQ.1)THEN
            PESO_Y=1.
          ELSE
            PESO_Y=REAL(II-KNOTY(IK-1))/REAL(KNOTY(IK)-KNOTY(IK-1))
          END IF
        ELSE                                               !eje Y, hacia arriba
          IF(IK.EQ.NKNOTY)THEN
            PESO_Y=1.
          ELSE
            PESO_Y=REAL(KNOTY(IK+1)-II)/REAL(KNOTY(IK+1)-KNOTY(IK))
          END IF
        END IF
C
        FPESO=PESO_X*PESO_Y
C
        END
C
C******************************************************************************
C
C dibuja la superficie de puntos
        SUBROUTINE PLOTSS(OVERP,LLIMITS,NPLOT,LCOLOR)
        IMPLICIT NONE
        LOGICAL OVERP,LLIMITS
        INTEGER NPLOT
        LOGICAL LCOLOR
        INCLUDE 'redlib.inc'
C
        INTEGER I,J
        INTEGER NX,NY
        INTEGER NSMAX_LOCAL,NCMAX_LOCAL
C XP,YP hay que dimensionarlos al mayor de NCMAX,NSMAX
        REAL XP(NCMAX),YP(NCMAX)
        REAL ZMINSS,ZMAXSS,OFFSS
        REAL SS(NCMAX,NSMAX)
        REAL X(NCMAX),Y(NSMAX),XX,YY
        REAL FACTOR
        LOGICAL IFPIXEL(NCMAX,NSMAX)
C
        COMMON/BLKSS1/NX,NY
        COMMON/BLKSS2/X,Y,SS
        COMMON/BLKSS3/ZMINSS,ZMAXSS,OFFSS
        COMMON/BLKSS4/IFPIXEL
C
C-------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        NSMAX_LOCAL=NSMAX
        NCMAX_LOCAL=NCMAX
        IF(NSMAX_LOCAL.GT.NCMAX_LOCAL)STOP'FATAL ERROR: NSMAX.GT.NCMAX'
C FACTOR determina como aparece de ampliada la estructura en el grafico
C que simula una representacion tridimensional de la superficie
        FACTOR=4.
C
        IF(.NOT.OVERP)THEN
          IF(LLIMITS)THEN
            ZMINSS=SS(1,1)
            ZMAXSS=ZMINSS
            DO I=1,NY
              DO J=1,NX
                IF(SS(J,I).LT.ZMINSS) ZMINSS=SS(J,I)
                IF(SS(J,I).GT.ZMAXSS) ZMAXSS=SS(J,I)
              END DO
            END DO
            OFFSS=ZMAXSS-ZMINSS
            ZMINSS=ZMINSS-OFFSS*.05
          END IF
          IF(NPLOT.EQ.1)THEN
            CALL PGVPORT(0.40,0.65,0.25,0.45)
          ELSEIF(NPLOT.EQ.2)THEN
            CALL PGVPORT(0.40,0.65,0.05,0.25)
          ELSE
            STOP 'FATAL ERROR in subroutine PLOTSS: invalid NPLOT.'
          END IF
          CALL PGWINDOW(X(1)-.05,X(NX)+(X(NX)-X(1))/3.+.05,
     +     ZMINSS,ZMAXSS+FACTOR*OFFSS)
          CALL PGSCH(0.7)
          CALL PGBOX(' ',0.0,0,'BTNS',0.0,0)
          CALL PGSCH(1.0)
        END IF
C
        DO I=1,NY
          DO J=1,NX
            IF(NY.EQ.1)THEN
              XP(J)=X(J)
              YP(J)=SS(J,I)
            ELSE
              XP(J)=X(J)+(Y(I)-Y(1))/(Y(NY)-Y(1))*(X(NX)-X(1))/3.
              YP(J)=SS(J,I)+(Y(I)-Y(1))/(Y(NY)-Y(1))*FACTOR*OFFSS
            END IF
          END DO
          CALL PGLINE(NX,XP,YP)
        END DO
        IF(NY.GT.1)THEN
          DO J=1,NX
            DO I=1,NY
              XP(I)=X(J)+(Y(I)-Y(1))/(Y(NY)-Y(1))*(X(NX)-X(1))/3.
              YP(I)=SS(J,I)+(Y(I)-Y(1))/(Y(NY)-Y(1))*FACTOR*OFFSS
            END DO
            CALL PGLINE(NY,XP,YP)
          END DO
        END IF
C dibujamos cruces en los puntos eliminados
        IF(LCOLOR) CALL PGSCI(2)
        DO J=1,NX
          DO I=1,NY
            IF(.NOT.IFPIXEL(J,I))THEN
              XX=X(J)+(Y(I)-Y(1))/(Y(NY)-Y(1))*(X(NX)-X(1))/3.
              YY=SS(J,I)+(Y(I)-Y(1))/(Y(NY)-Y(1))*FACTOR*OFFSS
              CALL PGPOINT(1,XX,YY,5)
            END IF
          END DO
        END DO
C
        END
C
C******************************************************************************
C Pide parametros para el ajuste de las superficies polinomicas. Si
C LDEF=.TRUE. se muestran valores por defecto.
        SUBROUTINE PIDEPOLY(NBINX,NBINY,BINMODE,BINSIGMA,GX,GY,LDEF)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INTEGER READILIM
        REAL READF
C
        INTEGER NMAXBINX,NMAXBINY            !valor maximo del binning en X e Y
        PARAMETER(NMAXBINX=101,NMAXBINY=101)   !modificar tambien en subrutinas
C
        INTEGER NBINX,NBINY
        INTEGER BINMODE
        REAL BINSIGMA
        INTEGER GX,GY
        LOGICAL LDEF
C
        CHARACTER*1 CMASK
        CHARACTER*50 CDUMMY
C
        COMMON/BLKPIDE1/NSCAN,NCHAN
        COMMON/BLKPIDE2/CMASK
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
        CDUMMY='@'
C
        IF(CMASK.EQ.'y')THEN
          NBINX=1
          NBINY=1
          GOTO 10
        END IF
C
C definimos el binning en la direccion X e Y para el ajuste inicial a
C superficies polinomicas de grado bajo
        WRITE(*,*)
        WRITE(*,101)'* Enter binning to fit initial '//
     +   'polynomial surface: '
        WRITE(*,100)'Binning in X-direction (odd) '
        IF(LDEF)WRITE(CDUMMY,*)NBINX
        IF(NCHAN.GT.NMAXBINX)THEN
          NBINX=READILIM(CDUMMY,1,NMAXBINX)
        ELSE
          NBINX=READILIM(CDUMMY,1,NCHAN)
        END IF
        WRITE(*,100)'Binniny in Y-direction (odd) '
        IF(LDEF)WRITE(CDUMMY,*)NBINY
        IF(NSCAN.GT.NMAXBINX)THEN
          NBINY=READILIM(CDUMMY,1,NMAXBINY)
        ELSE
          NBINY=READILIM(CDUMMY,1,NSCAN)
        END IF
C
10      CONTINUE
        WRITE(*,101)'* Enter binning mode:'
        WRITE(*,101)'  1=mean'
        WRITE(*,101)'  2=mean rejecting points with sigma'
        WRITE(*,101)'  3=median'
        WRITE(*,100)'Option '
        IF(LDEF)THEN
          WRITE(CDUMMY,*) BINMODE
          BINMODE=READILIM(CDUMMY,1,3)
        ELSE
          BINMODE=READILIM('3',1,3)
        END IF
        IF(BINMODE.EQ.2)THEN
          WRITE(*,100)'Times sigma to exclude points '
          WRITE(CDUMMY,*) BINSIGMA
          BINSIGMA=READF(CDUMMY)
        ELSE
          BINSIGMA=0.
        END IF
C
        WRITE(*,101)'* Enter polynomial degrees for surface fit:'
        WRITE(*,100)'Polynomial degree in X-direction '
        IF(LDEF)WRITE(CDUMMY,*)GX
        GX=READILIM(CDUMMY,0,19)
        WRITE(*,100)'Polynomial degree in Y-direction '
        IF(LDEF)WRITE(CDUMMY,*)GY
        GY=READILIM(CDUMMY,0,19)
C
100     FORMAT(A,$)
101     FORMAT(A)
        END
