C------------------------------------------------------------------------------
C Version 04-October-2007                                          file:plots.f
C @ ncl & fjg
C------------------------------------------------------------------------------
Comment
C
C Program: plots
C Classification: graphic display
C Description: General program to produce line plots and display images.
C
Comment
C
C Realiza plots (cortes en X e Y) de una imagen cualquiera. Permite
C dibujar la suma de varios scans/channels, asi como superponer plots
C
C El programa fue modificado para poder leer varios ficheros (que se van
C almacenando en diferentes buffers). A la hora de dibujar podemos elegir el
C buffer (o los buffers) a representar.
C
C NOTA: (2-Febrero-1995) descubro que al utilizar como numero de fichero 99 al
C hacer un OPEN,CLOSE, no funciona la salida grafica correctamente. Eso se
C debe a que ese numero debe ser el utilizado por PGPLOT para escribir
C al fichero PostScript.
C
C NOTA: (4-Octubre-2007) incluyo la opcion de leer trozos de un fichero
C especificando la region a leer (o a promediar) tras el nombre del
C fichero. Ejemplos:
C    fichero.u,3 -> lee solo el scan numero 3
C    fichero.u,14,18 -> lee una subimagen de 5 scans
C    fichero.u,14+18 -> lee un espectro promedio de los scans 14 a 18
C No funciona en la calculadora.
C
C NOTA: si NCMAX.LT.NSMAX hay que cambiar el dimensionado de XP,YP (los
C       arrays que contienen los puntos a dibujar), de modo que este
C       dimensionado siempre se realice con el numero mayor.
C
        PROGRAM PLOTS
        IMPLICIT NONE
C
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER TRUELEN
        INTEGER READILIM
        REAL READF
C
        INTEGER NBOTONES
        PARAMETER(NBOTONES=24)
C
        INTEGER I,J,L
        INTEGER NCBUFF                                 !numero de Buffer actual
        INTEGER NB,NBLOCAL,NBUFF
        INTEGER NN0,NN1,NN2
        INTEGER N1,N2
        INTEGER IXC1,IXC2
        INTEGER NS1_,NS2_
        INTEGER BMODE(NBOTONES),BMODE_0(NBOTONES)
        INTEGER MICOLOR(NMAXBUFF),MISLS(NMAXBUFF),MISLW(NMAXBUFF)
        INTEGER NSMAX_LOCAL,NCMAX_LOCAL
        INTEGER NTERM,IDN(MAX_ID_RED)
        INTEGER NSCAN2,NCHAN2
        INTEGER OLDCOLOR
        REAL A(NCMAX,NSMAX,NMAXBUFF)
        REAL XC,YC
        REAL XMIN,XMAX,YMIN,YMAX
        REAL XP(NCMAX),YP(NCMAX)
        REAL RVEL
        REAL STWV2,DISP2
        REAL ANGLE3D
        REAL YRMSTOL
        REAL PSEUDO_WEIGHT,PSEUDO_POWER
        REAL AMEAN,ASIGMA
        CHARACTER*1 CSEP_
        CHARACTER*1 CH
        CHARACTER*1 CTYPE,CFIT
        CHARACTER*50 CDUMMY
        CHARACTER*75 INFILE,INFILEBUFF(NMAXBUFF)
        CHARACTER*80 GLABEL
        CHARACTER*20 LABEL(NBOTONES),LABEL_0(NBOTONES)
        LOGICAL LCOLOR(MAX_ID_RED)
        LOGICAL LPLOTS
        LOGICAL LDEFBUFF(NMAXBUFF)        !indica si un buffer ha sido definido
        LOGICAL LUSEBUFF(NMAXBUFF)           !indica si un buffer es utilizable
        LOGICAL LOOP
        LOGICAL AUTO,OVER,SINGLE,ZOOM
        LOGICAL IFSCAN(NSMAX),IFCHAN(NCMAX)
        LOGICAL ASKXMIN,ASKXMAX,ASKYMIN,ASKYMAX
        LOGICAL NEWPLOT
        LOGICAL INSIDE
        LOGICAL LNORM
        LOGICAL LANYPLOT
        LOGICAL LBEXIST
        LOGICAL LOK
        LOGICAL PSEUDO_LUP
C
        COMMON/BLKNCBUFF/NCBUFF
        COMMON/BLKDATA1/NSCAN,NCHAN
        COMMON/BLKDATA2/A
        COMMON/BLKBUFF1/LDEFBUFF,LUSEBUFF
        COMMON/BLKBUFF2/INFILEBUFF
        COMMON/BLKTYPE/AUTO,OVER,SINGLE,ZOOM,CTYPE
        COMMON/BLKRANGE/NN0,NN1,NN2
        COMMON/BLKLIMITS/XMIN,XMAX,YMIN,YMAX
        COMMON/BLKASK/ASKXMIN,ASKXMAX,ASKYMIN,ASKYMAX
        COMMON/BLKIF/IFSCAN,IFCHAN
        COMMON/BLKCOLOR/MICOLOR,MISLS,MISLW,GLABEL
        COMMON/BLKLASER/LPLOTS
        COMMON/BLKPLOT/XP,YP
        COMMON/BLKLNORM/LNORM
        COMMON/BLKANGLE3D/ANGLE3D
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
C------------------------------------------------------------------------------
        DATA (LABEL_0(I),I=1,NBOTONES)/
     +   '[m]enu','[z]oom','[x] cut','[y] cut','ALL [@]','[n]ew File',
     +   '[a]uto','no-o[v]er','Xmin','Xmax','Ymin','Ymax',
     +   'cross[:]cut','[f]IT','no-no[r]m.',
     +     '[b]uffer','C:,S:,W:','[w]hole',
     +   '#[1]','#[2]','#[3]','#[4]','#[5]','#[6]'/
        DATA (BMODE_0(I),I=1,NBOTONES)/0,3,0,0,1,0, 
     +                                 1,3,0,0,0,0,
     +                                 0,3,1,0,0,0,
     +                                 3,3,3,3,3,3/
C------------------------------------------------------------------------------
        NEWPLOT=.TRUE. !avoid compilation warning
        OUTFILEX=OUTFILEX
C
        THISPROGRAM='plots'
        CALL WELCOME('04-October-2007')
C
        ANGLE3D=25.0
C
        NSMAX_LOCAL=NSMAX
        NCMAX_LOCAL=NCMAX
        IF(NSMAX_LOCAL.GT.NCMAX_LOCAL)STOP 'FATAL ERROR: NSMAX.GT.NCMAX'
C
        YRMSTOL=1.E-4
        PSEUDO_WEIGHT=100.0
        PSEUDO_POWER=2.0
        PSEUDO_LUP=.TRUE.
        RVEL=0.
C
        NCBUFF=0
        DO NBUFF=1,NMAXBUFF
          LDEFBUFF(NBUFF)=.FALSE.                              !Buffer definido
          LUSEBUFF(NBUFF)=.FALSE.                     !Buffer definido y usable
          DO I=1,LEN(INFILEBUFF(NBUFF))
            INFILEBUFF(NBUFF)(I:I)=' '
          END DO
        END DO
C
        LPLOTS=.TRUE.
C
        CALL RPGBEGIN(NTERM,IDN,LCOLOR)
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C Introducimos imagen inicial para tener ya definido NSCAN y NCHAN
        LOOP=.TRUE.
        DO WHILE(LOOP)
          WRITE(*,100)'Input file name'
          INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
          CALL HAYCOMA(INFILE,NSCAN,NS1_,NS2_,CSEP_)
          LOOP=(NS1_.EQ.0)
          IF(LOOP) CLOSE(20)
        END DO
        IF(CSEP_.EQ.'+')THEN
          NSCAN=1
        ELSE
          NSCAN=NS2_-NS1_+1
        END IF
        CALL NEWBUFF(NCBUFF)
        WRITE(CDUMMY,*)NCBUFF
        WRITE(*,100)'Buffer number '
        NCBUFF=READILIM(CDUMMY,1,NMAXBUFF)
        LDEFBUFF(NCBUFF)=.TRUE.
        IF(TRUELEN(OBJECT).GT.0)THEN
          INFILEBUFF(NCBUFF)=INFILE(1:TRUELEN(INFILE))//' ['//
     +     OBJECT(1:TRUELEN(OBJECT))//']'
        ELSE
          INFILEBUFF(NCBUFF)=INFILE(1:TRUELEN(INFILE))
        END IF
        LUSEBUFF(NCBUFF)=.TRUE.
        WRITE(*,100)'Reading file...'
        IF(NS1_.GT.1)THEN !ignoramos primeros NS1_-1 scans
          DO I=1,NS1_-1
            READ(20)(A(J,1,NCBUFF),J=1,NCHAN)
          END DO
        END IF
        DO I=NS1_,NS2_
          READ(20)(A(J,I-NS1_+1,NCBUFF),J=1,NCHAN)
        END DO
        CLOSE(20)
        IF(CSEP_.EQ.'+')THEN
          IF(NS2_.GT.NS1_)THEN
            DO J=1,NCHAN
              DO I=NS1_+1,NS2_
                A(J,1,NCBUFF)=A(J,1,NCBUFF)+A(J,I-NS1_+1,NCBUFF)
              END DO
              A(J,1,NCBUFF)=A(J,1,NCBUFF)/REAL(NS2_-NS1_+1)
            END DO
          END IF
        END IF
        WRITE(*,101)'OK! File read and closed.'
        CALL SHOWLIMITS(1,NSCAN,1,NCHAN,AMEAN,ASIGMA)
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C                              MENU INICIAL
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
3       CALL BUTTSYB(2)
C
4       DO I=1,NBOTONES
          CALL BUTTSEX(I,.FALSE.)
        END DO
        CALL BUTTON(1,'[p]lots',0)
        CALL BUTTON(2,'[l]ook',0)
        CALL BUTTON(3,'[n]ew File',0)
        CALL BUTTON(4,'[b]uffer',0)
        CALL BUTTON(5,'[m]ath',0)
        CALL BUTTON(6,'[q]uit!',0)
5       NB=0
        DO WHILE(NB.EQ.0)
          CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
          CALL IFBUTTON(XC,YC,NB)
          NBLOCAL=INDEX('plnbmq',CH)
          IF((NBLOCAL.NE.0).AND.(CH.NE.' '))THEN
            CALL BUTTQEX(NBLOCAL,LBEXIST)
            IF(LBEXIST) NB=NBLOCAL
          END IF
        END DO
C..............................................................................
        IF(NB.EQ.1)THEN
          CALL BUTTON(1,'[p]lots',5)
          IF(NCBUFF.EQ.0)THEN                        !tomamos el primero usable
            DO NBUFF=NMAXBUFF,1,-1
              IF(LUSEBUFF(NBUFF)) NCBUFF=NBUFF
            END DO
          END IF
          IF(NCBUFF.EQ.0)THEN   !si no hay usables, tomamos el primero definido
            DO NBUFF=NMAXBUFF,1,-1
              IF(LDEFBUFF(NBUFF)) NCBUFF=NBUFF
            END DO
            LUSEBUFF(NCBUFF)=.TRUE.
          END IF
          CALL PGERAS
          CALL BUTTSIT(.FALSE.)
          NEWPLOT=.TRUE.
          CALL BUTTSYB(4)
          GOTO 7
C..............................................................................
        ELSEIF(NB.EQ.2)THEN
          CALL BUTTON(2,'[l]ook',5)
          IF(NCBUFF.EQ.0)THEN                        !tomamos el primero usable
            DO NBUFF=NMAXBUFF,1,-1
              IF(LUSEBUFF(NBUFF)) NCBUFF=NBUFF
            END DO
          END IF
          IF(NCBUFF.EQ.0)THEN   !si no hay usables, tomamos el primero definido
            DO NBUFF=NMAXBUFF,1,-1
              IF(LDEFBUFF(NBUFF)) NCBUFF=NBUFF
            END DO
            LUSEBUFF(NCBUFF)=.TRUE.
          END IF
          CALL PGERAS
          CALL BUTTSIT(.TRUE.)
          CALL BUTTSYB(4)
          CALL SUBLOOK(STWV,DISP)
          CALL RPGERASB
          GOTO 3
C..............................................................................
        ELSEIF(NB.EQ.3)THEN
          CALL BUTTON(3,'[n]ew File',5)
          LOOP=.TRUE.
          DO WHILE(LOOP)
            WRITE(*,100)'Input file name'
            INFILE=INFILEX(20,'@',NSCAN2,NCHAN2,STWV2,DISP2,1,.FALSE.)
            CALL HAYCOMA(INFILE,NSCAN2,NS1_,NS2_,CSEP_)
            LOOP=(NS1_.EQ.0)
            IF(LOOP) CLOSE(20)
          END DO
          IF(CSEP_.EQ.'+')THEN
            NSCAN2=1
          ELSE
            NSCAN2=NS2_-NS1_+1
          END IF
          CALL CCSIZE(NSCAN,NCHAN,STWV,DISP,
     +     NSCAN2,NCHAN2,STWV2,DISP2,LOK)
          IF(.NOT.LOK)THEN
            CALL BUTTON(3,'[n]ew File',0)
            GOTO 4
          END IF
          CALL NEWBUFF(NCBUFF)
          WRITE(CDUMMY,*)NCBUFF
          WRITE(*,100)'Buffer number '
          NCBUFF=READILIM(CDUMMY,1,NMAXBUFF)
          LDEFBUFF(NCBUFF)=.TRUE.
          IF(TRUELEN(OBJECT).GT.0)THEN
            INFILEBUFF(NCBUFF)=INFILE(1:TRUELEN(INFILE))//' ['//
     +       OBJECT(1:TRUELEN(OBJECT))//']'
          ELSE
            INFILEBUFF(NCBUFF)=INFILE(1:TRUELEN(INFILE))
          END IF
          LUSEBUFF(NCBUFF)=.TRUE.
          WRITE(*,100)'Reading file...'
          IF(NS1_.GT.1)THEN !ignoramos primeros NS1_-1 scans
            DO I=1,NS1_-1
              READ(20)(A(J,1,NCBUFF),J=1,NCHAN)
            END DO
          END IF
          DO I=NS1_,NS2_
            READ(20)(A(J,I-NS1_+1,NCBUFF),J=1,NCHAN)
          END DO
          CLOSE(20)
          IF(CSEP_.EQ.'+')THEN
            IF(NS2_.GT.NS1_)THEN
              DO J=1,NCHAN
                DO I=NS1_+1,NS2_
                  A(J,1,NCBUFF)=A(J,1,NCBUFF)+A(J,I-NS1_+1,NCBUFF)
                END DO
                A(J,1,NCBUFF)=A(J,1,NCBUFF)/REAL(NS2_-NS1_+1)
              END DO
            END IF
          END IF
          WRITE(*,101)'OK! File read and closed.'
          CALL SHOWLIMITS(1,NSCAN,1,NCHAN,AMEAN,ASIGMA)
          CALL BUTTON(3,'[n]ew File',0)
          GOTO 4
C..............................................................................
        ELSEIF(NB.EQ.4)THEN
          CALL BUTTON(4,'[b]uffer',5)
          CALL SHOWBUFF
          CALL BUTTON(4,'[b]uffer',0)
          GOTO 5
C..............................................................................
        ELSEIF(NB.EQ.5)THEN
          CALL BUTTON(5,'[m]ath',5)
          CALL SUBIMATH
          CALL BUTTON(5,'[m]ath',0)
          GOTO 3
C..............................................................................
        ELSEIF(NB.EQ.6)THEN
          CALL BUTTON(6,'[q]uit!',5)
          CALL MY_PGEND                          !close all the graphic devices
          STOP
C..............................................................................
        END IF
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C                                PLOTS
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
7       WRITE(*,101)'--------------------------------------------------'
        WRITE(*,101)'* Remember:'
        WRITE(*,101)'<+> in SINGLE mode plots next X or Y cut'
        WRITE(*,101)'<-> in SINGLE mode plots previous X or Y cut'
        WRITE(*,101)'<*> in SINGLE mode plots a sequence of X or Y cuts'
        WRITE(*,101)'<c> plot buffer buttons with current buffer colors'
        WRITE(*,101)'<s> give statistics of currently plotted data'
        WRITE(*,101)'<l> plot indices, typical lines, etc'
        WRITE(*,101)'--------------------------------------------------'
C
        DO I=1,NBOTONES
          LABEL(I)=LABEL_0(I)
          BMODE(I)=BMODE_0(I)
        END DO
C
        AUTO=.TRUE.
        OVER=.FALSE.
        SINGLE=.FALSE.
        ASKXMIN=.TRUE.
        ASKXMAX=.TRUE.
        ASKYMIN=.TRUE.
        ASKYMAX=.TRUE.
        DO NBUFF=1,NMAXBUFF
          MICOLOR(NBUFF)=NBUFF
          MISLS(NBUFF)=1
          MISLW(NBUFF)=1
        END DO
        GLABEL=' '
        LNORM=.FALSE.
        CFIT='1'
        LANYPLOT=.FALSE.            !determina si hay algun grafico dibujado
C
C actualizamos buffers segun esten definidos o no y en utilizacion o no
ccc10      DO NBUFF=1,NMAXBUFF
        DO NBUFF=1,NMAXBUFF
          IF(LDEFBUFF(NBUFF))THEN
            IF(LUSEBUFF(NBUFF))THEN
              IF(NBUFF.EQ.NCBUFF)THEN
                BMODE(NBUFF+18)=5              !boton apretado, texto en blanco
              ELSE
                BMODE(NBUFF+18)=1                              !texto en blanco
              END IF
            ELSE
              BMODE(NBUFF+18)=0                                 !texto en negro
            END IF
          END IF
        END DO
C
        WRITE(CDUMMY,'(A2,I2,A3,I1,A3,I2)')
     +   'C:',MICOLOR(NCBUFF),',S:',MISLS(NCBUFF),',W:',MISLW(NCBUFF)
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        LABEL(17)=CDUMMY(1:L)
C
        DO I=1,NBOTONES
          IF(BMODE(I).GE.0) CALL BUTTON(I,LABEL(I),0)
          IF(BMODE(I).GT.0) CALL BUTTON(I,LABEL(I),BMODE(I))
        END DO
C
        NB=0
C------------------------------------------------------------------------------
20      CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
        CALL IFBUTTON(XC,YC,NB)
        NBLOCAL=INDEX('mzxy@nav    :frb w123456',CH)
        IF((NBLOCAL.NE.0).AND.(CH.NE.' '))THEN
          CALL BUTTQEX(NBLOCAL,LBEXIST)
          IF(LBEXIST) NB=NBLOCAL
        END IF
C------------------------------------------------------------------------------
C opciones no establecidas con botones
C------------------------------------------------------------------------------
        IF(CH.EQ.'+')THEN
          IF(NCBUFF.EQ.0)THEN
            WRITE(*,101)'WARNING: select buffer number first!'
            GOTO 20
          END IF
          IF(NN0.EQ.0)THEN
            WRITE(*,101)'WARNING: plot single cut first!'
            GOTO 20
          END IF
          IF(SINGLE)THEN
            IF(CTYPE.EQ.'X')THEN
              IF(NSCAN.GT.1)THEN
                IFSCAN(NN0)=.FALSE.
                NN0=NN0+1
                IF(NN0.GT.NSCAN) NN0=1
                IFSCAN(NN0)=.TRUE.
                CALL SUBPLOT(STWV,DISP,.FALSE.)
                GOTO 20
              END IF
            ELSE !CTYPE.EQ.'Y'
              IF(NCHAN.GT.1)THEN
                IFCHAN(NN0)=.FALSE.
                NN0=NN0+1
                IF(NN0.GT.NCHAN) NN0=1
                IFCHAN(NN0)=.TRUE.
                CALL SUBPLOT(STWV,DISP,.FALSE.)
                GOTO 20
              END IF
            END IF
          END IF
C------------------------------------------------------------------------------
        ELSEIF(CH.EQ.'-')THEN
          IF(NCBUFF.EQ.0)THEN
            WRITE(*,101)'WARNING: select buffer number first!'
            GOTO 20
          END IF
          IF(NN0.EQ.0)THEN
            WRITE(*,101)'WARNING: plot single cut first!'
            GOTO 20
          END IF
          IF(SINGLE)THEN
            IF(CTYPE.EQ.'X')THEN
              IF(NSCAN.GT.1)THEN
                IFSCAN(NN0)=.FALSE.
                NN0=NN0-1
                IF(NN0.LT.1) NN0=NSCAN
                IFSCAN(NN0)=.TRUE.
                CALL SUBPLOT(STWV,DISP,.FALSE.)
                GOTO 20
              END IF
            ELSE !CTYPE.EQ.'Y'
              IF(NCHAN.GT.1)THEN
                IFCHAN(NN0)=.FALSE.
                NN0=NN0-1
                IF(NN0.LT.1) NN0=NCHAN
                IFCHAN(NN0)=.TRUE.
                CALL SUBPLOT(STWV,DISP,.FALSE.)
                GOTO 20
              END IF
            END IF
          END IF
C------------------------------------------------------------------------------
        ELSEIF(CH.EQ.'*')THEN
          IF(NCBUFF.EQ.0)THEN
            WRITE(*,101)'WARNING: select buffer number first!'
            GOTO 20
          END IF
          IF(NN0.EQ.0)THEN
            WRITE(*,101)'WARNING: plot single cut first!'
            GOTO 20
          END IF
          IF(SINGLE)THEN
            IF(CTYPE.EQ.'X')THEN
              IF(NSCAN.GT.1)THEN
                DO I=1,NSCAN
                  IFSCAN(NN0)=.FALSE.
                  NN0=NN0+1
                  IF(NN0.GT.NSCAN) NN0=1
                  IFSCAN(NN0)=.TRUE.
                  CALL SUBPLOT(STWV,DISP,.FALSE.)
                END DO
                GOTO 20
              END IF
            ELSE !CTYPE.EQ.'Y'
              IF(NCHAN.GT.1)THEN
                DO J=1,NCHAN
                  IFCHAN(NN0)=.FALSE.
                  NN0=NN0+1
                  IF(NN0.GT.NCHAN) NN0=1
                  IFCHAN(NN0)=.TRUE.
                  CALL SUBPLOT(STWV,DISP,.FALSE.)
                END DO
                GOTO 20
              END IF
            END IF
          END IF
C------------------------------------------------------------------------------
        ELSEIF(CH.EQ.'c')THEN
          IF(NCBUFF.EQ.0)THEN
            WRITE(*,101)'WARNING: select buffer number first!'
            GOTO 20
          END IF
          DO NBUFF=1,NMAXBUFF
            IF(LDEFBUFF(NBUFF))THEN
              IF(LUSEBUFF(NBUFF))THEN
                CALL PGQCI(OLDCOLOR)
                CALL BUTTON(NBUFF+18,LABEL(NBUFF+18),-MICOLOR(NBUFF)-1)
                CALL PGSCI(OLDCOLOR)
              END IF
            END IF
          END DO
          CALL PGSCI(1)
          GOTO 20
C------------------------------------------------------------------------------
        ELSEIF(CH.EQ.'s')THEN
          IF(NCBUFF.EQ.0)THEN
            WRITE(*,101)'WARNING: select buffer number first!'
            GOTO 20
          END IF
          CALL SUBPLOT(STWV,DISP,.TRUE.)
          GOTO 20
C------------------------------------------------------------------------------
        ELSEIF(CH.EQ.'l')THEN
          IF((STWV.NE.0.).AND.(DISP.NE.0.).AND.(CTYPE.EQ.'X'))THEN
            IF(NCBUFF.EQ.0)THEN
              WRITE(*,101)'WARNING: select buffer number first!'
            ELSE
              WRITE(*,100)'Radial velocity (km/sec) '
              WRITE(CDUMMY,*) RVEL
              RVEL=READF(CDUMMY)
              CALL SUBPMORE(STWV,DISP,RVEL)
            END IF
          END IF
          GOTO 20
C------------------------------------------------------------------------------
        END IF
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
          IF(NB.EQ.0)THEN
            IF(.NOT.INSIDE(XC,YC))THEN
              WRITE(*,101)'ERROR: cursor out of plot.'
            ELSE
              WRITE(CDUMMY,*)NINT(XC),',',YC
              CALL RMBLANK(CDUMMY,CDUMMY,L)
              IF((STWV.NE.0.).AND.(DISP.NE.0.).AND.(CTYPE.EQ.'X'))THEN
                WRITE(*,100)'Cursor at (x,y): '//CDUMMY(1:L)
                WRITE(CDUMMY,*)STWV+DISP*(XC-1.)
                CALL RMBLANK(CDUMMY,CDUMMY,L)
                WRITE(*,101)'    ---> wavelength: '//CDUMMY(1:L)
              ELSE
                WRITE(*,101)'Cursor at (x,y): '//CDUMMY(1:L)
              END IF
            END IF
C------------------------------------------------------------------------------
          ELSEIF(NB.EQ.1)THEN
            CALL RPGERASB
            GOTO 3
C------------------------------------------------------------------------------
          ELSEIF(NB.EQ.2)THEN
            IF(NCBUFF.EQ.0)THEN
              WRITE(*,101)'WARNING: select buffer number first!'
              GOTO 20
            END IF
            CALL BUTTON(2,LABEL(2),5)
            WRITE(*,101)'Press mouse button...'
            IF(LCOLOR(1)) CALL PGSCI(5)
            CALL RPGBAND(6,0,0.,0.,XC,YC,CH)
            IF(LCOLOR(1)) CALL PGSCI(1)
            IF((CH.NE.'A').AND.(CH.NE.'D').AND.(CH.NE.'X')) THEN
              WRITE(*,101)'ERROR: mouse buttom has not been detected.'
              CALL BUTTON(2,LABEL(2),BMODE(2))
              GOTO 20
            END IF
            IXC1=INT(XC+0.5)
            IF(IXC1.LT.NN1) IXC1=NN1
            IF(IXC1.GT.NN2) IXC1=NN2
            WRITE(*,100)'Point #1, cursor at x= '
            WRITE(*,*)INT(XC+0.5)
            WRITE(*,101)'Press mouse button...'
            IF(LCOLOR(1)) CALL PGSCI(5)
            CALL RPGBAND(4,0,REAL(IXC1),0.,XC,YC,CH)
            IF(LCOLOR(1)) CALL PGSCI(1)
            IF((CH.NE.'A').AND.(CH.NE.'D').AND.(CH.NE.'X')) THEN
              WRITE(*,101)'ERROR: mouse buttom has not been detected.'
              CALL BUTTON(2,LABEL(2),BMODE(2))
              GOTO 20
            END IF
            IXC2=INT(XC+0.5)
            IF(IXC2.LT.NN1) IXC2=NN1
            IF(IXC2.GT.NN2) IXC2=NN2
            WRITE(*,100)'Point #2, cursor at x= '
            WRITE(*,*)INT(XC+0.5)
            NN1=MIN0(IXC1,IXC2)
            NN2=MAX0(IXC1,IXC2)
            CALL BUTTON(2,LABEL(2),BMODE(2))
            ZOOM=.TRUE.
            CALL SUBPLOT(STWV,DISP,.FALSE.)
            ZOOM=.FALSE.
C------------------------------------------------------------------------------
          ELSEIF(NB.EQ.3)THEN
            IF(NCBUFF.EQ.0)THEN
              WRITE(*,101)'WARNING: select buffer number first!'
              GOTO 20
            END IF
            CALL BUTTON(3,LABEL(3),5)
            NN1=1
            NN2=NCHAN
            CTYPE='X'
            IF(NSCAN.EQ.1)THEN
              IFSCAN(1)=.TRUE.
              NN0=1
              GOTO 32
            END IF
            IF(SINGLE)THEN
              WRITE(*,100)'Scan number '
              NN0=READILIM('@',1,NSCAN)
            ELSE
              NN0=0
              IF(LABEL(5).EQ.'ALL [@]')THEN
                DO I=1,NSCAN
                  IFSCAN(I)=.TRUE.
                END DO
              ELSE
                DO I=1,NSCAN
                  IFSCAN(I)=.FALSE.
                END DO
30              WRITE(*,'(A,I5,A)')'(valid range: 1,',NSCAN,')'
                WRITE(*,100)'First and last scan (0,0=EXIT) '
                CALL READ2I('0,0',N1,N2)
                IF((N1.EQ.0).AND.(N2.EQ.0))GOTO 31
                IF((N1.LT.1).OR.(N2.GT.NSCAN).OR.(N1.GT.N2))THEN
                  WRITE(*,101)'Invalid numbers. Try again.'
                  GOTO 30
                END IF
                DO I=N1,N2
                  IFSCAN(I)=.TRUE.
                END DO
                GOTO 30
31              CONTINUE
              END IF
            END IF
32          CALL SUBPLOT(STWV,DISP,.FALSE.)
            LANYPLOT=.TRUE.
            IF(.NOT.OVER)THEN
              BMODE(2)=0
              CALL BUTTON(2,LABEL(2),BMODE(2))
            END IF
            CALL BUTTON(3,LABEL(3),0)
            IF(NEWPLOT)THEN
              NEWPLOT=.FALSE.
              BMODE(8)=0
              CALL BUTTON(8,LABEL(8),BMODE(8))
              CALL BUTTON(8,LABEL(8),1)
              CALL BUTTON(14,LABEL(14),0)
            END IF
C------------------------------------------------------------------------------
          ELSEIF(NB.EQ.4)THEN
            IF(NCBUFF.EQ.0)THEN
              WRITE(*,101)'WARNING: select buffer number first!'
              GOTO 20
            END IF
            CALL BUTTON(4,LABEL(4),5)
            NN1=1
            NN2=NSCAN
            CTYPE='Y'
            IF(NCHAN.EQ.1)THEN
              IFSCAN(1)=.TRUE.
              NN0=1
              GOTO 42
            END IF
            IF(SINGLE)THEN
              WRITE(*,100)'Channel number '
              NN0=READILIM('@',1,NCHAN)
            ELSE
              NN0=0
              IF(LABEL(5).EQ.'ALL [@]')THEN
                DO I=1,NCHAN
                  IFCHAN(I)=.TRUE.
                END DO
              ELSE
                DO I=1,NCHAN
                  IFCHAN(I)=.FALSE.
                END DO
40              WRITE(*,'(A,I5,A)')'(valid range: 1,',NCHAN,')'
                WRITE(*,100)'First and last channel (0,0=EXIT) '
                CALL READ2I('0,0',N1,N2)
                IF((N1.EQ.0).AND.(N2.EQ.0))GOTO 41
                IF((N1.LT.1).OR.(N2.GT.NCHAN).OR.(N1.GT.N2))THEN
                  WRITE(*,101)'Invalid numbers. Try again.'
                  GOTO 40
                END IF
                DO I=N1,N2
                  IFCHAN(I)=.TRUE.
                END DO
                GOTO 40
41              CONTINUE
              END IF
            END IF
42          CALL SUBPLOT(STWV,DISP,.FALSE.)
            LANYPLOT=.TRUE.
            IF(.NOT.OVER)THEN
              BMODE(2)=0
              CALL BUTTON(2,LABEL(2),BMODE(2))
            END IF
            CALL BUTTON(4,LABEL(4),0)
            IF(NEWPLOT)THEN
              NEWPLOT=.FALSE.
              BMODE(8)=0
              CALL BUTTON(8,LABEL(8),BMODE(8))
              CALL BUTTON(8,LABEL(8),1)
              CALL BUTTON(14,LABEL(14),0)
            END IF
C------------------------------------------------------------------------------
          ELSEIF(NB.EQ.5)THEN
            IF(LABEL(5).EQ.'SINGLE [@]')THEN
              LABEL(5)='ADDED [@]'
              SINGLE=.FALSE.
            ELSEIF(LABEL(5).EQ.'ADDED [@]')THEN
              LABEL(5)='ALL [@]'
            ELSEIF(LABEL(5).EQ.'ALL [@]')THEN
              LABEL(5)='SINGLE [@]'
              SINGLE=.TRUE.
            END IF
            CALL BUTTON(5,LABEL(5),0)
            CALL BUTTON(5,LABEL(5),1)
C------------------------------------------------------------------------------
          ELSEIF(NB.EQ.6)THEN
            CALL BUTTON(6,LABEL(6),5)
            LOOP=.TRUE.
            DO WHILE(LOOP)
              WRITE(*,100)'Input file name'
              INFILE=INFILEX(20,'@',NSCAN2,NCHAN2,STWV2,DISP2,1,.FALSE.)
              CALL HAYCOMA(INFILE,NSCAN2,NS1_,NS2_,CSEP_)
              LOOP=(NS1_.EQ.0)
              IF(LOOP) CLOSE(20)
            END DO
            IF(CSEP_.EQ.'+')THEN
              NSCAN2=1
            ELSE
              NSCAN2=NS2_-NS1_+1
            END IF
            CALL CCSIZE(NSCAN,NCHAN,STWV,DISP,
     +       NSCAN2,NCHAN2,STWV2,DISP2,LOK)
            IF(.NOT.LOK)THEN
              CALL BUTTON(6,LABEL(6),0)
              GOTO 20
            END IF
            BMODE(18+NCBUFF)=1                       !ya no es imagen principal
            CALL BUTTON(18+NCBUFF,LABEL(18+NCBUFF),0)
            CALL BUTTON(18+NCBUFF,LABEL(18+NCBUFF),1)
            CALL NEWBUFF(NCBUFF)
            WRITE(CDUMMY,*)NCBUFF
            WRITE(*,100)'Buffer number '
            NCBUFF=READILIM(CDUMMY,1,NMAXBUFF)
            LDEFBUFF(NCBUFF)=.TRUE.
            IF(TRUELEN(OBJECT).GT.0)THEN
              INFILEBUFF(NCBUFF)=INFILE(1:TRUELEN(INFILE))//' ['//
     +         OBJECT(1:TRUELEN(OBJECT))//']'
            ELSE
              INFILEBUFF(NCBUFF)=INFILE(1:TRUELEN(INFILE))
            END IF
            LUSEBUFF(NCBUFF)=.TRUE.
            WRITE(*,100)'Reading file...'
            IF(NS1_.GT.1)THEN !ignoramos primeros NS1_-1 scans
              DO I=1,NS1_-1
                READ(20)(A(J,1,NCBUFF),J=1,NCHAN)
              END DO
            END IF
            DO I=NS1_,NS2_
              READ(20)(A(J,I-NS1_+1,NCBUFF),J=1,NCHAN)
            END DO
            CLOSE(20)
            IF(CSEP_.EQ.'+')THEN
              IF(NS2_.GT.NS1_)THEN
                DO J=1,NCHAN
                  DO I=NS1_+1,NS2_
                  A(J,1,NCBUFF)=A(J,1,NCBUFF)+A(J,I-NS1_+1,NCBUFF)
                  END DO
                  A(J,1,NCBUFF)=A(J,1,NCBUFF)/REAL(NS2_-NS1_+1)
                END DO
              END IF
            END IF
            WRITE(*,101)'OK! File read and closed.'
            CALL BUTTON(18+NCBUFF,LABEL(18+NCBUFF),0)
            IF(LUSEBUFF(NCBUFF))THEN
              CALL BUTTON(18+NCBUFF,LABEL(18+NCBUFF),5)
            END IF
            CALL SHOWLIMITS(1,NSCAN,1,NCHAN,AMEAN,ASIGMA)
            CALL BUTTON(6,LABEL(6),0)
C------------------------------------------------------------------------------
          ELSEIF(NB.EQ.7)THEN
            IF(AUTO)THEN
              IF(ASKXMIN.AND.ASKXMAX.AND.ASKYMIN.AND.ASKYMAX)THEN
                LABEL(7)='no[a]uto'
                AUTO=.FALSE.
              END IF
            ELSE
              LABEL(7)='[a]uto'
              AUTO=.TRUE.
            END IF
            CALL BUTTON(7,LABEL(7),0)
            CALL BUTTON(7,LABEL(7),1)
            ASKXMIN=.TRUE.
            ASKXMAX=.TRUE.
            ASKYMIN=.TRUE.
            ASKYMAX=.TRUE.
            DO I=9,12
              BMODE(I)=0
              CALL BUTTON(I,LABEL(I),BMODE(I))
            END DO
C------------------------------------------------------------------------------
          ELSEIF(NB.EQ.8)THEN
            IF(LABEL(8).EQ.'no-o[v]er')THEN
              LABEL(8)='o[v]er'
              OVER=.TRUE.
              BMODE(2)=3
              CALL BUTTON(2,LABEL(2),BMODE(2))
              BMODE(7)=3
              CALL BUTTON(7,LABEL(7),BMODE(7))
              DO I=9,12
                CALL BUTTON(I,LABEL(I),3)
              END DO
              BMODE(18)=3
              CALL BUTTON(18,LABEL(18),BMODE(18))
            ELSE
              LABEL(8)='no-o[v]er'
              OVER=.FALSE.
              BMODE(2)=0
              CALL BUTTON(2,LABEL(2),BMODE(2))
              BMODE(7)=0
              CALL BUTTON(7,LABEL(7),BMODE(7))
              CALL BUTTON(7,LABEL(7),1)
              DO I=9,12
                CALL BUTTON(I,LABEL(I),BMODE(I))
              END DO
              BMODE(18)=0
              CALL BUTTON(18,LABEL(18),BMODE(18))
            END IF
            CALL BUTTON(8,LABEL(8),0)
            CALL BUTTON(8,LABEL(8),1)
C------------------------------------------------------------------------------
          ELSEIF(NB.EQ.9)THEN
            CALL BUTTON(9,LABEL(9),5)
            WRITE(*,100)'Enter Xmin'
            XMIN=READF('@')
            ASKXMIN=.FALSE.
            BMODE(9)=1
            CALL BUTTON(9,LABEL(9),0)
            CALL BUTTON(9,LABEL(9),1)
C------------------------------------------------------------------------------
          ELSEIF(NB.EQ.10)THEN
            CALL BUTTON(10,LABEL(10),5)
            WRITE(*,100)'Enter Xmax'
            XMAX=READF('@')
            ASKXMAX=.FALSE.
            BMODE(10)=1
            CALL BUTTON(10,LABEL(10),0)
            CALL BUTTON(10,LABEL(10),1)
C------------------------------------------------------------------------------
          ELSEIF(NB.EQ.11)THEN
            CALL BUTTON(11,LABEL(11),5)
            WRITE(*,100)'Enter Ymin'
            YMIN=READF('@')
            ASKYMIN=.FALSE.
            BMODE(11)=1
            CALL BUTTON(11,LABEL(11),0)
            CALL BUTTON(11,LABEL(11),1)
C------------------------------------------------------------------------------
          ELSEIF(NB.EQ.12)THEN
            CALL BUTTON(12,LABEL(12),5)
            WRITE(*,100)'Enter Ymax'
            YMAX=READF('@')
            ASKYMAX=.FALSE.
            BMODE(12)=1
            CALL BUTTON(12,LABEL(12),0)
            CALL BUTTON(12,LABEL(12),1)
C------------------------------------------------------------------------------
          ELSEIF(NB.EQ.13)THEN
            CALL BUTTON(13,LABEL(13),5)
            IF(NCBUFF.EQ.0)THEN
              WRITE(*,101)'WARNING: select buffer number first!'
              CALL BUTTON(13,LABEL(13),0)
              GOTO 20
            END IF
            IF(LANYPLOT)THEN                     !si hay algun grafico dibujado
              IF((NSCAN.GT.1).AND.(NCHAN.GT.1))THEN    !si no, no tiene sentido
                WRITE(*,101)'Press mouse button...'
                IF(LCOLOR(1)) CALL PGSCI(5)
                CALL RPGBAND(6,0,0.,0.,XC,YC,CH)
                IF(LCOLOR(1)) CALL PGSCI(1)
                IF((CH.NE.'A').AND.(CH.NE.'D').AND.(CH.NE.'X')) THEN
                  WRITE(*,101)'ERROR: mouse buttom has not been '//
     +             'detected.'
                  CALL BUTTON(13,LABEL(13),0)
                  GOTO 20
                END IF
                IXC1=INT(XC+0.5)
                IF(IXC1.LT.NN1) IXC1=NN1
                IF(IXC1.GT.NN2) IXC1=NN2
                WRITE(*,110)'Point #1, cursor at x= ',IXC1
                WRITE(*,101)'Press mouse button...'
                IF(LCOLOR(1)) CALL PGSCI(5)
                CALL RPGBAND(4,0,REAL(IXC1),0.,XC,YC,CH)
                IF(LCOLOR(1)) CALL PGSCI(1)
                IF((CH.NE.'A').AND.(CH.NE.'D').AND.(CH.NE.'X')) THEN
                  WRITE(*,101)'ERROR: mouse buttom has not been '//
     +             'detected.'
                  CALL BUTTON(13,LABEL(13),0)
                  GOTO 20
                END IF
                IXC2=INT(XC+0.5)
                IF(IXC2.LT.NN1) IXC2=NN1
                IF(IXC2.GT.NN2) IXC2=NN2
                WRITE(*,110)'Point #2, cursor at x= ',IXC2
                NN1=MIN0(IXC1,IXC2)
                NN2=MAX0(IXC1,IXC2)
                IF(LABEL(5).EQ.'SINGLE [@]')THEN
                  LABEL(5)='ADDED [@]'
                  SINGLE=.FALSE.
                  CALL BUTTON(5,LABEL(5),0)
                  CALL BUTTON(5,LABEL(5),1)
                ELSEIF(LABEL(5).EQ.'ADDED [@]')THEN
                ELSEIF(LABEL(5).EQ.'ALL [@]')THEN
                  LABEL(5)='ADDED [@]'
                  CALL BUTTON(5,LABEL(5),0)
                  CALL BUTTON(5,LABEL(5),1)
                END IF
                IF(CTYPE.EQ.'X')THEN
                  CTYPE='Y'
                  DO J=1,NCHAN
                    IFCHAN(J)=.FALSE.
                  END DO
                  DO J=NN1,NN2
                    IFCHAN(J)=.TRUE.
                  END DO
                ELSE !CTYPE.EQ.'Y'
                  CTYPE='X'
                  DO J=1,NSCAN
                    IFSCAN(J)=.FALSE.
                  END DO
                  DO J=NN1,NN2
                    IFSCAN(J)=.TRUE.
                  END DO
                END IF
                NN1=1
                IF(CTYPE.EQ.'X')THEN
                  NN2=NCHAN
                ELSE
                  NN2=NSCAN
                END IF
                CALL SUBPLOT(STWV,DISP,.FALSE.)
                CALL BUTTON(13,LABEL(13),0)
                GOTO 20
              END IF
            ELSE
              CALL BUTTON(13,LABEL(13),0)
              GOTO 20
            END IF
C------------------------------------------------------------------------------
          ELSEIF(NB.EQ.14)THEN
            IF(NCBUFF.EQ.0)THEN
              WRITE(*,101)'WARNING: select buffer number first!'
              GOTO 20
            END IF
            CALL BUTTON(14,LABEL(14),5)
            WRITE(*,101)'---------------------------------'
            WRITE(*,101)'(1) fit polynomial'
            WRITE(*,101)'(2) fit gaussian+continuum'
            WRITE(*,101)'(x) fit gaussian+continuum (cte)'
            WRITE(*,101)'(3) fit gaussian+cte'
            WRITE(*,101)'(c) fit Cauchy function+continuum'
            WRITE(*,101)'(p) fit pseudo-continuum'
            WRITE(*,101)'---------------------------------'
            WRITE(*,100)'Option '
            CFIT(1:1)=READC(CFIT,'123cpx')
            CALL SUBFIT(CFIT,NCHAN,STWV,DISP,YRMSTOL,
     +       PSEUDO_WEIGHT,PSEUDO_POWER,PSEUDO_LUP)
            CALL BUTTON(14,LABEL(14),0)
C------------------------------------------------------------------------------
          ELSEIF(NB.EQ.15)THEN
            IF(LABEL(15).EQ.'no[r]m.')THEN
              LABEL(15)='no-no[r]m.'
              LNORM=.FALSE.
              CALL BUTTON(15,LABEL(15),0)
              CALL BUTTON(15,LABEL(15),1)
            ELSE
              LABEL(15)='no[r]m.'
              LNORM=.TRUE.
              CALL BUTTON(15,LABEL(15),0)
              CALL BUTTON(15,LABEL(15),1)
            END IF
C------------------------------------------------------------------------------
          ELSEIF(NB.EQ.16)THEN
            CALL BUTTON(16,LABEL(16),5)
            CALL SHOWBUFF
            CALL BUTTON(16,LABEL(16),0)
C------------------------------------------------------------------------------
          ELSEIF(NB.EQ.17)THEN
            CALL BUTTON(17,LABEL(17),5)
            WRITE(*,*)
            WRITE(*,101)'--------------------------------------------'
            WRITE(*,101)'                 COLOR'
            WRITE(*,101)'--------------------------------------------'
            WRITE(*,101)' 0: black        1: white       2: red'
            WRITE(*,101)' 3: green        4: blue        5: cyan'
            WRITE(*,101)' 6: magenta      7: yellow      8: orange'
            WRITE(*,101)' 9: 3+7         10: 3+5        11: 4+5'
            WRITE(*,101)'12: grey        13: light grey 14: dark grey'
            WRITE(*,101)'--------------------------------------------'
            WRITE(*,'(A,I1,A,$)')'Buffer #',NCBUFF,' >>'
ccc50          WRITE(*,100)'New color (PGPLOT number) '
            WRITE(*,100)'New color (PGPLOT number) '
            WRITE(CDUMMY,*)MICOLOR(NCBUFF)
            MICOLOR(NCBUFF)=READILIM(CDUMMY,0,14)
            WRITE(*,*)
            WRITE(*,101)'--------------------------------------------'
            WRITE(*,101)'               LINE STYLE'
            WRITE(*,101)'--------------------------------------------'
            WRITE(*,101)'1: full line                2: dashed'
            WRITE(*,101)'3: dot-dash-dot-dash        4: dotted'
            WRITE(*,101)'5: dash-dot-dot-dot'
            WRITE(*,101)'--------------------------------------------'
            WRITE(*,'(A,I1,A,$)')'Buffer #',NCBUFF,' >>'
ccc51          WRITE(*,100)'New line style (PGPLOT number) '
            WRITE(*,100)'New line style (PGPLOT number) '
            WRITE(CDUMMY,*)MISLS(NCBUFF)
            MISLS(NCBUFF)=READILIM(CDUMMY,1,5)
            WRITE(*,*)
            WRITE(*,101)'--------------------------------------------'
            WRITE(*,101)'               LINE WIDTH'
            WRITE(*,101)'--------------------------------------------'
            WRITE(*,101)'Line  width  is  specified by  the number of'
            WRITE(*,101)'strokes  to  be  used  (in  the range: 1-21)'
            WRITE(*,101)'--------------------------------------------'
            WRITE(*,'(A,I1,A,$)')'Buffer #',NCBUFF,' >>'
ccc52          WRITE(*,100)'New line width (PGPLOT number) '
            WRITE(*,100)'New line width (PGPLOT number) '
            WRITE(CDUMMY,*)MISLW(NCBUFF)
            MISLW(NCBUFF)=READILIM(CDUMMY,1,21)
C
            WRITE(CDUMMY,'(A2,I2,A3,I1,A3,I2)')'C:',MICOLOR(NCBUFF),
     +       ',S:',MISLS(NCBUFF),',W:',MISLW(NCBUFF)
            CALL RMBLANK(CDUMMY,CDUMMY,L)
            LABEL(17)=CDUMMY(1:L)
            CALL BUTTON(17,LABEL(17),0)
C------------------------------------------------------------------------------
          ELSEIF(NB.EQ.18)THEN
            IF(NCBUFF.EQ.0)THEN
              WRITE(*,101)'WARNING: select buffer number first!'
              GOTO 20
            END IF
            NN1=1
            IF(CTYPE.EQ.'X')THEN
              NN2=NCHAN
            ELSE
              NN2=NSCAN
            END IF
            CALL SUBPLOT(STWV,DISP,.FALSE.)
C------------------------------------------------------------------------------
          ELSEIF((NB.GE.19).AND.(NB.LE.24))THEN
            IF(NB-18.EQ.NCBUFF)THEN                  !si es NCBUFF se desactiva
              LUSEBUFF(NB-18)=.FALSE.
              CALL BUTTON(NB,LABEL(NB),0)
              NCBUFF=0
            ELSE                                               !si no es NCBUFF
              IF(LUSEBUFF(NB-18))THEN      !si se usaba, se convierte en NCBUFF
                IF(NCBUFF.NE.0)THEN     !si existe, se quita el anterior NCBUFF
                  CALL BUTTON(NCBUFF+18,LABEL(NCBUFF+18),0)
                  CALL BUTTON(NCBUFF+18,LABEL(NCBUFF+18),1)
                END IF
                NCBUFF=NB-18
                CALL BUTTON(NB,LABEL(NB),5)
              ELSE                      !si no se usaba, se convierte en usable
                LUSEBUFF(NB-18)=.TRUE.
                CALL BUTTON(NB,LABEL(NB),0)
                CALL BUTTON(NB,LABEL(NB),1)
              END IF
            END IF
C actualizamos botones
            IF(NCBUFF.EQ.0)THEN
              WRITE(CDUMMY,'(A2,I2,A3,I1,A3,I2)')
     +         'C:',0,',S:',0,',W:',0
              CALL RMBLANK(CDUMMY,CDUMMY,L)
              LABEL(17)=CDUMMY(1:L)
              BMODE(17)=3
              CALL BUTTON(17,LABEL(17),0)
              CALL BUTTON(17,LABEL(17),3)
            ELSE
              WRITE(CDUMMY,'(A2,I2,A3,I1,A3,I2)')'C:',MICOLOR(NCBUFF),
     +         ',S:',MISLS(NCBUFF),',W:',MISLW(NCBUFF)
              CALL RMBLANK(CDUMMY,CDUMMY,L)
              LABEL(17)=CDUMMY(1:L)
              BMODE(17)=0
              CALL BUTTON(17,LABEL(17),0)
            END IF
C------------------------------------------------------------------------------
          END IF
C------------------------------------------------------------------------------
        GOTO 20
C
100     FORMAT(A,$)
101     FORMAT(A)
110     FORMAT(A,I6)
        END
C
C******************************************************************************
C
        SUBROUTINE SHOWLIMITS(NS1,NS2,NC1,NC2,AMEAN,ASIGMA)
        IMPLICIT NONE
        INTEGER NS1,NS2,NC1,NC2                             !limites a recorrer
        REAL AMEAN,ASIGMA
C
        INCLUDE 'redlib.inc'
C
        INTEGER NCBUFF
        INTEGER I,J,L,LL
        INTEGER IMIN,IMAX,JMIN,JMAX
        INTEGER NPIX
        REAL A(NCMAX,NSMAX,NMAXBUFF)
        REAL DATAMIN,DATAMAX
        CHARACTER*50 CDUMMY
        DOUBLE PRECISION MEAN,SIGMA
C
        COMMON/BLKNCBUFF/NCBUFF
        COMMON/BLKDATA1/NSCAN,NCHAN
        COMMON/BLKDATA2/A
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        IMIN=NS1
        IMAX=NS2
        JMIN=NC1
        JMAX=NC2
        DATAMIN=A(NC1,NS1,NCBUFF)
        DATAMAX=A(NC1,NS1,NCBUFF)
        MEAN=0.D0
        DO I=NS1,NS2
          DO J=NC1,NC2
            MEAN=MEAN+DBLE(A(J,I,NCBUFF))
            IF(A(J,I,NCBUFF).LT.DATAMIN)THEN
              IMIN=I
              JMIN=J
              DATAMIN=A(J,I,NCBUFF)
            END IF
            IF(A(J,I,NCBUFF).GT.DATAMAX)THEN
              IMAX=I
              JMAX=J
              DATAMAX=A(J,I,NCBUFF)
            END IF
          END DO
        END DO
        NPIX=(NS2-NS1+1)*(NC2-NC1+1)
        MEAN=MEAN/DBLE(NPIX)
        SIGMA=0.D0
        DO I=NS1,NS2
          DO J=NC1,NC2
            SIGMA=SIGMA+(DBLE(A(J,I,NCBUFF))-MEAN)*
     +                  (DBLE(A(J,I,NCBUFF))-MEAN)
          END DO
        END DO
        SIGMA=DSQRT(SIGMA/DBLE(NPIX-1))
        WRITE(*,*)
        WRITE(*,'(A,I5,A,I5)')'> From Scan    #',NS1,' to ',NS2
        WRITE(*,'(A,I5,A,I5)')'> From Channel #',NC1,' to ',NC2
        WRITE(*,'(A,I10)')'> Total number of pixels: ',NPIX
        WRITE(CDUMMY,*)DATAMIN
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,'(A,A,$)')'> Minimum: ',CDUMMY(1:L)
        DO LL=1,15-L
          WRITE(*,100)' '
        END DO
        WRITE(*,'(A,I5,A,I5,A)')' in pixel: ',JMIN,',',IMIN,
     +   '  (channel,scan)'
        WRITE(CDUMMY,*)DATAMAX
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,'(A,A,$)')'> Maximum: ',CDUMMY(1:L)
        DO LL=1,15-L
          WRITE(*,100)' '
        END DO
        WRITE(*,'(A,I5,A,I5,A)')' in pixel: ',JMAX,',',IMAX,
     +   '  (channel,scan)'
        WRITE(*,100)'> Mean   : '
        WRITE(*,*)REAL(MEAN)
        WRITE(*,100)'> Sigma  : '
        WRITE(*,*)REAL(SIGMA)
        WRITE(*,*)
        AMEAN=REAL(MEAN)
        ASIGMA=REAL(SIGMA)
100     FORMAT(A,$)
        END
C
C******************************************************************************
C
        SUBROUTINE SUBPLOT(STWV,DISP,LSTAT)
        IMPLICIT NONE
C
        INCLUDE 'redlib.inc'
        INTEGER TRUELEN
        REAL READF
C
        LOGICAL LSTAT                 !indica si hay que hacer o no estadistica
C
        INTEGER NCBUFF,NBUFF
        INTEGER NPTOS
        INTEGER NN0,NN1,NN2
        INTEGER I,J,L
        INTEGER NADDED
        INTEGER MICOLOR(NMAXBUFF),MISLS(NMAXBUFF),MISLW(NMAXBUFF)
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        REAL XP(NCMAX),YP(NCMAX),YPBUFF(NCMAX,NMAXBUFF)
        REAL A(NCMAX,NSMAX,NMAXBUFF)
        REAL XMIN,XMAX,YMIN,YMAX,DX,DY
        REAL XMINL,XMAXL
        REAL FACTNORM(NMAXBUFF)
        CHARACTER*1 CTYPE
        CHARACTER*50 CDUMMY
        CHARACTER*75 INFILEBUFF(NMAXBUFF)
        CHARACTER*80 GLABEL,YLABEL
        LOGICAL LCOLOR(MAX_ID_RED)
        LOGICAL LDEFBUFF(NMAXBUFF)        !indica si un buffer ha sido definido
        LOGICAL LUSEBUFF(NMAXBUFF)           !indica si un buffer es utilizable
        LOGICAL LPLOTS
        LOGICAL AUTO,OVER,SINGLE,ZOOM
        LOGICAL ASKXMIN,ASKXMAX,ASKYMIN,ASKYMAX
        LOGICAL IFSCAN(NSMAX),IFCHAN(NCMAX)
        LOGICAL LW,LNORM,FIRSTYP
C
        COMMON/BLKNCBUFF/NCBUFF
        COMMON/BLKBUFF1/LDEFBUFF,LUSEBUFF
        COMMON/BLKBUFF2/INFILEBUFF
        COMMON/BLKDATA1/NSCAN,NCHAN
        COMMON/BLKDATA2/A
        COMMON/BLKTYPE/AUTO,OVER,SINGLE,ZOOM,CTYPE
        COMMON/BLKRANGE/NN0,NN1,NN2
        COMMON/BLKLIMITS/XMIN,XMAX,YMIN,YMAX
        COMMON/BLKASK/ASKXMIN,ASKXMAX,ASKYMIN,ASKYMAX
        COMMON/BLKIF/IFSCAN,IFCHAN
        COMMON/BLKCOLOR/MICOLOR,MISLS,MISLW,GLABEL
        COMMON/BLKLASER/LPLOTS
        COMMON/BLKPLOT/XP,YP
        COMMON/BLKLNORM/LNORM
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
C------------------------------------------------------------------------------
        IF((STWV.NE.0.).AND.(DISP.NE.0.).AND.(CTYPE.EQ.'X'))THEN
          LW=.TRUE.
        ELSE
          LW=.FALSE.
        END IF
        IF(OVER)THEN
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGSCI(0)
            CALL PGMTEXT('T',2.5,1.,1.,GLABEL)
            CALL PGSCI(1)
          END DO
        END IF
C
        IF(SINGLE)THEN
          NADDED=1
          IF(CTYPE.EQ.'X')THEN
            WRITE(GLABEL,'(A,I5,A)')'(scan #',NN0,')'
          ELSE
            WRITE(GLABEL,'(A,I5,A)')'(channel #',NN0,')'
          END IF
        ELSE
          NADDED=0
          IF(CTYPE.EQ.'X')THEN
            DO I=1,NSCAN
              IF(IFSCAN(I)) NADDED=NADDED+1
            END DO
            IF(NADDED.EQ.0)THEN
              WRITE(*,101)'ERROR: No. added scans = 0'
              RETURN
            END IF
            WRITE(GLABEL,'(A,I5,A)')'(#',NADDED,' scans added)'
          ELSE
            DO I=1,NCHAN
              IF(IFCHAN(I)) NADDED=NADDED+1
            END DO
            IF(NADDED.EQ.0)THEN
              WRITE(*,101)'ERROR: No. added channels = 0'
              RETURN
            END IF
            WRITE(GLABEL,'(A,I5,A)')'(#',NADDED,' channels added)'
          END IF
        END IF
C------------------------------------------------------------------------------
        DO I=NN1,NN2
          XP(I-NN1+1)=REAL(I)
        END DO
C------------------------------------------------------------------------------
        IF(CTYPE.EQ.'X')THEN                                      !CTYPE.EQ.'X'
          DO NBUFF=1,NMAXBUFF
            IF(LUSEBUFF(NBUFF))THEN
              IF(SINGLE)THEN
                DO I=NN1,NN2
                  YPBUFF(I-NN1+1,NBUFF)=A(I,NN0,NBUFF)
                END DO
              ELSE
                DO I=NN1,NN2
                  YPBUFF(I-NN1+1,NBUFF)=0.
                END DO
                DO I=NN1,NN2
                  DO J=1,NSCAN
                    IF(IFSCAN(J)) YPBUFF(I-NN1+1,NBUFF)=
     +               YPBUFF(I-NN1+1,NBUFF)+A(I,J,NBUFF)
                  END DO
                END DO
              END IF
            END IF
          END DO
C------------------------------------------------------------------------------
        ELSE                                                      !CTYPE.EQ.'Y'
          DO NBUFF=1,NMAXBUFF
            IF(LUSEBUFF(NBUFF))THEN
              IF(SINGLE)THEN
                DO I=NN1,NN2
                  YPBUFF(I-NN1+1,NBUFF)=A(NN0,I,NBUFF)
                END DO
              ELSE
                DO I=NN1,NN2
                  YPBUFF(I-NN1+1,NBUFF)=0.
                END DO
                DO I=NN1,NN2
                  DO J=1,NCHAN
                    IF(IFCHAN(J)) YPBUFF(I-NN1+1,NBUFF)=
     +               YPBUFF(I-NN1+1,NBUFF)+A(J,I,NBUFF)
                  END DO
                END DO
              END IF
            END IF
          END DO
C------------------------------------------------------------------------------
        END IF
C------------------------------------------------------------------------------
        NPTOS=NN2-NN1+1
C------------------------------------------------------------------------------
        IF(LNORM)THEN
          DO NBUFF=1,NMAXBUFF
            IF(LUSEBUFF(NBUFF))THEN
              FACTNORM(NBUFF)=0.
              DO I=NN1,NN2
                FACTNORM(NBUFF)=FACTNORM(NBUFF)+YPBUFF(I-NN1+1,NBUFF)
              END DO
              FACTNORM(NBUFF)=FACTNORM(NBUFF)/REAL(NN2-NN1+1)
              IF(FACTNORM(NBUFF).EQ.0.0)THEN
                WRITE(*,100)'WARNING: normalization factor=0 in'
                WRITE(*,'(A,I1)')' Buffer #',NBUFF
                WRITE(*,101)'> Normalization factor fixed to 1.'
                WRITE(*,*)
                FACTNORM(NBUFF)=1.
              END IF
            END IF
          END DO
        ELSE
          DO NBUFF=1,NMAXBUFF
            IF(LUSEBUFF(NBUFF)) FACTNORM(NBUFF)=1.
          END DO
        END IF
C
        DO NBUFF=1,NMAXBUFF
          IF(LUSEBUFF(NBUFF))THEN
            DO I=1,NPTOS
              YPBUFF(I,NBUFF)=YPBUFF(I,NBUFF)/REAL(NADDED)
              YPBUFF(I,NBUFF)=YPBUFF(I,NBUFF)/FACTNORM(NBUFF)
            END DO
          END IF
        END DO
C------------------------------------------------------------------------------
        IF(OVER)THEN
          DO NBUFF=1,NMAXBUFF
            IF(LUSEBUFF(NBUFF))THEN
              DO I=1,NPTOS
                YP(I)=YPBUFF(I,NBUFF)
              END DO
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                IF(LCOLOR(ITERM)) CALL PGSCI(MICOLOR(NBUFF))
                CALL PGSLS(MISLS(NBUFF))
                CALL PGSLW(MISLW(NBUFF))
                CALL PGBIN(NPTOS,XP,YP,.TRUE.)
              END DO
              GLABEL='(Last '//GLABEL(2:)
            END IF
          END DO
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
            CALL PGSLS(1)
            CALL PGSLW(1)
            CALL PGMTEXT('T',2.5,1.,1.,GLABEL)
          END DO
C..............................................................................
        ELSE
          IF((AUTO).OR.(ZOOM))THEN
c.............
            IF(ASKXMIN)THEN
              XMIN=XP(1)
              DO I=2,NPTOS
                IF(XP(I).LT.XMIN) XMIN=XP(I)
              END DO
            END IF
c......
            IF(ASKXMAX)THEN
              XMAX=XP(1)
              DO I=2,NPTOS
                IF(XP(I).GT.XMAX) XMAX=XP(I)
              END DO
            END IF
c......
            DX=XMAX-XMIN
            IF(ASKXMIN)THEN
              XMIN=XMIN-DX*.05
            END IF
            IF(ASKXMAX)THEN
              XMAX=XMAX+DX*.05
            END IF
c......
            IF(ASKYMIN)THEN 
              FIRSTYP=.TRUE.
              DO NBUFF=1,NMAXBUFF
                IF(LUSEBUFF(NBUFF))THEN
                  DO I=1,NPTOS
                    YP(I)=YPBUFF(I,NBUFF)
                  END DO
                  IF(FIRSTYP)THEN
                    YMIN=YP(1)
                    FIRSTYP=.FALSE.
                  END IF
                  DO I=1,NPTOS
                    IF(YP(I).LT.YMIN) YMIN=YP(I)
                  END DO
                END IF
              END DO
            END IF
c......
            IF(ASKYMAX)THEN
              FIRSTYP=.TRUE.
              DO NBUFF=1,NMAXBUFF
                IF(LUSEBUFF(NBUFF))THEN
                  DO I=1,NPTOS
                    YP(I)=YPBUFF(I,NBUFF)
                  END DO
                  IF(FIRSTYP)THEN
                    YMAX=YP(1)
                    FIRSTYP=.FALSE.
                  END IF
                  DO I=1,NPTOS
                    IF(YP(I).GT.YMAX) YMAX=YP(I)
                  END DO
                END IF
              END DO
            END IF
c......
            DY=YMAX-YMIN
            IF(ASKYMIN)THEN
              YMIN=YMIN-DY*.05
            END IF
            IF(ASKYMAX)THEN
              YMAX=YMAX+DY*.05
            END IF
c.............
          ELSE
            IF(ASKXMIN)THEN
              WRITE(*,100)'Enter Xmin'
              XMIN=READF('@')
            END IF
            IF(ASKXMAX)THEN
              WRITE(*,100)'Enter Xmax'
              XMAX=READF('@')
            END IF
            IF(ASKYMIN)THEN
              WRITE(*,100)'Enter Ymin'
              YMIN=READF('@')
            END IF
            IF(ASKYMAX)THEN
              WRITE(*,100)'Enter Ymax'
              YMAX=READF('@')
            END IF
          END IF
          IF(XMIN.GE.XMAX)THEN
            XMIN=XMIN-1.0
            XMAX=XMAX+1.0
            WRITE(*,101)'WARNING: XMIN.GE.XMAX in SUBPLOT.'
          END IF
          IF(YMIN.GE.YMAX)THEN
            WRITE(*,101)'WARNING: YMIN.GE.YMAX in SUBPLOT.'
            YMIN=YMIN-1.0
            YMAX=YMAX+1.0
          END IF
          IF(LPLOTS)THEN
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              IF(ITERM.EQ.1)THEN
                CALL RPGERASW(0.,1.,0.,0.80)
              ELSE
              END IF
            END DO
            IF(LW)THEN
              XMINL=(XMIN-1.)*DISP+STWV
              XMAXL=(XMAX-1.)*DISP+STWV
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                IF(ITERM.EQ.1)THEN
                  CALL RPGENV(XMIN,XMAX,YMIN,YMAX,0,-2)
                ELSE
                  CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,-2)
                END IF
                CALL PGSLW(3)
                CALL PGSWIN(XMINL,XMAXL,YMIN,YMAX)
                CALL PGBOX('CMTS',0.,0,'BCNST',0.,0)
                CALL PGSWIN(XMIN,XMAX,YMIN,YMAX)
                CALL PGBOX('BNTS',0.,0,' ',0.,0)
                CALL PGSLW(1)
              END DO
            ELSE
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                CALL PGSLW(3)
                IF(ITERM.EQ.1)THEN
                  CALL RPGENV(XMIN,XMAX,YMIN,YMAX,0,0)
                ELSE
                  CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,0)
                END IF
                CALL PGBOX('M',0.,0,'M',0.,0)
                CALL PGSLW(1)
              END DO
            END IF
          ELSE
            IF(LW)THEN
              XMINL=(XMIN-1.)*DISP+STWV
              XMAXL=(XMAX-1.)*DISP+STWV
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                CALL PGSLW(3)
                CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,-2)
                CALL PGSWIN(XMINL,XMAXL,YMIN,YMAX)
                CALL PGBOX('CMTS',0.,0,'BCNST',0.,0)
                CALL PGSWIN(XMIN,XMAX,YMIN,YMAX)
                CALL PGBOX('BNTS',0.,0,' ',0.,0)
                CALL PGSLW(1)
              END DO
            ELSE
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                CALL PGSLW(3)
                CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,0)
                CALL PGBOX('M',0.,0,' ',0.,0)
                CALL PGSLW(1)
              END DO
            END IF
          END IF
C
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGIDEN_RED
          END DO
C
          IF(LNORM)THEN
            YLABEL='No. counts (normalized between '
            WRITE(CDUMMY,*)NN1
            CALL RMBLANK(CDUMMY,CDUMMY,L)
            YLABEL=YLABEL(1:TRUELEN(YLABEL))//' '//
     +       CDUMMY(1:L)//'-'
            WRITE(CDUMMY,*)NN2
            CALL RMBLANK(CDUMMY,CDUMMY,L)
            YLABEL=YLABEL(1:TRUELEN(YLABEL))//CDUMMY(1:L)//')'
          ELSE
            IF(SINGLE)THEN
              YLABEL='No. counts'
            ELSE
              YLABEL='No. counts (averaged)'
            END IF
          END IF
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGSLW(3)
            IF(CTYPE.EQ.'X')THEN
              CALL PGLABEL('channel',YLABEL,CHAR(32))
            ELSE
              CALL PGLABEL('scan',YLABEL,CHAR(32))
            END IF
            CALL PGMTEXT('T',2.5,1.,1.,GLABEL)
            IF(NCBUFF.NE.0) CALL PGMTEXT('T',2.5,0.,0.,'file: '//
     +       INFILEBUFF(NCBUFF)(1:TRUELEN(INFILEBUFF(NCBUFF))))
            CALL PGSLW(1)
          END DO
C
          DO NBUFF=1,NMAXBUFF
            IF(LUSEBUFF(NBUFF))THEN
              DO I=1,NPTOS
                YP(I)=YPBUFF(I,NBUFF)
              END DO
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                IF(LCOLOR(ITERM)) CALL PGSCI(MICOLOR(NBUFF))
                CALL PGSLS(MISLS(NBUFF))
                CALL PGSLW(MISLW(NBUFF))
                CALL PGBIN(NPTOS,XP,YP,.TRUE.)
              END DO
              IF(LSTAT)THEN
                WRITE(*,130)
                L=TRUELEN(INFILEBUFF(NBUFF))
                WRITE(*,'(A,I1,A,$)')'> Current buffer is #',NBUFF,': '
                WRITE(*,101)INFILEBUFF(NBUFF)(1:L)
                CALL GIVESTAT(NPTOS,XP,YP)
              END IF
            END IF
          END DO
          IF(LSTAT) WRITE(*,130)
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
            CALL PGSLS(1)
            CALL PGSLW(1)
          END DO
C..............................................................................
        END IF
C------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
130     FORMAT(79('-'))
        END
C
C******************************************************************************
C
        SUBROUTINE SUBFIT(CFIT,NCHAN,STWV,DISP,YRMSTOL,
     +   PSEUDO_WEIGHT,PSEUDO_POWER,PSEUDO_LUP)
        IMPLICIT NONE
        CHARACTER*1 CFIT
        INCLUDE 'redlib.inc'
        INCLUDE 'futils.inc'
        INTEGER READILIM
        REAL READF
C
        INTEGER NCBUFF
        INTEGER I,K
        INTEGER NP,NTERMS
        INTEGER NN0,NN1,NN2
        INTEGER MICOLOR(NMAXBUFF),MISLS(NMAXBUFF),MISLW(NMAXBUFF)
        INTEGER MICOLOR0
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        REAL XP(NCMAX),YP(NCMAX)
        REAL XF(NCMAX),YF(NCMAX),SIGMAY(NCMAX)
        REAL XX(1000),YY(1000)
        REAL XMIN,XMAX,YMIN,YMAX
        REAL A(20),CHISQR,X,POL
        REAL X0,SIGMA,AMP,Y0
        REAL EX0,ESIGMA,EAMP
        REAL EEX0,EESIGMA,EEAMP,EEY0
        REAL YRMSTOL
        REAL PSEUDO_WEIGHT,PSEUDO_POWER
        REAL XFMIN,XFMAX
        CHARACTER*1 CMODE,COK,CLUP
        CHARACTER*50 CDUMMY
        CHARACTER*80 GLABEL
        LOGICAL LCOLOR(MAX_ID_RED)
        LOGICAL PSEUDO_LUP
        COMMON/BLKNCBUFF/NCBUFF
        COMMON/BLKPLOT/XP,YP
        COMMON/BLKFITG1/NP
        COMMON/BLKFITG2/XF,YF
        COMMON/BLKCOLOR/MICOLOR,MISLS,MISLW,GLABEL
        COMMON/BLKRANGE/NN0,NN1,NN2
        COMMON/BLKLIMITS/XMIN,XMAX,YMIN,YMAX
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        MICOLOR0=MICOLOR(NCBUFF)+1
        CMODE='?'
C
        DO I=1,NCHAN
          SIGMAY(I)=0.
        END DO
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM)) CALL PGSCI(MICOLOR(NCBUFF))
          CALL PGBIN(NN2-NN1+1,XP,YP,.TRUE.)
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
        END DO
C
C------------------------------------------------------------------------------
        IF((CFIT.EQ.'1').OR.(CFIT.EQ.'p'))THEN
          CALL SELPFIT(NCHAN,CMODE)
          IF(NP.EQ.0)THEN
            WRITE(*,101)'ERROR: Number of points for fit = 0'
            RETURN
          END IF
          WRITE(*,110)'Number of points for fit = ',NP
10        WRITE(*,100)'Polynomial degree (-1=EXIT) '
          NTERMS=READILIM('-1',-1,19)
          IF(NTERMS.EQ.-1) RETURN
          IF(NTERMS.GT.NP-1)THEN
            WRITE(*,101)'ERROR: not enough points for this degree.'
            GOTO 10
          END IF
          NTERMS=NTERMS+1
          IF(CFIT.EQ.'1')THEN
            CALL POLFIT(XF,YF,SIGMAY,NP,NTERMS,0,A,CHISQR)
          ELSEIF(CFIT.EQ.'p')THEN
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
            CALL PSEUDOFIT(XF,YF,NP,NTERMS,YRMSTOL,
     +       PSEUDO_WEIGHT,PSEUDO_POWER,PSEUDO_LUP,A)
          ELSE
            WRITE(*,101) 'CFIT='//CFIT
            STOP 'FATAL ERROR: invalid CFIT option in SUBFIT.'
          END IF
          WRITE(*,*)
          WRITE(*,101)'Coefficients from fit: (y=a0+a1*x+a2*x*x+...)'
          DO I=1,NTERMS
            IF(I.LT.11)THEN
              WRITE(*,'(A,I1,A,$)')'> a(',I-1,') : '
            ELSE
              WRITE(*,'(A,I2,A,$)')'> a(',I-1,'): '
            END IF
            WRITE(*,*)A(I)
          END DO
          DO I=1,1000
            X=REAL(I-1)/999.*(XMAX-XMIN)+XMIN
            POL=A(NTERMS)
            DO K=NTERMS-1,1,-1
              POL=POL*X+A(K)
            END DO
            XX(I)=X
            YY(I)=POL
          END DO
          MICOLOR0=MICOLOR0+1
          IF(MICOLOR0.EQ.15) MICOLOR0=1
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(MICOLOR0)
            CALL PGLINE(1000,XX,YY)
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END DO
          GOTO 10
C------------------------------------------------------------------------------
        ELSEIF((CFIT.EQ.'2').OR.(CFIT.EQ.'x').OR.(CFIT.EQ.'c'))THEN
C caculamos el continuo
          WRITE(*,*)
          WRITE(*,101)'Enter points to fit continuum:'
          CALL SELPFIT(NCHAN,CMODE)
          IF(NP.EQ.0)THEN
            WRITE(*,101)'ERROR: Number of points for fit = 0'
            RETURN
          END IF
          WRITE(*,110)'Number of points for fit = ',NP
          IF(CFIT.EQ.'x')THEN
            NTERMS=0
            GOTO 22
          END IF
20        WRITE(*,100)'Polynomial degree '
          NTERMS=READILIM('1',0,19)
          IF(NTERMS.GT.NP-1)THEN
            WRITE(*,101)'ERROR: not enough points for this degree.'
            GOTO 20
          END IF
          IF((NTERMS.LT.0).OR.(NTERMS.GT.19))THEN
            WRITE(*,110)'ERROR: polynomial degree must be < 20.'
            GOTO 20
          END IF
22        NTERMS=NTERMS+1
          CALL POLFIT(XF,YF,SIGMAY,NP,NTERMS,0,A,CHISQR)
          WRITE(*,*)
          WRITE(*,101)'Coefficients from fit: (y=a0+a1*x+a2*x*x+...)'
          DO I=1,NTERMS
            IF(I.LT.10)THEN
              WRITE(*,'(A,I1,A,$)')'> a(',I,') : '
            ELSE
              WRITE(*,'(A,I2,A,$)')'> a(',I,'): '
            END IF
            WRITE(*,*)A(I)
          END DO
C dibujamos solo en la region de puntos seleccionada para hacer el ajuste
          XFMIN=XF(1)
          XFMAX=XFMIN
          DO I=2,NP
            IF(XF(I).LT.XFMIN) XFMIN=XF(I)
            IF(XF(I).GT.XFMAX) XFMAX=XF(I)
          END DO
          DO I=1,1000
            X=REAL(I-1)/999.*(XFMAX-XFMIN)+XFMIN
            POL=A(NTERMS)
            DO K=NTERMS-1,1,-1
              POL=POL*X+A(K)
            END DO
            XX(I)=X
            YY(I)=POL
          END DO
          MICOLOR0=MICOLOR0+1
          IF(MICOLOR0.EQ.15) MICOLOR0=1
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(MICOLOR0)
            CALL PGLINE(1000,XX,YY)
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END DO
          IF(CFIT.NE.'x')THEN
            WRITE(*,100)'Is this fit ok (y/n) '
            COK(1:1)=READC('y','yn')
            IF(COK.EQ.'n') GOTO 20
          END IF
C calculamos la funcion (gaussiana/cauchy)
          WRITE(*,*)
          IF((CFIT.EQ.'2').OR.(CFIT.EQ.'x'))THEN
            WRITE(*,101)'Enter points to fit gaussian:'
          ELSEIF(CFIT.EQ.'c')THEN
            WRITE(*,101)'Enter points to fit Cauchy function:'
          END IF
          CALL SELPFIT(NCHAN,CMODE)
          IF(NP.EQ.0)THEN
            WRITE(*,101)'ERROR: Number of points for fit = 0'
            RETURN
          END IF
          WRITE(*,110)'Number of points for fit = ',NP
          IF(NP.LT.3)THEN
            WRITE(*,101)'ERROR: No. of points for fit < 3'
            WRITE(*,*)
            RETURN
          END IF
C restamos el continuo calculado con el polinomio
          DO I=1,NP
            POL=A(NTERMS)
            DO K=NTERMS-1,1,-1
              POL=POL*XF(I)+A(K)
            END DO
            YF(I)=YF(I)-POL
          END DO
          IF((CFIT.EQ.'2').OR.(CFIT.EQ.'x'))THEN
            WRITE(CDUMMY,*) YRMSTOL
            WRITE(*,100)'YRMSTOL for DOWNHILL '
            YRMSTOL=READF(CDUMMY)
            CALL GAUSSFIT(NP,XF,YF,YF,X0,SIGMA,AMP,
     +       EX0,ESIGMA,EAMP,EEX0,EESIGMA,EEAMP,YRMSTOL,0)
            WRITE(*,101)'Coefficients from fit: '//
     +       'y=a*exp[-(x-x0)^2/(2 sig^2)]'
            WRITE(*,100)'> a   (+ rmsDOWNHILL): '
            WRITE(*,*)AMP,EEAMP
            WRITE(*,100)'> x0  (+ rmsDOWNHILL): '
            WRITE(*,*)X0,EEX0
            WRITE(*,100)'> sig (+ rmsDOWNHILL): '
            WRITE(*,*)SIGMA,EESIGMA
            WRITE(*,100)' x0 - 1 sig, x0 + 1 sig : '
            WRITE(*,*)X0-SIGMA,X0+SIGMA
            WRITE(*,100)' x0 - 2 sig, x0 + 2 sig : '
            WRITE(*,*)X0-2.*SIGMA,X0+2.*SIGMA
            WRITE(*,100)' x0 - 3 sig, x0 + 3 sig : '
            WRITE(*,*)X0-3.*SIGMA,X0+3.*SIGMA
            WRITE(*,100)' x0 - 4 sig, x0 + 4 sig : '
            WRITE(*,*)X0-4.*SIGMA,X0+4.*SIGMA
            WRITE(*,100)' x0 - 6 sig, x0 + 6 sig : '
            WRITE(*,*)X0-6.*SIGMA,X0+6.*SIGMA
            WRITE(*,100)' x0 - 8 sig, x0 + 8 sig : '
            WRITE(*,*)X0-8.*SIGMA,X0+8.*SIGMA
            IF((STWV.NE.0.).AND.(DISP.NE.0.))THEN
              WRITE(*,100)'Wavelength at x=x0 : '
              WRITE(*,*)STWV+DISP*(X0-1.)
            END IF
          ELSEIF(CFIT.EQ.'c')THEN
            WRITE(CDUMMY,*) YRMSTOL
            WRITE(*,100)'YRMSTOL for DOWNHILL '
            YRMSTOL=READF(CDUMMY)
            CALL CAUCHYFIT(X0,SIGMA,AMP,EEX0,EESIGMA,EEAMP,YRMSTOL)
            WRITE(*,101)'Coefficients from fit: '//
     +       'y=a/[sig^2 + (x-x0)^2]'
            WRITE(*,100)'> a   (+rms DOWNHILL): '
            WRITE(*,*)AMP
            WRITE(*,100)'> x0  (+rms DOWNHILL): '
            WRITE(*,*)X0
            WRITE(*,100)'> sig (+rms DOWNHILL): '
            WRITE(*,*)SIGMA
        print*,x0,sigma
          END IF
C dibujamos solo en la region de puntos seleccionada para hacer el ajuste
C incluyendo las regiones del continuo
          DO I=1,NP
            IF(XF(I).LT.XFMIN) XFMIN=XF(I)
            IF(XF(I).GT.XFMAX) XFMAX=XF(I)
          END DO
          DO I=1,1000
            XX(I)=REAL(I-1)/999.*(XFMAX-XFMIN)+XFMIN
          END DO
          IF((CFIT.EQ.'2').OR.(CFIT.EQ.'x'))THEN
            DO I=1,1000
              YY(I)=AMP*EXP(-(XX(I)-X0)*(XX(I)-X0)/(2.*SIGMA*SIGMA))
            END DO
          ELSEIF(CFIT.EQ.'c')THEN
            DO I=1,1000
              YY(I)=AMP/(SIGMA*SIGMA+(XX(I)-X0)*(XX(I)-X0))
            END DO
          END IF
          DO I=1,1000
            POL=A(NTERMS)
            DO K=NTERMS-1,1,-1
              POL=POL*XX(I)+A(K)
            END DO
            YY(I)=YY(I)+POL
          END DO
          MICOLOR0=MICOLOR0+1
          IF(MICOLOR0.EQ.15) MICOLOR0=1
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(MICOLOR0)
            CALL PGLINE(1000,XX,YY)
            IF((CFIT.EQ.'2').OR.(CFIT.EQ.'x'))THEN
              IF(LCOLOR(ITERM)) CALL PGSCI(MICOLOR0+1)
              CALL PGSLS(2)
              CALL PGMOVE(X0-SIGMA,YMIN)
              CALL PGDRAW(X0-SIGMA,YMAX)
              CALL PGMOVE(X0+SIGMA,YMIN)
              CALL PGDRAW(X0+SIGMA,YMAX)
              CALL PGMOVE(X0-2.*SIGMA,YMIN)
              CALL PGDRAW(X0-2.*SIGMA,YMAX)
              CALL PGMOVE(X0+2.*SIGMA,YMIN)
              CALL PGDRAW(X0+2.*SIGMA,YMAX)
              CALL PGMOVE(X0-3.*SIGMA,YMIN)
              CALL PGDRAW(X0-3.*SIGMA,YMAX)
              CALL PGMOVE(X0+3.*SIGMA,YMIN)
              CALL PGDRAW(X0+3.*SIGMA,YMAX)
              CALL PGSLS(1)
            END IF
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END DO
C------------------------------------------------------------------------------
        ELSEIF(CFIT.EQ.'3')THEN
          WRITE(*,*)
          WRITE(*,101)'Enter points to fit gaussian+cte:'
30        CALL SELPFIT(NCHAN,CMODE)
          WRITE(*,110)'Number of points for fit = ',NP
          IF(NP.LT.4)THEN
            WRITE(*,101)'ERROR: No. of points for fit < 4'
            WRITE(*,*)
            RETURN
          END IF
          WRITE(CDUMMY,*) YRMSTOL
          WRITE(*,100)'YRMSTOL for DOWNHILL '
          YRMSTOL=READF(CDUMMY)
          CALL GAUSCFIT(X0,SIGMA,AMP,Y0,EEX0,EESIGMA,EEAMP,EEY0,YRMSTOL)
          WRITE(*,101)'Coefficients from fit: '//
     +     'y=y0+a*exp[-(x-x0)^2/(2 sig^2)]'
          WRITE(*,100)'> a   (+rms DOWNHILL): '
          WRITE(*,*)AMP,EEAMP
          WRITE(*,100)'> x0  (+rms DOWNHILL): '
          WRITE(*,*)X0,EEX0
          WRITE(*,100)'> sig (+rms DOWNHILL): '
          WRITE(*,*)SIGMA,EESIGMA
          WRITE(*,100)'> y0  (+rms DOWNHILL): '
          WRITE(*,*) Y0,EEY0
          IF((STWV.NE.0.).AND.(DISP.NE.0.))THEN
            WRITE(*,100)'Wavelength at x=x0 : '
            WRITE(*,*)STWV+DISP*(X0-1.)
          END IF
        print*,x0,sigma
C dibujamos solo en la region de puntos seleccionada para hacer el ajuste
C incluyendo las regiones del continuo
          XFMIN=XF(1)
          XFMAX=XFMIN
          DO I=2,NP
            IF(XF(I).LT.XFMIN) XFMIN=XF(I)
            IF(XF(I).GT.XFMAX) XFMAX=XF(I)
          END DO
          DO I=1,1000
            XX(I)=REAL(I-1)/999.*(XFMAX-XFMIN)+XFMIN
            YY(I)=AMP*EXP(-(XX(I)-X0)*(XX(I)-X0)/(2.*SIGMA*SIGMA))+Y0
          END DO
          MICOLOR0=MICOLOR0+1
          IF(MICOLOR0.EQ.15) MICOLOR0=1
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(MICOLOR0)
            CALL PGLINE(1000,XX,YY)
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END DO
          GOTO 30
C------------------------------------------------------------------------------
        ELSE
          WRITE(*,101)'FATAL ERROR: imposible option in subroutine '//
     +     'SUBFIT.'
          STOP
        END IF
100     FORMAT(A,$)
101     FORMAT(A)
110     FORMAT(A,I6)
        END
C
C******************************************************************************
C Selecciona, de los puntos dibujados XP(),YP(), aquellos que se utilizaran
C para realizar un ajuste XF(),YF(). NP es el numero de puntos tomados para
C calcular dicho ajuste.
        SUBROUTINE SELPFIT(NCHAN,CMODE)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'futils.inc'
        CHARACTER*1 CMODE
C
        INTEGER NCBUFF
        INTEGER I,L
        INTEGER NP
        INTEGER NC0,NC1,NC2
        INTEGER NN0,NN1,NN2
        INTEGER MICOLOR(NMAXBUFF),MISLS(NMAXBUFF),MISLW(NMAXBUFF)
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        REAL XC,YC
        REAL XP(NCMAX),YP(NCMAX)
        REAL XF(NCMAX),YF(NCMAX)
        CHARACTER*1 CH
        CHARACTER*50 CDUMMY
        CHARACTER*80 GLABEL
        LOGICAL LFIT(NCMAX)
        LOGICAL LCOLOR(MAX_ID_RED)
        COMMON/BLKNCBUFF/NCBUFF
        COMMON/BLKPLOT/XP,YP
        COMMON/BLKFITG1/NP
        COMMON/BLKFITG2/XF,YF
        COMMON/BLKRANGE/NN0,NN1,NN2
        COMMON/BLKCOLOR/MICOLOR,MISLS,MISLW,GLABEL
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        IF(CMODE.EQ.'?')THEN
          WRITE(*,100)'Select region by keyboard or mourse (k/m) '
          CMODE(1:1)=READC('m','km')
        END IF
        DO I=1,NCHAN
          LFIT(I)=.FALSE.
        END DO
        WRITE(CDUMMY,'(I10,A1,I10)')NN1,',',NN2
        CALL RMBLANK(CDUMMY,CDUMMY,L)
C
        IF(CMODE.EQ.'m') GOTO 20
C
10      WRITE(*,101)'Valid range: '//CDUMMY(1:L)
11      WRITE(*,100)'Enter channel region to be employed (0,0=EXIT) '
        CALL READ2I('0,0',NC1,NC2)
        IF((NC1.EQ.0).AND.(NC2.EQ.0))GOTO 60
        IF((NC1.LT.NN1).OR.(NC2.GT.NN2).OR.(NC1.GT.NC2))THEN
          WRITE(*,101)'Invalid numbers. Try again.'
          GOTO 10
        END IF
        NP=0
        DO I=NC1,NC2
          LFIT(I)=.TRUE.
          NP=NP+1
          XF(NP)=XP(I-NN1+1)
          YF(NP)=YP(I-NN1+1)
        END DO
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM)) CALL PGSCI(MICOLOR(NCBUFF)+1)
          CALL PGBIN(NP,XF,YF,.TRUE.)
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
        END DO
        GOTO 11
C
20      WRITE(*,101)'Enter X to exit from selection mode'
        WRITE(*,100)'Press mouse button (point #1)...'
        IF(LCOLOR(1)) CALL PGSCI(5)
        CALL RPGBAND(6,0,0.,0.,XC,YC,CH)
        IF(LCOLOR(1)) CALL PGSCI(1)
        IF((CH.NE.'A').AND.(CH.NE.'D').AND.(CH.NE.'X'))THEN
          WRITE(*,101)'ERROR: mouse button has not been detected.'
          GOTO 20
        END IF
        IF(CH.EQ.'X') GOTO 60
        NC1=INT(XC+0.5)
        IF(NC1.LT.NN1) NC1=NN1
        IF(NC1.GT.NN2) NC1=NN2
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM)) CALL PGSCI(MICOLOR(NCBUFF)+1)
          CALL PGMOVE(XP(NC1-NN1+1)-.5,YP(NC1-NN1+1))
          CALL PGDRAW(XP(NC1-NN1+1)+.5,YP(NC1-NN1+1))
        END DO
        WRITE(*,110)'  X-position: ',NC1
22      WRITE(*,100)'Press mouse button (point #2)...'
        IF(LCOLOR(1)) CALL PGSCI(5)
        CALL RPGBAND(4,0,REAL(NC1),0.,XC,YC,CH)
        IF(LCOLOR(1)) CALL PGSCI(1)
        IF((CH.NE.'A').AND.(CH.NE.'D').AND.(CH.NE.'X'))THEN
          WRITE(*,101)'ERROR: mouse button has not been detected.'
          GOTO 22
        END IF
        IF(CH.EQ.'X') GOTO 60
        NC2=INT(XC+0.5)
        IF(NC2.LT.NN1) NC2=NN1
        IF(NC2.GT.NN2) NC2=NN2
        WRITE(*,110)'  X-position: ',NC2
        IF(NC1.GT.NC2)THEN
          NC0=NC1
          NC1=NC2
          NC2=NC0
        END IF
        NP=0
        DO I=NC1,NC2
          LFIT(I)=.TRUE.
          NP=NP+1
          XF(NP)=XP(I-NN1+1)
          YF(NP)=YP(I-NN1+1)
        END DO
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM)) CALL PGSCI(MICOLOR(NCBUFF)+1)
          CALL PGBIN(NP,XF,YF,.TRUE.)
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
        END DO
        GOTO 20
C
60      NP=0
        DO I=NN1,NN2
          IF(LFIT(I))THEN
            NP=NP+1
            XF(NP)=XP(I-NN1+1)
            YF(NP)=YP(I-NN1+1)
          END IF
        END DO
100     FORMAT(A,$)
101     FORMAT(A)
110     FORMAT(A,I6)
        END
C
C******************************************************************************
C
        LOGICAL FUNCTION INSIDE(XC,YC)
        IMPLICIT NONE
        REAL XC,YC
        REAL XMIN,XMAX,YMIN,YMAX
        COMMON/BLKLIMITS/XMIN,XMAX,YMIN,YMAX
C
        INSIDE=.TRUE.
        IF(XC.GT.XMAX) INSIDE=.FALSE.
        IF(XC.LT.XMIN) INSIDE=.FALSE.
        IF(YC.GT.YMAX) INSIDE=.FALSE.
        IF(YC.LT.YMIN) INSIDE=.FALSE.
        END
C
C******************************************************************************
C******************************************************************************
C
C Permite dibujar la imagen completa o hacer zoom sobre ella
C Los parametros que se pasan a traves del COMMON son
C A(j,i,n) - la matriz imagen con i:scans, j:canales, n:no. buffer
C NSCAN:     el numero de scans de la imagen
C NCHAN:     el numero de canales de la imagen
        SUBROUTINE SUBLOOK(STWV,DISP)
        IMPLICIT NONE
C
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER TRUELEN
        INTEGER READILIM
        REAL READF
C
        INTEGER NBOTONES
        PARAMETER (NBOTONES=24)
        INTEGER NMAXCONTOUR
        PARAMETER (NMAXCONTOUR=100)
C
        INTEGER NCBUFF,NBUFF,NUSEDBUFF
        INTEGER NC1,NC2,NS1,NS2
        INTEGER NX1,NX2,NY1,NY2,NX0,NY0,NCOLOR
        INTEGER I,J,LFG1,LFG2,LBG1,LBG2
        INTEGER II,JJ
        INTEGER IXC1,IXC2,IYC1,IYC2
        INTEGER NB,NBLOCAL
        INTEGER BMODE(NBOTONES)
        INTEGER JUST
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        INTEGER NCONTOUR
        INTEGER NSCAN2,NCHAN2
        INTEGER NS1_,NS2_
        REAL AGRAY(NCMAX,NSMAX)
        REAL A(NCMAX,NSMAX,NMAXBUFF)
        REAL ARRAY3D(NCMAX,NSMAX),ANGLE3D
        REAL AMEAN,ASIGMA
        REAL BG,FG
        REAL TR(6),XC,YC
        REAL XMIN,XMAX,YMIN,YMAX,XMINL,XMAXL
        REAL CONTOUR(NMAXCONTOUR),LEVELMIN,LEVELMAX
        REAL STWV2,DISP2
        CHARACTER*1 CSEP_
        CHARACTER*1 CH,CCONTOUR
        CHARACTER*50 CFG1,CFG2,CBG1,CBG2
        CHARACTER*20 LABEL(NBOTONES)
        CHARACTER*50 CDUMMY
        CHARACTER*75 INFILE,INFILEBUFF(NMAXBUFF)
        LOGICAL LDEFBUFF(NMAXBUFF),LUSEBUFF(NMAXBUFF)
        LOGICAL INSIDE1
        LOGICAL LEXIT
        LOGICAL LW
        LOGICAL LBEXIST
        LOGICAL LCOLOR(MAX_ID_RED)
        LOGICAL LOK
        LOGICAL LOOP
C
        COMMON/BLKNCBUFF/NCBUFF
        COMMON/BLKDATA1/NSCAN,NCHAN
        COMMON/BLKDATA2/A
        COMMON/BLKBUFF1/LDEFBUFF,LUSEBUFF
        COMMON/BLKBUFF2/INFILEBUFF
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
        COMMON/BLKAGRAY/AGRAY
        COMMON/BLKANGLE3D/ANGLE3D
C
C------------------------------------------------------------------------------
        DATA (LABEL(I),I=1,NBOTONES)/
     +   '[m]enu','[z]oom (m)','zoom [k]','[w]hole',
     +   '[s]et BG/FG','[n]ew File',
     +   '3[d] plot','s[c]ale(x=\\b/y)','[b]uffer',
     +   '[+]>>>','<<<[-]','min[,]max',
     +   'x single','y single','x added','y added','x all','y all',
     +   '#[1]','#[2]','#[3]','#[4]','#[5]','#[6]'/
        DATA (BMODE(I),I=1,NBOTONES)/0, 0, 0, 0, 0, 0,
     +                               0, 1, 0, 0, 0, 0,
     +                               0, 0, 0, 0, 0, 0,
     +                               3, 3, 3, 3, 3, 3/
        OUTFILEX=OUTFILEX
C
        LW=((STWV.NE.0.0).AND.(DISP.NE.0.0))
C actualizamos buffers segun esten definidos o no y en utilizacion o no
        DO NBUFF=1,NMAXBUFF
          IF(LDEFBUFF(NBUFF))THEN
            IF(LUSEBUFF(NBUFF))THEN
              IF(NBUFF.EQ.NCBUFF)THEN
                BMODE(NBUFF+18)=5              !boton apretado, texto en blanco
              ELSE
                BMODE(NBUFF+18)=1                              !texto en blanco
              END IF
            ELSE
              BMODE(NBUFF+18)=0                                 !texto en negro
            END IF
          END IF
        END DO
C
        JUST=0
C
ccc5       DO I=1,NBOTONES
        DO I=1,NBOTONES
          IF(BMODE(I).GE.0) CALL BUTTON(I,LABEL(I),0)
          IF(BMODE(I).GT.0) CALL BUTTON(I,LABEL(I),BMODE(I))
        END DO
C
        TR(1)=0.
        TR(2)=1.
        TR(3)=0.
        TR(4)=0.
        TR(5)=0.
        TR(6)=1.
C
ccc10      NC1=1
        NC1=1
        NC2=NCHAN
        NS1=1
        NS2=NSCAN
        CALL SHOWLIMITS(NS1,NS2,NC1,NC2,AMEAN,ASIGMA)     !subrutina en plots.f
        DO I=NS1,NS2
          DO J=NC1,NC2
            AGRAY(J,I)=A(J,I,NCBUFF)
          END DO
        END DO
ccc     BG=AGRAY(NC1,NS1)
ccc        FG=BG
ccc        DO II=NS1,NS2
ccc          DO JJ=NC1,NC2
ccc            BG=AMIN1(BG,AGRAY(JJ,II))
ccc            FG=AMAX1(FG,AGRAY(JJ,II))
ccc          END DO
ccc        END DO
        BG=AMEAN-5.*ASIGMA
        FG=AMEAN+5.*ASIGMA
        WRITE(*,100)'Background: '
        WRITE(*,*)BG
        WRITE(*,100)'Foreground: '
        WRITE(*,*)FG
        WRITE(*,101)'----------------------------------'
        WRITE(*,101)'* Remember:'
        WRITE(*,101)' <*> plot contour map'
        WRITE(*,100)' <D> 3D-plot of current image area'
        WRITE(*,101)' (with hidden-line suppression)'
        WRITE(*,101)'----------------------------------'
C
16      DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(ITERM.EQ.1)THEN
            CALL RPGERASW(0.,1.,0.,0.80)
          ELSE
          END IF
        END DO
        XMIN=REAL(NC1)-0.8
        XMAX=REAL(NC2)+0.8
        YMIN=REAL(NS1)-0.8
        YMAX=REAL(NS2)+0.8
        IF(LW)THEN
          XMINL=(XMIN-1.)*DISP+STWV
          XMAXL=(XMAX-1.)*DISP+STWV
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(ITERM.EQ.1)THEN
              CALL RPGENV(XMIN,XMAX,YMIN,YMAX,JUST,-2)
            ELSE
              CALL PGENV(XMIN,XMAX,YMIN,YMAX,JUST,-2)
            END IF
            CALL PGSLW(3)
            CALL PGSWIN(XMINL,XMAXL,YMIN,YMAX)
            CALL PGBOX('CIMTS',0.,0,'BICNST',0.,0)
            CALL PGSWIN(XMIN,XMAX,YMIN,YMAX)
            CALL PGBOX('BINTS',0.,0,'IMTS',0.,0)
            CALL PGSLW(1)
          END DO
        ELSE
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGSLW(3)
            IF(ITERM.EQ.1)THEN
              CALL RPGENV(XMIN,XMAX,YMIN,YMAX,JUST,0)
            ELSE
              CALL PGENV(XMIN,XMAX,YMIN,YMAX,JUST,0)
            END IF
            CALL PGBOX('ITSM',0.,0,'ITSM',0.,0)
            CALL PGSLW(1)
          END DO
        END IF
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGIDEN_RED
          CALL PGSLW(3)
          CALL PGLABEL('channel','scan',CHAR(32))
          CALL PGMTEXT('T',2.5,0.5,0.5,INFILEBUFF(NCBUFF))
          CALL PGSLW(1)
        END DO
        WRITE(CFG1,*)FG
        CALL RMBLANK(CFG1,CFG1,LFG1)
        WRITE(CBG1,*)BG
        CALL RMBLANK(CBG1,CBG1,LBG1)
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGMTEXT('T',2.5,1.0,1.0,'FG: '//CFG1(1:LFG1))
          CALL PGMTEXT('T',2.5,0.0,0.0,'BG: '//CBG1(1:LBG1))
        END DO
        CALL SHOWLIMITS(NS1,NS2,NC1,NC2,AMEAN,ASIGMA)     !subrutina en plots.f
17      DO I=NS1,NS2
          DO J=NC1,NC2
            AGRAY(J,I)=A(J,I,NCBUFF)
          END DO
        END DO
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGGRAY(AGRAY,NCMAX,NSMAX,NC1,NC2,NS1,NS2,FG,BG,TR)
        END DO
        NCOLOR=2
C------------------------------------------------------------------------------
C
18      CONTINUE
        NB=0
C------------------------------------------------------------------------------
20      CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
        CALL IFBUTTON(XC,YC,NB)
        NBLOCAL=INDEX('mzkwsndcb+-,      123456',CH)
        IF((NBLOCAL.NE.0).AND.(CH.NE.' '))THEN
          CALL BUTTQEX(NBLOCAL,LBEXIST)
          IF(LBEXIST) NB=NBLOCAL
        END IF
C------------------------------------------------------------------------------
C opciones no establecidas con botones
C------------------------------------------------------------------------------
        IF(CH.EQ.'*')THEN
          WRITE(*,100)'Minimum level: '
          LEVELMIN=READF('@')
          WRITE(*,100)'Maximum level: '
          LEVELMAX=READF('@')
          WRITE(*,100)'No. of contour levels '
          NCONTOUR=READILIM('2',2,NMAXCONTOUR)
          DO I=1,NCONTOUR
            CONTOUR(I)=LEVELMIN+REAL(I-1)/REAL(NCONTOUR-1)*
     +       (LEVELMAX-LEVELMIN)
          END DO
          WRITE(*,101)'Contour subroutine:'
          WRITE(*,101)'(1) PGCONS'
          WRITE(*,101)'(2) PGCONT'
          WRITE(*,100)'Option (1/2) '
          CCONTOUR(1:1)=READC('1','12')
          NCOLOR=NCOLOR+1
          IF(NCOLOR.GT.12) NCOLOR=2
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(NCOLOR)
            IF(CCONTOUR.EQ.'1')THEN
              CALL PGCONS(AGRAY,NCMAX,NSMAX,NC1,NC2,NS1,NS2,
     +         CONTOUR,NCONTOUR,TR)
            ELSE
              CALL PGCONT(AGRAY,NCMAX,NSMAX,NC1,NC2,NS1,NS2,
     +         CONTOUR,NCONTOUR,TR)
            END IF
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END DO
          GOTO 20
C------------------------------------------------------------------------------
        ELSEIF(CH.EQ.'D')THEN
          WRITE(*,100)'Angle '
          WRITE(CDUMMY,*)ANGLE3D
          ANGLE3D=READF(CDUMMY)
          CALL PLOT3DH(ARRAY3D,NC1,NC2,NC2-NC1+1,NS1,NS2,NS2-NS1+1,
     +     BG,FG,ANGLE3D)
          WRITE(*,100)'Press any button to continue...'
          CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
          WRITE(*,101)'   ...OK!'
          GOTO 16
        END IF
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
        IF(NB.EQ.1)THEN
          CALL BUTTON(NB,LABEL(NB),5)
          RETURN
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.2)THEN
          CALL BUTTON(NB,LABEL(NB),5)
          WRITE(*,101)'Press cursor at two corners of the imaginary '
     +     //'BOX to be zoomed'
          CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
          IF((CH.NE.'A').AND.(CH.NE.'D').AND.(CH.NE.'X')) THEN
            WRITE(*,101)'ERROR: mouse buttom has not been detected.'
            CALL BUTTON(NB,LABEL(NB),0)
            GOTO 18
          END IF
          IXC1=INT(XC+0.5)
          IF(IXC1.LT.NC1) IXC1=NC1
          IF(IXC1.GT.NC2) IXC1=NC2
          IYC1=INT(YC+0.5)
          IF(IYC1.LT.NS1) IYC1=NS1
          IF(IYC1.GT.NS2) IYC1=NS2
          WRITE(*,112)'Cursor at ',IXC1,IYC1
C
          IF(LCOLOR(1)) CALL PGSCI(5)
          CALL RPGBAND(2,0,REAL(IXC1),REAL(IYC1),XC,YC,CH)
          IF(LCOLOR(1)) CALL PGSCI(1)
          IF((CH.NE.'A').AND.(CH.NE.'D').AND.(CH.NE.'X')) THEN
            WRITE(*,101)'ERROR: mouse buttom has not been detected.'
            CALL BUTTON(NB,LABEL(NB),0)
            GOTO 18
          END IF
          IXC2=INT(XC+0.5)
          IF(IXC2.LT.NC1) IXC2=NC1
          IF(IXC2.GT.NC2) IXC2=NC2
          IYC2=INT(YC+0.5)
          IF(IYC2.LT.NS1) IYC2=NS1
          IF(IYC2.GT.NS2) IYC2=NS2
          WRITE(*,112)'Cursor at ',IXC2,IYC2
C
          IF((IXC1.EQ.IXC2).AND.(IYC1.EQ.IYC2))THEN
            WRITE(*,101)'ERROR: invalid corners.'
            CALL BUTTON(NB,LABEL(NB),0)
            GOTO 18
          END IF
C
          NC1=MIN0(IXC1,IXC2)
          NC2=MAX0(IXC1,IXC2)
          NS1=MIN0(IYC1,IYC2)
          NS2=MAX0(IYC1,IYC2)
          CALL BUTTON(NB,LABEL(NB),0)
          GOTO 16
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.3)THEN
          CALL BUTTON(NB,LABEL(NB),5)
          WRITE(*,101)'Enter the coordinates of two corners of the '//
     +     'imaginary box to be zoomed.'
          WRITE(*,100)'First point  (channel,scan)'
          CALL READ2I('@',IXC1,IYC1)
          IF(.NOT.INSIDE1(IXC1,1,NCHAN,IYC1,1,NSCAN))THEN
            WRITE(*,101)'ERROR: limits out of plot.'
            CALL BUTTON(NB,LABEL(NB),0)
            GOTO 18
          END IF
          WRITE(*,100)'Second point (channel,scan)'
          CALL READ2I('@',IXC2,IYC2)
          IF(.NOT.INSIDE1(IXC2,1,NCHAN,IYC2,1,NSCAN))THEN
            WRITE(*,101)'ERROR: limits out of plot.'
            CALL BUTTON(NB,LABEL(NB),0)
            GOTO 18
          END IF
          NC1=MIN0(IXC1,IXC2)
          NC2=MAX0(IXC1,IXC2)
          NS1=MIN0(IYC1,IYC2)
          NS2=MAX0(IYC1,IYC2)
          CALL BUTTON(NB,LABEL(NB),0)
          GOTO 16
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.4)THEN
          CALL BUTTON(NB,LABEL(NB),5)
          IF(NCBUFF.EQ.0)THEN
            WRITE(*,101)'ERROR: no image selected!'
            CALL BUTTON(NB,LABEL(NB),0)
            GOTO 18
          END IF
          NC1=1
          NC2=NCHAN
          NS1=1
          NS2=NSCAN
          CALL BUTTON(NB,LABEL(NB),0)
          GOTO 16
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.5)THEN
          CALL BUTTON(NB,LABEL(NB),5)
          WRITE(*,100)'BACKGROUND    : '
          WRITE(*,*)BG
          WRITE(*,100)'FOREGROUND    : '
          WRITE(*,*)FG
          WRITE(CFG1,*)FG
          CALL RMBLANK(CFG1,CFG1,LFG1)
          WRITE(CBG1,*)BG
          CALL RMBLANK(CBG1,CBG1,LBG1)
          WRITE(*,100)'NEW BACKGROUND '
          WRITE(CDUMMY,*)BG
          BG=READF(CDUMMY)
          WRITE(*,100)'NEW FOREGROUND '
          WRITE(CDUMMY,*)FG
          FG=READF(CDUMMY)
          WRITE(CFG2,*)FG
          CALL RMBLANK(CFG2,CFG2,LFG2)
          WRITE(CBG2,*)BG
          CALL RMBLANK(CBG2,CBG2,LBG2)
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGSCI(0)
            CALL PGMTEXT('T',2.5,1.0,1.0,'FG: '//CFG1(1:LFG1))
            CALL PGMTEXT('T',2.5,0.0,0.0,'BG: '//CBG1(1:LBG1))
            CALL PGSCI(1)
            CALL PGMTEXT('T',2.5,1.0,1.0,'FG: '//CFG2(1:LFG2))
            CALL PGMTEXT('T',2.5,0.0,0.0,'BG: '//CBG2(1:LBG2))
          END DO
          CALL BUTTON(NB,LABEL(NB),0)
          GOTO 17
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.6)THEN
          CALL BUTTON(NB,LABEL(NB),5)
          LOOP=.TRUE.
          DO WHILE(LOOP)
            WRITE(*,100)'Input file name'
            INFILE=INFILEX(20,'@',NSCAN2,NCHAN2,STWV2,DISP2,1,.FALSE.)
            CALL HAYCOMA(INFILE,NSCAN2,NS1_,NS2_,CSEP_)
            LOOP=(NS1_.EQ.0)
            IF(LOOP) CLOSE(20)
          END DO
          IF(CSEP_.EQ.'+')THEN
            NSCAN2=1
          ELSE
            NSCAN2=NS2_-NS1_+1
          END IF
          CALL CCSIZE(NSCAN,NCHAN,STWV,DISP,
     +     NSCAN2,NCHAN2,STWV2,DISP2,LOK)
          IF(.NOT.LOK)THEN
            CALL BUTTON(NB,LABEL(NB),0)
            GOTO 20
          END IF
          BMODE(18+NCBUFF)=1                         !ya no es imagen principal
          CALL BUTTON(18+NCBUFF,LABEL(18+NCBUFF),0)
          CALL BUTTON(18+NCBUFF,LABEL(18+NCBUFF),BMODE(18+NCBUFF))
          CALL NEWBUFF(NCBUFF)
          WRITE(CDUMMY,*)NCBUFF
          WRITE(*,100)'Buffer numer '
          NCBUFF=READILIM(CDUMMY,1,NMAXBUFF)
          LDEFBUFF(NCBUFF)=.TRUE.
          IF(TRUELEN(OBJECT).GT.0)THEN
            INFILEBUFF(NCBUFF)=INFILE(1:TRUELEN(INFILE))//' ['//
     +       OBJECT(1:TRUELEN(OBJECT))//']'
          ELSE
            INFILEBUFF(NCBUFF)=INFILE(1:TRUELEN(INFILE))
          END IF
          LUSEBUFF(NCBUFF)=.TRUE.
          IF(NS1_.GT.1)THEN !ignoramos primeros NS1_-1 scans
            DO I=1,NS1_-1
              READ(20) (A(J,1,NCBUFF),J=1,NCHAN)
            END DO
          END IF
          DO I=NS1_,NS2_
            READ(20) (A(J,I-NS1_+1,NCBUFF),J=1,NCHAN)
          END DO
          CLOSE(20)
          IF(CSEP_.EQ.'+')THEN
            IF(NS2_.GT.NS1_)THEN
              DO J=1,NCHAN
                DO I=NS1_+1,NS2_
                A(J,1,NCBUFF)=A(J,1,NCBUFF)+A(J,I-NS1_+1,NCBUFF)
                END DO
                A(J,1,NCBUFF)=A(J,1,NCBUFF)/REAL(NS2_-NS1_+1)
              END DO
            END IF
          END IF
          BMODE(NCBUFF+18)=5                   !boton apretado, texto en blanco
          CALL BUTTON(NCBUFF+18,LABEL(NCBUFF+18),0)
          CALL BUTTON(NCBUFF+18,LABEL(NCBUFF+18),BMODE(NCBUFF+18))
          CALL SHOWLIMITS(NS1,NS2,NC1,NC2,AMEAN,ASIGMA)   !subrutina en plots.f
          CALL BUTTON(NB,LABEL(NB),0)
          GOTO 17
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.7)THEN
          CALL BUTTON(NB,LABEL(NB),5)
          CALL PLOT3D(NC1,NC2,NS1,NS2,BG,FG)
          WRITE(*,100)'Press any button to continue...'
          CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
          WRITE(*,101)'   ...OK!'
          CALL BUTTON(NB,LABEL(NB),0)
          GOTO 16
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.8)THEN
          IF(JUST.EQ.0)THEN
            JUST=1
            CALL BUTTON(8,'S[c]ale(x=y)',0)
            CALL BUTTON(8,'S[c]ale(x=y)',1)
          ELSE
            JUST=0
            CALL BUTTON(8,'  S[c]ale(x=\\b/y)',0)
            CALL BUTTON(8,'  S[c]ale(x=\\b/y)',1)
          END IF
          GOTO 16
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.9)THEN
          CALL BUTTON(NB,LABEL(NB),5)
          CALL SHOWBUFF
          CALL BUTTON(NB,LABEL(NB),0)
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.10)THEN
          IF(NCBUFF.EQ.0)THEN
            WRITE(*,101)'ERROR: select buffer number first!'
            GOTO 20
          END IF
          CALL BUTTON(NB,LABEL(NB),5)
          NUSEDBUFF=0                            !numero de buffers utilizables
          DO NBUFF=1,NMAXBUFF
            IF(LUSEBUFF(NBUFF)) NUSEDBUFF=NUSEDBUFF+1
          END DO
          IF(NUSEDBUFF.GT.1)THEN                  !solo si hay mas de un buffer
            CALL BUTTON(NCBUFF+18,LABEL(NCBUFF+18),0)
            CALL BUTTON(NCBUFF+18,LABEL(NCBUFF+18),1)
            LEXIT=.FALSE.
            DO WHILE(.NOT.LEXIT)
              NCBUFF=NCBUFF+1
              IF(NCBUFF.GT.NMAXBUFF) NCBUFF=1
              IF(LUSEBUFF(NCBUFF))THEN
                LEXIT=.TRUE.
                CALL BUTTON(NCBUFF+18,LABEL(NCBUFF+18),5)
              END IF
            END DO
          ELSE
            CALL BUTTON(NB,LABEL(NB),0)
            GOTO 20
          END IF
          CALL BUTTON(NB,LABEL(NB),0)
          GOTO 17
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.11)THEN
          CALL BUTTON(NB,LABEL(NB),5)
          IF(NCBUFF.EQ.0)THEN
            WRITE(*,101)'ERROR: select buffer number first!'
            CALL BUTTON(NB,LABEL(NB),0)
            GOTO 20
          END IF
          NUSEDBUFF=0                            !numero de buffers utilizables
          DO NBUFF=1,NMAXBUFF
            IF(LUSEBUFF(NBUFF)) NUSEDBUFF=NUSEDBUFF+1
          END DO
          IF(NUSEDBUFF.GT.1)THEN                  !solo si hay mas de un buffer
            CALL BUTTON(NCBUFF+18,LABEL(NCBUFF+18),0)
            CALL BUTTON(NCBUFF+18,LABEL(NCBUFF+18),1)
            LEXIT=.FALSE.
            DO WHILE(.NOT.LEXIT)
              NCBUFF=NCBUFF-1
              IF(NCBUFF.LT.1) NCBUFF=NMAXBUFF
              IF(LUSEBUFF(NCBUFF))THEN
                LEXIT=.TRUE.
                CALL BUTTON(NCBUFF+18,LABEL(NCBUFF+18),5)
              END IF
            END DO
          ELSE
            CALL BUTTON(NB,LABEL(NB),0)
            GOTO 20
          END IF
          CALL BUTTON(NB,LABEL(NB),0)
          GOTO 17
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.12)THEN
          CALL BUTTON(NB,LABEL(NB),5)
          IF(NCBUFF.EQ.0)THEN
            WRITE(*,101)'ERROR: no image selected!'
            CALL BUTTON(NB,LABEL(NB),0)
            GOTO 18
          END IF
          WRITE(CFG1,*)FG
          CALL RMBLANK(CFG1,CFG1,LFG1)
          WRITE(CBG1,*)BG
          CALL RMBLANK(CBG1,CBG1,LBG1)
          BG=AGRAY(NC1,NS1)
          FG=BG
          DO II=NS1,NS2
            DO JJ=NC1,NC2
              BG=AMIN1(BG,AGRAY(JJ,II))
              FG=AMAX1(FG,AGRAY(JJ,II))
            END DO
          END DO
          WRITE(*,100)'Background: '
          WRITE(*,*)BG
          WRITE(*,100)'Foreground: '
          WRITE(*,*)FG
          WRITE(CFG2,*)FG
          CALL RMBLANK(CFG2,CFG2,LFG2)
          WRITE(CBG2,*)BG
          CALL RMBLANK(CBG2,CBG2,LBG2)
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGSCI(0)
            CALL PGMTEXT('T',2.5,1.0,1.0,'FG: '//CFG1(1:LFG1))
            CALL PGMTEXT('T',2.5,0.0,0.0,'BG: '//CBG1(1:LBG1))
            CALL PGSCI(1)
            CALL PGMTEXT('T',2.5,1.0,1.0,'FG: '//CFG2(1:LFG2))
            CALL PGMTEXT('T',2.5,0.0,0.0,'BG: '//CBG2(1:LBG2))
          END DO
          CALL BUTTON(NB,LABEL(NB),0)
          GOTO 17
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.13)THEN                                     !x cut, single
          CALL BUTTON(NB,LABEL(NB),5)
          IF(NCBUFF.EQ.0)THEN
            WRITE(*,101)'ERROR: no image selected!'
            CALL BUTTON(NB,LABEL(NB),0)
            GOTO 18
          END IF
          WRITE(*,100)'Press mouse button...'
          IF(LCOLOR(1)) CALL PGSCI(5)
          CALL RPGBAND(5,0,0.,0.,XC,YC,CH)
          IF(LCOLOR(1)) CALL PGSCI(1)
          WRITE(*,100)'  ..thanks'
          IF((CH.NE.'A').AND.(CH.NE.'D').AND.(CH.NE.'X')) THEN
            WRITE(*,*)
            WRITE(*,101)'ERROR: mouse buttom has not been detected.'
            CALL BUTTON(NB,LABEL(NB),0)
            GOTO 18
          END IF
          NY1=INT(YC+0.5)
          IF(NY1.LT.NS1) NY1=NS1
          IF(NY1.GT.NS2) NY1=NS2
          WRITE(*,110)'  -> Scan #',NY1
          NY2=NY1
          NCOLOR=NCOLOR+1
          IF(NCOLOR.GT.12) NCOLOR=2
          CALL PLOTLOCAL('x',NY1,NY2,XMIN,XMAX,YMIN,YMAX,NCOLOR,FG,BG)
          CALL BUTTON(NB,LABEL(NB),0)
          GOTO 18
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.14)THEN                                     !y cut, single
          CALL BUTTON(NB,LABEL(NB),5)
          IF(NCBUFF.EQ.0)THEN
            WRITE(*,101)'ERROR: no image selected!'
            CALL BUTTON(NB,LABEL(NB),0)
            GOTO 18
          END IF
          WRITE(*,100)'Press mouse button...'
          IF(LCOLOR(1)) CALL PGSCI(5)
          CALL RPGBAND(6,0,0.,0.,XC,YC,CH)
          IF(LCOLOR(1)) CALL PGSCI(1)
          WRITE(*,100)'  ..thanks'
          IF((CH.NE.'A').AND.(CH.NE.'D').AND.(CH.NE.'X')) THEN
            WRITE(*,*)
            WRITE(*,101)'ERROR: mouse buttom has not been detected.'
            CALL BUTTON(NB,LABEL(NB),0)
            GOTO 18
          END IF
          NX1=INT(XC+0.5)
          IF(NX1.LT.NC1) NX1=NC1
          IF(NX1.GT.NC2) NX1=NC2
          WRITE(*,110)'  -> Channel #',NX1
          NX2=NX1
          NCOLOR=NCOLOR+1
          IF(NCOLOR.GT.12) NCOLOR=2
          CALL PLOTLOCAL('y',NX1,NX2,XMIN,XMAX,YMIN,YMAX,NCOLOR,FG,BG)
          CALL BUTTON(NB,LABEL(NB),0)
          GOTO 18
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.15)THEN                                      !x cut, added
          CALL BUTTON(NB,LABEL(NB),5)
          IF(NCBUFF.EQ.0)THEN
            WRITE(*,101)'ERROR: no image selected!'
            CALL BUTTON(NB,LABEL(NB),0)
            GOTO 18
          END IF
          WRITE(*,100)'Press mouse button...'
          IF(LCOLOR(1)) CALL PGSCI(5)
          CALL RPGBAND(5,0,0.,0.,XC,YC,CH)
          IF(LCOLOR(1)) CALL PGSCI(1)
          WRITE(*,100)'  ..thanks'
          IF((CH.NE.'A').AND.(CH.NE.'D').AND.(CH.NE.'X')) THEN
            WRITE(*,*)
            WRITE(*,101)'ERROR: mouse buttom has not been detected.'
            CALL BUTTON(NB,LABEL(NB),0)
            GOTO 18
          END IF
          NY1=INT(YC+0.5)
          IF(NY1.LT.NS1) NY1=NS1
          IF(NY1.GT.NS2) NY1=NS2
          WRITE(*,110)'  -> Scan #',NY1
          WRITE(*,100)'Press mouse button...'
          IF(LCOLOR(1)) CALL PGSCI(5)
          CALL RPGBAND(3,0,0.,REAL(NY1),XC,YC,CH)
          IF(LCOLOR(1)) CALL PGSCI(1)
          WRITE(*,100)'  ..thanks'
          IF((CH.NE.'A').AND.(CH.NE.'D').AND.(CH.NE.'X')) THEN
            WRITE(*,*)
            WRITE(*,101)'ERROR: mouse buttom has not been detected.'
            CALL BUTTON(NB,LABEL(NB),0)
            GOTO 18
          END IF
          NY2=INT(YC+0.5)
          IF(NY2.LT.NS1) NY2=NS1
          IF(NY2.GT.NS2) NY2=NS2
          WRITE(*,110)'  -> Scan #',NY2
          IF(NY1.GT.NY2)THEN
            NY0=NY1
            NY1=NY2
            NY2=NY1
          END IF
          NCOLOR=NCOLOR+1
          IF(NCOLOR.GT.12) NCOLOR=2
          CALL PLOTLOCAL('x',NY1,NY2,XMIN,XMAX,YMIN,YMAX,NCOLOR,FG,BG)
          CALL BUTTON(NB,LABEL(NB),0)
          GOTO 18
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.16)THEN                                      !y cut, added
          CALL BUTTON(NB,LABEL(NB),5)
          IF(NCBUFF.EQ.0)THEN
            WRITE(*,101)'ERROR: no image selected!'
            CALL BUTTON(NB,LABEL(NB),0)
            GOTO 18
          END IF
          WRITE(*,100)'Press mouse button...'
          IF(LCOLOR(1)) CALL PGSCI(5)
          CALL RPGBAND(6,0,0.,0.,XC,YC,CH)
          IF(LCOLOR(1)) CALL PGSCI(1)
          WRITE(*,100)'  ..thanks'
          IF((CH.NE.'A').AND.(CH.NE.'D').AND.(CH.NE.'X')) THEN
            WRITE(*,*)
            WRITE(*,101)'ERROR: mouse buttom has not been detected.'
            CALL BUTTON(NB,LABEL(NB),0)
            GOTO 18
          END IF
          NX1=INT(XC+0.5)
          IF(NX1.LT.NC1) NX1=NC1
          IF(NX1.GT.NC2) NX1=NC2
          WRITE(*,110)'  -> Channel #',NX1
          WRITE(*,100)'Press mouse button...'
          IF(LCOLOR(1)) CALL PGSCI(5)
          CALL RPGBAND(4,0,REAL(NX1),0.,XC,YC,CH)
          IF(LCOLOR(1)) CALL PGSCI(1)
          WRITE(*,100)'  ..thanks'
          IF((CH.NE.'A').AND.(CH.NE.'D').AND.(CH.NE.'X')) THEN
            WRITE(*,*)
            WRITE(*,101)'ERROR: mouse buttom has not been detected.'
            CALL BUTTON(NB,LABEL(NB),0)
            GOTO 18
          END IF
          NX2=INT(XC+0.5)
          IF(NX2.LT.NC1) NX2=NC1
          IF(NX2.GT.NC2) NX2=NC2
          WRITE(*,110)'  -> Channel #',NX2
          IF(NX1.GT.NX2)THEN
            NX0=NX1
            NX1=NX2
            NX2=NX1
          END IF
          NCOLOR=NCOLOR+1
          IF(NCOLOR.GT.12) NCOLOR=2
          CALL PLOTLOCAL('y',NX1,NX2,XMIN,XMAX,YMIN,YMAX,NCOLOR,FG,BG)
          CALL BUTTON(NB,LABEL(NB),0)
          GOTO 18
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.17)THEN                                        !x cut, all
          CALL BUTTON(NB,LABEL(NB),5)
          IF(NCBUFF.EQ.0)THEN
            WRITE(*,101)'ERROR: no image selected!'
            CALL BUTTON(NB,LABEL(NB),0)
            GOTO 18
          END IF
          NY1=NS1
          NY2=NS2
          WRITE(*,110)'  -> Scan #',NY1
          WRITE(*,110)'  -> Scan #',NY2
          NCOLOR=NCOLOR+1
          IF(NCOLOR.GT.12) NCOLOR=2
          CALL PLOTLOCAL('x',NY1,NY2,XMIN,XMAX,YMIN,YMAX,NCOLOR,FG,BG)
          CALL BUTTON(NB,LABEL(NB),0)
          GOTO 18
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.18)THEN                                        !y cut, all
          CALL BUTTON(NB,LABEL(NB),5)
          IF(NCBUFF.EQ.0)THEN
            WRITE(*,101)'ERROR: no image selected!'
            CALL BUTTON(NB,LABEL(NB),0)
            GOTO 18
          END IF
          NX1=NC1
          NX2=NC2
          WRITE(*,110)'  -> Scan #',NX1
          WRITE(*,110)'  -> Scan #',NX2
          NCOLOR=NCOLOR+1
          IF(NCOLOR.GT.12) NCOLOR=2
          CALL PLOTLOCAL('y',NX1,NX2,XMIN,XMAX,YMIN,YMAX,NCOLOR,FG,BG)
          CALL BUTTON(NB,LABEL(NB),0)
          GOTO 18
C------------------------------------------------------------------------------
        ELSEIF((NB.GE.19).AND.(NB.LE.24))THEN
          IF(NB-18.EQ.NCBUFF)THEN                    !si es NCBUFF se desactiva
            LUSEBUFF(NB-18)=.FALSE.
            CALL BUTTON(NB,LABEL(NB),0)
            NCBUFF=0
          ELSE                                                 !si no es NCBUFF
            IF(LUSEBUFF(NB-18))THEN        !si se usaba, se convierte en NCBUFF
              IF(NCBUFF.NE.0)THEN       !si existe, se quita el anterior NCBUFF
                CALL BUTTON(NCBUFF+18,LABEL(NCBUFF+18),0)
                CALL BUTTON(NCBUFF+18,LABEL(NCBUFF+18),1)
              END IF
              NCBUFF=NB-18
              CALL BUTTON(NB,LABEL(NB),5)
            ELSE                        !si no se usaba, se convierte en usable
              LUSEBUFF(NB-18)=.TRUE.
              CALL BUTTON(NB,LABEL(NB),0)
              CALL BUTTON(NB,LABEL(NB),1)
            END IF
          END IF
C------------------------------------------------------------------------------
        ELSE
          JJ=INT(XC+0.5)
          II=INT(YC+0.5)
          IF(INSIDE1(JJ,1,NCHAN,II,1,NSCAN))THEN
            WRITE(*,111)'Cursor at ',JJ,II,'      Pixel value: '
            WRITE(*,*)AGRAY(JJ,II)
          ELSE
            WRITE(*,101)'ERROR: Cursor out of plot.'
          END IF
        END IF
C------------------------------------------------------------------------------
        GOTO 20
C
100     FORMAT(A,$)
101     FORMAT(A)
110     FORMAT(A,I6)
112     FORMAT(A,I5,2X,I5)
111     FORMAT(A,I5,2X,I5,A,$)
C
        END
C
C******************************************************************************
C Hace un plot sobre la imagen de SUBLOOK
C
        SUBROUTINE PLOTLOCAL(CAXIS,NX1,NX2,XMIN,XMAX,YMIN,YMAX,
     +   NCOLOR,FG,BG)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
C
        CHARACTER*1 CAXIS
        INTEGER NX1,NX2
        REAL XMIN,XMAX,YMIN,YMAX
        INTEGER NCOLOR
        REAL FG,BG
C
        INTEGER I,J
        INTEGER NCBUFF
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        REAL XPLOT(NCMAX),YPLOT(NCMAX)               !el mayor de NCMAX y NSMAX
        LOGICAL LCOLOR(MAX_ID_RED)
C
        REAL A(NCMAX,NSMAX,NMAXBUFF)
        COMMON/BLKDATA1/NSCAN,NCHAN
        COMMON/BLKDATA2/A
        COMMON/BLKNCBUFF/NCBUFF
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        IF(CAXIS.EQ.'x')THEN
          DO J=1,NCHAN
            XPLOT(J)=REAL(J)
            YPLOT(J)=0.
          END DO
          DO I=NX1,NX2
            DO J=1,NCHAN
              YPLOT(J)=YPLOT(J)+A(J,I,NCBUFF)
            END DO 
          END DO
          DO J=1,NCHAN
            YPLOT(J)=YPLOT(J)/REAL(NX2-NX1+1)
          END DO
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGWINDOW(XMIN,XMAX,BG,FG)
            IF(LCOLOR(ITERM)) CALL PGSCI(NCOLOR)
            CALL PGBIN(NCHAN,XPLOT,YPLOT,.TRUE.)
            CALL PGWINDOW(XMIN,XMAX,YMIN,YMAX)
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END DO
C..............................................................................
        ELSE !CAXIS.EQ.'y'
          DO I=1,NSCAN
            XPLOT(I)=REAL(I)
            YPLOT(I)=0.
          END DO
          DO J=NX1,NX2
            DO I=1,NSCAN
              YPLOT(I)=YPLOT(I)+A(J,I,NCBUFF)
            END DO 
          END DO
          DO I=1,NSCAN
            YPLOT(I)=YPLOT(I)/REAL(NX2-NX1+1)
          END DO
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGWINDOW(YMIN,YMAX,BG,FG)
            IF(LCOLOR(ITERM)) CALL PGSCI(NCOLOR)
            CALL PGBIN(NSCAN,XPLOT,YPLOT,.TRUE.)
            CALL PGWINDOW(XMIN,XMAX,YMIN,YMAX)
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END DO
        END IF
C------------------------------------------------------------------------------
        END
C
C******************************************************************************
C
C Si el punto J0,I0 esta dentro del rectangulo definido por J1,J2,I1,I2
C la funcion devuelve .TRUE.
        LOGICAL FUNCTION INSIDE1(J0,J1,J2,I0,I1,I2)
        IMPLICIT NONE
        INTEGER J0,J1,J2
        INTEGER I0,I1,I2
C
        INSIDE1=.TRUE.
        IF(J0.LT.J1) GOTO 70
        IF(J0.GT.J2) GOTO 70
        IF(I0.LT.I1) GOTO 70
        IF(I0.GT.I2) GOTO 70
        RETURN
70      CONTINUE
        INSIDE1=.FALSE.
        RETURN
        END
C
C******************************************************************************
C Muestra el status actual de los diferentes buffers
C
        SUBROUTINE SHOWBUFF
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INTEGER TRUELEN
C
        INTEGER NCBUFF
        INTEGER NBUFF
        INTEGER L
        CHARACTER*75 INFILEBUFF(NMAXBUFF)
        LOGICAL LDEFBUFF(NMAXBUFF),LUSEBUFF(NMAXBUFF)
C
        COMMON/BLKDATA1/NSCAN,NCHAN
        COMMON/BLKNCBUFF/NCBUFF
        COMMON/BLKBUFF1/LDEFBUFF,LUSEBUFF
        COMMON/BLKBUFF2/INFILEBUFF
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        WRITE(*,*)
        WRITE(*,101)'>>> Buffer status: '
        WRITE(*,100)'    NSCAN: '
        WRITE(*,*)NSCAN
        WRITE(*,100)'    NCHAN: '
        WRITE(*,*)NCHAN
        IF(NCBUFF.EQ.0)THEN
          WRITE(*,101)'>>> WARNING: No current buffer defined!'
        ELSE
          WRITE(*,112)'>>> Current buffer is #',NCBUFF
        END IF
        DO NBUFF=1,NMAXBUFF
          WRITE(*,'(A,I1,$)')'#',NBUFF
          IF(LDEFBUFF(NBUFF))THEN
            L=TRUELEN(INFILEBUFF(NBUFF))
            IF(LUSEBUFF(NBUFF))THEN
              WRITE(*,100)':    [ACTIVE]: '//INFILEBUFF(NBUFF)(1:L)
            ELSE
              WRITE(*,100)': [NO ACTIVE]: '//INFILEBUFF(NBUFF)(1:L)
            END IF
            WRITE(*,*)
          ELSE
            WRITE(*,'(A)')': [undefined]'
          END IF
        END DO
        WRITE(*,*)
C
100     FORMAT(A,$)
101     FORMAT(A)
112     FORMAT(A,I1)
        END
C
C******************************************************************************
C Devuelve un numero de buffer todavia no definido. Si todos los buffers estan
C definidos, devuelve 1
C
        SUBROUTINE NEWBUFF(NBUFF)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INTEGER NBUFF
C
        LOGICAL LDEFBUFF(NMAXBUFF),LUSEBUFF(NMAXBUFF)
C
        COMMON/BLKBUFF1/LDEFBUFF,LUSEBUFF
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        DO NBUFF=1,NMAXBUFF
          IF(.NOT.LDEFBUFF(NBUFF)) RETURN
        END DO
        NBUFF=1
C
        END
C
C******************************************************************************
C Realiza una estadistica sencilla
C
        SUBROUTINE GIVESTAT(NPTOS,XP,YP)
        IMPLICIT NONE
C
        INTEGER NPTOS
        REAL XP(NPTOS),YP(NPTOS)
C
        INTEGER J,JMIN,JMAX,L
        REAL DATAMIN,DATAMAX
        DOUBLE PRECISION MEAN,SIGMA
        CHARACTER*50 CDUMMY
C------------------------------------------------------------------------------
        MEAN=0.D0
        DATAMIN=YP(1)
        DATAMAX=DATAMIN
        JMIN=1
        JMAX=1
        DO J=1,NPTOS
          MEAN=MEAN+DBLE(YP(J))
          IF(YP(J).LT.DATAMIN)THEN
            JMIN=J
            DATAMIN=YP(J)
          END IF
          IF(YP(J).GT.DATAMAX)THEN
            JMAX=J
            DATAMAX=YP(J)
          END IF
        END DO
        MEAN=MEAN/DBLE(NPTOS)
        SIGMA=0.D0
        IF(NPTOS.GT.1)THEN
          DO J=1,NPTOS
            SIGMA=SIGMA+(MEAN-DBLE(YP(J)))*(MEAN-DBLE(YP(J)))
          END DO
          SIGMA=DSQRT(SIGMA/DBLE(NPTOS-1))
        END IF
        WRITE(*,100)'> Total no. of points: '
        WRITE(*,*) NPTOS
        WRITE(CDUMMY,*)DATAMIN
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,'(A,A,A,$)')'> Minimum: ',CDUMMY(1:L),' pixel: '
        WRITE(*,*)XP(JMIN)
        WRITE(CDUMMY,*)DATAMAX
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,'(A,A,A,$)')'> Maximum: ',CDUMMY(1:L),' pixel: '
        WRITE(*,*)XP(JMAX)
        WRITE(*,100)'> Mean   : '
        WRITE(*,*)REAL(MEAN)
        WRITE(*,100)'> Sigma  : '
        WRITE(*,*)REAL(SIGMA)
C
100     FORMAT(A,$)
        END
C
C******************************************************************************
C Dibuja en proyeccion 3D el trozo de imagen NC1,NC2,NS1,NS2 tomando como
C valores minimo y maximo BG y FG respectivamente.
C
        SUBROUTINE PLOT3D(NC1,NC2,NS1,NS2,BG,FG)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INTEGER NC1,NC2,NS1,NS2
        REAL BG,FG
        REAL AGRAY(NCMAX,NSMAX)
C
        INTEGER I,J,II,JJ
        INTEGER NX,NY
        INTEGER NSMAX_LOCAL,NCMAX_LOCAL
        INTEGER NTERM,IDN(MAX_ID_RED)
        INTEGER ITERM
C XP,YP hay que dimensionarlos al mayor de NCMAX,NSMAX
        REAL XP(NCMAX),YP(NCMAX)
        REAL ZMINSS,ZMAXSS,OFFSS
        REAL X(NCMAX),Y(NSMAX)
        REAL FACTOR
C
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKAGRAY/AGRAY
C-------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        NSMAX_LOCAL=NSMAX
        NCMAX_LOCAL=NCMAX
        IF(NSMAX_LOCAL.GT.NCMAX_LOCAL)STOP'FATAL ERROR: NSMAX.GT.NCMAX'
C FACTOR determina como aparece de ampliada la estructura en el grafico
C que simula una representacion tridimensional de la superficie
        FACTOR=4.
        ZMINSS=BG
        ZMAXSS=FG
        OFFSS=ZMAXSS-ZMINSS
        ZMINSS=ZMINSS-OFFSS*.05
        DO J=NC1,NC2
          X(J-NC1+1)=REAL(J)
        END DO
        DO I=NS1,NS2
          Y(I-NS1+1)=REAL(I)
        END DO
        NX=NC2-NC1+1
        NY=NS2-NS1+1
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(ITERM.EQ.1)THEN
            CALL RPGERASW(0.,1.,0.,0.80)
            CALL RPGENV(X(1)-.05,X(NX)+(X(NX)-X(1))/3.+.05,
     +       ZMINSS,ZMAXSS+FACTOR*OFFSS,0,-2)
          ELSE
            CALL PGENV(X(1)-.05,X(NX)+(X(NX)-X(1))/3.+.05,
     +       ZMINSS,ZMAXSS+FACTOR*OFFSS,0,-2)
          END IF
        END DO
C
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
c.......
          DO I=1,NY
            II=I+NS1-1
            DO J=1,NX
              JJ=J+NC1-1
              IF(NY.EQ.1)THEN
                XP(J)=X(J)
                YP(J)=AGRAY(JJ,II)
              ELSE
                XP(J)=X(J)+(Y(I)-Y(1))/(Y(NY)-Y(1))*(X(NX)-X(1))/3.
                YP(J)=AGRAY(JJ,II)+(Y(I)-Y(1))/(Y(NY)-Y(1))*FACTOR*OFFSS
              END IF
            END DO
            CALL PGLINE(NX,XP,YP)
          END DO
          IF(NY.GT.1)THEN
            DO J=1,NX
              JJ=J+NC1-1
              DO I=1,NY
                II=I+NS1-1
                XP(I)=X(J)+(Y(I)-Y(1))/(Y(NY)-Y(1))*(X(NX)-X(1))/3.
                YP(I)=AGRAY(JJ,II)+(Y(I)-Y(1))/(Y(NY)-Y(1))*FACTOR*OFFSS
              END DO
              CALL PGLINE(NY,XP,YP)
            END DO
          END IF
c.......
        END DO
C
        END
C
C******************************************************************************
C Dibuja en proyeccion 3D el trozo de imagen NC1,NC2,NS1,NS2 tomando como
C valores minimo y maximo BG y FG respectivamente, y teniendo en cuenta las
C superficies ocultas (utiliza la subrutina FREDDY de PGPLOT).
C
        SUBROUTINE PLOT3DH(ARRAY3D,NC1,NC2,KX,NS1,NS2,NY,BG,FG,ANGLE)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INTEGER NC1,NC2,NS1,NS2,KX,NY
        REAL ARRAY3D(KX,NY)
        REAL BG,FG,ANGLE
C
        INTEGER I,J
        INTEGER ITERM,NTERM,IDN(MAX_ID_RED)
        REAL SIZE
        REAL AGRAY(NCMAX,NSMAX)
C
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKAGRAY/AGRAY
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
       SIZE=1.0
C
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          DO I=NS1,NS2
            DO J=NC1,NC2
              ARRAY3D(J-NC1+1,I-NS1+1)=AGRAY(J,I)
            END DO
          END DO
          IF(ITERM.EQ.1)THEN
            CALL RPGERASW(0.,1.,0.,0.80)
            CALL RPGENV(0.,SIZE,0.,SIZE,0,-2)
          ELSE
            CALL PGENV(0.,SIZE,0.,SIZE,0,-2)
          END IF
          CALL FREDDY(ARRAY3D,KX,NY,SIZE,ANGLE,BG,FG)
        END DO
C
        END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE FREDDY(ARRAY,KX,NY,SIZE,ANGLE,BG,FG)
      INTEGER KX, NY
      REAL ARRAY(KX,NY), SIZE, ANGLE
C
C Draws isometric plot of array
C
      REAL FMAX,FMIN,DELTAX,DELTAY,DELTAV,SINE,PEAK,X,DX,HEIGHT
      INTEGER I,J,KI,KJ,NX,MX,MY,STEP,LEFT,RIGHT,IT,MN,INCX
      LOGICAL VISBLE
      COMMON /FREDCM/ DELTAX,X,STEP,LEFT,RIGHT,IT,NX,VISBLE
C
      MN = KX*NY
      NX = KX
C     Check array size:
      IF(NX.LT.2 .OR. NY.LT.2) RETURN
        FMAX=FG
        FMIN=BG
ccc      FMAX = ARRAY(1,1)
ccc      FMIN = FMAX
ccc      DO 20 J=1,NY
ccc          DO 10 I=1,NX
ccc              FMIN = AMIN1(ARRAY(I,J),FMIN)
ccc              FMAX = AMAX1(ARRAY(I,J),FMAX)
ccc   10     CONTINUE
ccc   20 CONTINUE
      DELTAX = SIZE/(NX+NY)
      SINE = SIN(ANGLE/58.)
      DELTAY = DELTAX*SINE
      HEIGHT = SIZE*(1.-ABS(SINE))
      DELTAV = HEIGHT
      FMAX = FMAX-FMIN
      IF(FMAX.LT.0.0001) FMAX = DELTAV
      DELTAV = DELTAV/FMAX
      MX = NX+1
      MY = NY+1
      STEP = MX
C
C Start PGPLOT buffering.
C
      CALL PGBBUF
C
C Work our way down the Y axis, then up the X axis,
C calculating the Y plotter coordinates for each
C column of the plot, doing the hidden-line suppression
C at the same time.
C
      DO 50 J=1,NY
          KJ = MY-J
          KI = 1
C               ( KI,KJ are coordinates of bottom of column)
          ARRAY(KI,KJ) = DELTAY*(KI+KJ) + DELTAV*(ARRAY(KI,KJ)-FMIN)
   30     PEAK = ARRAY(KI,KJ)
   40     KI = KI+1
          KJ = KJ+1
          IF(KI.GT.NX .OR. KJ.GT.NY) GOTO 50
          ARRAY(KI,KJ) = DELTAY*(KI+KJ) + DELTAV*(ARRAY(KI,KJ)-FMIN)
          IF(ARRAY(KI,KJ).GT.PEAK) GOTO 30
          IF(ARRAY(KI,KJ).LE.PEAK) ARRAY(KI,KJ) = -ABS(ARRAY(KI,KJ))
          GOTO 40
   50 CONTINUE
C
C Now to work our way up the X axis
C
      DO 80 I=2,NX
          KI = I
          KJ = 1
          ARRAY(KI,KJ) = DELTAY*(KI+KJ)+DELTAV*(ARRAY(KI,KJ)-FMIN)
   60     PEAK = ARRAY(KI,KJ)
   70     KI = KI+1
          KJ = KJ+1
          IF(KI.GT.NX .OR. KJ.GT.NY) GOTO 80
          ARRAY(KI,KJ) = DELTAY*(KI+KJ)+DELTAV*(ARRAY(KI,KJ)-FMIN)
          IF(ARRAY(KI,KJ).GT.PEAK) GOTO 60
          IF(ARRAY(KI,KJ).LE.PEAK) ARRAY(KI,KJ) = -ABS(ARRAY(KI,KJ))
          GOTO 70
   80 CONTINUE
C
C Draw a line along the bottom of the vertical faces
C
      CALL PGMOVE(DELTAX*(NX+NY-2), DELTAY*(MX))
      CALL PGDRAW(DELTAX*(NY-1),    DELTAY*2)
      CALL PGDRAW(0.0,              DELTAY*MY)
C
C Array is now ready for plotting.  If a point is
C positive, then it is to be plotted at that Y
C coordinate; if it is negative, then it is
C invisible, but at minus that Y coordinate (the point
C where the line heading towards it disappears has to
C be determined by finding the intersection of it and
C the cresting line).
C
C Plot rows:
C
      DO 110 J=1,NY,2
          KJ = MY-J
          DX = DELTAX*(J-2)
          X = DX+DELTAX
          CALL PGMOVE(X,DELTAY*(KJ+1))
          CALL PGDRAW(X,ARRAY(1,KJ))
          VISBLE = .TRUE.
          DO 90 I=2,NX
              RIGHT = I+NX*(KJ-1)
              LEFT = RIGHT-1
              IT = RIGHT
              X = DX+DELTAX*I
              CALL FREDGO(ARRAY,MN)
   90     CONTINUE
C
C Now at far end of row so come back
C
          KJ = KJ-1
          IF(KJ.LE.0) GOTO 170
          VISBLE = ARRAY(NX,KJ).GE.0.0
          DX = DELTAX*(NX+J)
          IF(VISBLE) CALL PGMOVE(DX-DELTAX,ARRAY(NX,KJ))
          DELTAX = -DELTAX
          DO 100 I=2,NX
              KI = MX-I
              LEFT = KI+NX*(KJ-1)
              RIGHT = LEFT+1
              IT = LEFT
              X = DX+DELTAX*I
              CALL FREDGO(ARRAY,MN)
  100     CONTINUE
C
          X = DX+DELTAX*NX
          IF(.NOT.VISBLE) CALL PGMOVE(X,ARRAY(1,KJ))
          CALL PGDRAW(X,DELTAY*(KJ+1))
C               (set DELTAX positive for return trip)
          DELTAX = -DELTAX
  110 CONTINUE
C
C Now do the columns:
C as we fell out of the last DO-loop we do the
C columns in ascending-X order
C
      INCX = 1
      KI = 1
C               (set DELTAX -ve since scanning R to L)
  120 DX = DELTAX*(KI+NY-1)
      DELTAX = -DELTAX
      X = DX+DELTAX
      CALL PGMOVE(X,ARRAY(1,1))
  130 VISBLE = .TRUE.
      DO 140 J=2,NY
          LEFT = KI+NX*(J-1)
          RIGHT = LEFT-NX
          IT = LEFT
          X = DX+DELTAX*J
          CALL FREDGO(ARRAY,MN)
  140 CONTINUE
C
C At far end, increment X and check still inside array
C
      KI = KI+INCX
      IF(KI.LE.0 .OR. KI.GT.NX) GOTO 180
      VISBLE = ARRAY(KI,NY).GE.0.0
      DELTAX = -DELTAX
      DX = DELTAX*(KI-2)
      X = DX+DELTAX
      IF(VISBLE) CALL PGMOVE(X,ARRAY(KI,NY))
      DO 150 J=2,NY
          KJ = MY-J
          RIGHT = KI+NX*(KJ-1)
          LEFT = RIGHT+NX
          IT = RIGHT
          X = DX+DELTAX*J
          CALL FREDGO(ARRAY,MN)
  150 CONTINUE
C
      X = DX+DELTAX*NY
      IF(.NOT.VISBLE) CALL PGMOVE(X,ARRAY(KI,1))
      IF(KI.EQ.1) GOTO 180
      CALL PGDRAW(X,DELTAY*(KI+1))
      KI = KI+INCX
      IF(KI.GT.NX) GOTO 180
      IF(KI.EQ.1) GOTO 120
  160 DELTAX = -DELTAX
      DX = DELTAX*(1-KI-NY)
      X = DX+DELTAX
      CALL PGMOVE(X,DELTAY*(KI+1))
      CALL PGDRAW(X,ARRAY(KI,1))
      GOTO 130
C
C Do columns backwards because ended rows at far end of X
C
  170 KI = NX
      INCX = -1
      DX = DELTAX*(KI+NY)
      GOTO 160
C
C
  180 CALL PGEBUF
      END
C-----------------------------------------------------------------------
      SUBROUTINE FREDGO(ARRAY,MN)
      INTEGER MN
      REAL ARRAY(MN)
C
      INTEGER STEP,LEFT,RIGHT,IT,NX
      LOGICAL VISBLE
      REAL AL,AR,BL,EM,XX,X,Y,DELTAX
      COMMON /FREDCM/ DELTAX,X,STEP,LEFT,RIGHT,IT,NX,VISBLE
C
C Test visibility
C
      IF(ARRAY(IT).LT.0.0) GOTO 80
C
C This point is visible - was last?
C
      IF(VISBLE) GOTO 50
C
C No: calculate point where this line vanishes
C
   10 IF(LEFT.LE.NX .OR. MOD(LEFT-1,NX).EQ.0 .OR.
     1     RIGHT.LE.NX .OR. MOD(RIGHT-1,NX).EQ.0) GOTO 100
      AL = ABS(ARRAY(LEFT))
      AR = ABS(ARRAY(RIGHT))
      IF(ARRAY(LEFT).LT.0.0) GOTO 70
C               Right-hand point is crested
   20 RIGHT = RIGHT-STEP
      IF(ARRAY(RIGHT).LT.0.0) GOTO 20
C               Left-hand end of cresting line is either
C               RIGHT+NX or RIGHT-1
      LEFT = RIGHT+NX
      IF(ARRAY(LEFT).LT.0.0) LEFT = RIGHT-1
C
C               RIGHT and LEFT index into the endpoints of the
C               cresting line
   30 BL = ABS(ARRAY(LEFT))
      EM = ABS(ARRAY(RIGHT))-BL
      XX = EM-AR+AL
      IF(ABS(XX).LT.0.0001) GOTO 60
      XX = (AL-BL)/XX
   40 Y = EM*XX+BL
      IF(DELTAX.GT.0.0) XX = 1.0-XX
      XX = X-XX*DELTAX
      IF(VISBLE) GOTO 90
C               Drawing a line from an invisible point
C               to a visible one
      CALL PGMOVE(XX,Y)
      VISBLE = .TRUE.
   50 CALL PGDRAW(X,ARRAY(IT))
      RETURN
C
   60 XX = 0.5
      GOTO 40
C
C Left-hand point crested
C
   70 LEFT = LEFT-STEP
      IF(ARRAY(LEFT).LT.0.0) GOTO 70
C
C Right-hand end of cresting line is either LEFT+1 or LEFT-NX
C
      RIGHT = LEFT+1
      IF(ARRAY(RIGHT).LT.0.0) RIGHT = LEFT-NX
      GOTO 30
C
C This point is invisible; if last one was too, then forget it;
C else draw a line towards it
C
   80 IF(.NOT.VISBLE) RETURN
      GOTO 10
C
   90 CALL PGDRAW(XX,Y)
  100 VISBLE = .FALSE.
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C
C******************************************************************************
C
        SUBROUTINE SUBPMORE(STWV,DISP,RVEL)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INTEGER READILIM
C
        REAL RVEL
C
        REAL C                                             !velocidad de la luz
        PARAMETER (C=299792.46)
C
        INTEGER IOPC
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        REAL RCVEL1
        LOGICAL LCOLOR(MAX_ID_RED)
C
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        RCVEL1=1.0+RVEL/C
        RCVEL1=RCVEL1/SQRT(1.-(RVEL*RVEL)/(C*C))        !correccion relativista
C
        WRITE(*,101)'(1) plot indices'
        WRITE(*,101)'(2) plot typical lines'
        WRITE(*,101)'(3) plot photometric bands'
        WRITE(*,101)'(0) exit'
        WRITE(*,100)'Option '
        IOPC=READILIM('0',0,3)
C
        IF(IOPC.EQ.0)THEN
          RETURN
        ELSEIF(IOPC.EQ.1)THEN
          CALL SUBPINDEX(STWV,DISP,RCVEL1)
        ELSEIF(IOPC.EQ.2)THEN
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL SUBPLINES(STWV,DISP,RCVEL1,LCOLOR(ITERM))
          END DO
        ELSEIF(IOPC.EQ.3)THEN
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL SUBPBANDS(STWV,DISP,RCVEL1,LCOLOR(ITERM))
          END DO
        END IF
C
100     FORMAT(A,$)
101     FORMAT(A)
        END
C
C******************************************************************************
C
C Dibuja las bandas de los indices seleccionados
        SUBROUTINE SUBPINDEX(STWV,DISP,RCVEL1)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INTEGER READILIM
C
        REAL RCVEL1
C
        INTEGER K,IWL
        INTEGER NINDEX
        INTEGER NI,NINDEXT,ITI,NIND1,NIND2
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        INTEGER NCONTI,NABSOR
        REAL WLMIN,WLMAX
        REAL WLMIN0,WLMAX0
        REAL XMIN,XMAX,YMIN,YMAX
        REAL X1,X2
        REAL WV(NWVMAX),FWV(NWVMAX/4)
        REAL DY
        CHARACTER*8 CLABEL
        LOGICAL LANY,LINDOK(NINDMAX)
        LOGICAL LCOLOR(MAX_ID_RED)
C
        COMMON/BLKDATA1/NSCAN,NCHAN
        COMMON/BLKLIMITS/XMIN,XMAX,YMIN,YMAX
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
C------------------------------------------------------------------------------
C evitamos algunos warnings de compilacion
        NCONTI=0
        NABSOR=0
C------------------------------------------------------------------------------
        WLMIN=STWV-DISP/2.
        WLMAX=STWV+REAL(NCHAN-1)*DISP+DISP/2.0
C
        DY=YMAX-YMIN
        CALL SELINDEX(0,WV,FWV,NINDEXT,CLABEL)
        LANY=.FALSE.
        DO NINDEX=1,NINDEXT
          CALL SELINDEX(NINDEX,WV,FWV,ITI,CLABEL)
          IF((ITI.EQ.1).OR.(ITI.EQ.2))THEN
            DO K=1,6
              WV(K)=WV(K)*RCVEL1
            END DO
          ELSEIF((ITI.EQ.3).OR.(ITI.EQ.4).OR.(ITI.EQ.5))THEN
            DO K=1,4
              WV(K)=WV(K)*RCVEL1
            END DO
          ELSEIF((ITI.GE.101).OR.(ITI.LE.9999))THEN
            NCONTI=(ITI/100)
            NABSOR=ITI-NCONTI*100
            DO K=1,2*(NCONTI+NABSOR)
              WV(K)=WV(K)*RCVEL1
            END DO
          ELSE
            WRITE(*,101)'FATAL ERROR: invalid ITI value.'
          END IF
          IF((ITI.EQ.1).OR.(ITI.EQ.2))THEN
            WLMIN0=WV(1)
            DO IWL=2,6
              WLMIN0=AMIN1(WLMIN0,WV(IWL))
            END DO
            WLMAX0=WV(1)
            DO IWL=2,6
              WLMAX0=AMAX1(WLMAX0,WV(IWL))
            END DO
            LINDOK(NINDEX)=((WLMIN0.GE.WLMIN).AND.(WLMAX0.LE.WLMAX))
          ELSEIF((ITI.EQ.3).OR.(ITI.EQ.4).OR.(ITI.EQ.5))THEN
            LINDOK(NINDEX)=((WV(1).GE.WLMIN).AND.(WV(4).LE.WLMAX))
          ELSEIF((ITI.GE.101).OR.(ITI.LE.9999))THEN
            WLMIN0=WV(1)
            DO IWL=2,2*(NCONTI+NABSOR)
              WLMIN0=AMIN1(WLMIN0,WV(IWL))
            END DO
            WLMAX0=WV(1)
            DO IWL=2,2*(NCONTI+NABSOR)
              WLMAX0=AMAX1(WLMAX0,WV(IWL))
            END DO
            LINDOK(NINDEX)=((WLMIN0.GE.WLMIN).AND.(WLMAX0.LE.WLMAX))
          END IF
          IF(LINDOK(NINDEX))THEN
            LANY=.TRUE.
          END IF
        END DO
C
        IF(LANY)THEN
          CALL SHINDEX(LINDOK,0)
          WRITE(*,100)'Index bands to be plotted'
          NINDEX=READILIM('@',-1,NINDEXT)
          IF(NINDEX.GE.0)THEN
            IF(NINDEX.EQ.0)THEN
              NIND1=1
              NIND2=NINDEXT
            ELSE
              NIND1=NINDEX
              NIND2=NINDEX
            END IF
            DO NI=NIND1,NIND2
              IF(LINDOK(NI))THEN
                CALL SELINDEX(NI,WV,FWV,ITI,CLABEL)
c..............................................................................
                IF((ITI.EQ.1).OR.(ITI.EQ.2))THEN
                  DO K=1,6
                    WV(K)=WV(K)*RCVEL1
                  END DO
                  DO ITERM=NTERM,1,-1
                    CALL PGSLCT(IDN(ITERM))
                    IF(LCOLOR(ITERM)) CALL PGSCI(4)
                    X1=(WV(1)-STWV)/DISP+1.
                    X2=(WV(2)-STWV)/DISP+1.
                    CALL PGMOVE(X1,YMIN)
                    CALL PGDRAW(X1,YMAX)
                    CALL PGMOVE(X2,YMIN)
                    CALL PGDRAW(X2,YMAX)
                    CALL PGRECT(X1,X2,YMIN,YMIN+DY/100.)
                    CALL PGRECT(X1,X2,YMAX,YMAX-DY/100.)
                    IF(LCOLOR(ITERM)) CALL PGSCI(3)
                    X1=(WV(3)-STWV)/DISP+1.
                    X2=(WV(4)-STWV)/DISP+1.
                    CALL PGMOVE(X1,YMIN)
                    CALL PGDRAW(X1,YMAX)
                    CALL PGMOVE(X2,YMIN)
                    CALL PGDRAW(X2,YMAX)
                    CALL PGRECT(X1,X2,YMIN,YMIN+DY/100.)
                    CALL PGRECT(X1,X2,YMAX,YMAX-DY/100.)
                    CALL PGPTEXT((X1+X2)/2.,YMAX-DY/20.,0.,.5,CLABEL)
                    IF(LCOLOR(ITERM)) CALL PGSCI(2)
                    X1=(WV(5)-STWV)/DISP+1.
                    X2=(WV(6)-STWV)/DISP+1.
                    CALL PGMOVE(X1,YMIN)
                    CALL PGDRAW(X1,YMAX)
                    CALL PGMOVE(X2,YMIN)
                    CALL PGDRAW(X2,YMAX)
                    CALL PGRECT(X1,X2,YMIN,YMIN+DY/100.)
                    CALL PGRECT(X1,X2,YMAX,YMAX-DY/100.)
                    IF(LCOLOR(ITERM)) CALL PGSCI(1)
                  END DO
c..............................................................................
                ELSEIF((ITI.EQ.3).OR.(ITI.EQ.4).OR.(ITI.EQ.5))THEN
                  DO K=1,4
                    WV(K)=WV(K)*RCVEL1
                  END DO
                  DO ITERM=NTERM,1,-1
                    CALL PGSLCT(IDN(ITERM))
                    IF(LCOLOR(ITERM)) CALL PGSCI(4)
                    X1=(WV(1)-STWV)/DISP+1.
                    X2=(WV(2)-STWV)/DISP+1.
                    CALL PGMOVE(X1,YMIN)
                    CALL PGDRAW(X1,YMAX)
                    CALL PGMOVE(X2,YMIN)
                    CALL PGDRAW(X2,YMAX)
                    CALL PGRECT(X1,X2,YMIN,YMIN+DY/100.)
                    CALL PGRECT(X1,X2,YMAX,YMAX-DY/100.)
                    CALL PGPTEXT((X1+X2)/2.,YMAX-DY/20.,0.,.5,CLABEL)
                    IF(LCOLOR(ITERM)) CALL PGSCI(2)
                    X1=(WV(3)-STWV)/DISP+1.
                    X2=(WV(4)-STWV)/DISP+1.
                    CALL PGMOVE(X1,YMIN)
                    CALL PGDRAW(X1,YMAX)
                    CALL PGMOVE(X2,YMIN)
                    CALL PGDRAW(X2,YMAX)
                    CALL PGRECT(X1,X2,YMIN,YMIN+DY/100.)
                    CALL PGRECT(X1,X2,YMAX,YMAX-DY/100.)
                    CALL PGPTEXT((X1+X2)/2.,YMAX-DY/20.,0.,.5,CLABEL)
                    IF(LCOLOR(ITERM)) CALL PGSCI(1)
                  END DO
c..............................................................................
                ELSEIF((ITI.GE.101).OR.(ITI.LE.9999))THEN
                  NCONTI=(ITI/100)
                  NABSOR=ITI-NCONTI*100
                  DO K=1,2*(NCONTI+NABSOR)
                    WV(K)=WV(K)*RCVEL1
                  END DO
                  DO ITERM=NTERM,1,-1
                    CALL PGSLCT(IDN(ITERM))
                    DO K=1,NCONTI+NABSOR
                      IF(K.LE.NCONTI)THEN
                        IF(LCOLOR(ITERM)) CALL PGSCI(5)
                      ELSE
                        IF(LCOLOR(ITERM)) CALL PGSCI(3)
                      END IF
                      X1=(WV(2*K-1)-STWV)/DISP+1.
                      X2=(WV(2*K)-STWV)/DISP+1.
                      CALL PGMOVE(X1,YMIN)
                      IF(K.LE.NCONTI)THEN
                        CALL PGDRAW(X1,YMAX)
                      ELSE
                        CALL PGDRAW(X1,YMAX-(YMAX-YMIN)/40.)
                      END IF
                      CALL PGMOVE(X2,YMIN)
                      IF(K.LE.NCONTI)THEN
                        CALL PGDRAW(X2,YMAX)
                      ELSE
                        CALL PGDRAW(X2,YMAX-(YMAX-YMIN)/40.)
                      END IF
                      CALL PGRECT(X1,X2,YMIN,YMIN+DY/100.)
                      IF(K.LE.NCONTI)THEN
                        CALL PGRECT(X1,X2,YMAX,YMAX-DY/100.)
                      ELSE
                        CALL PGRECT(X1,X2,YMAX-(YMAX-YMIN)/40.,
     +                   YMAX-(YMAX-YMIN)/40.-(YMAX-YMIN)/100.)
                      END IF
                      IF(K.GT.NCONTI) 
     +                 CALL PGPTEXT((X1+X2)/2.,YMAX-DY/10.,0.,.5,CLABEL)
                      IF(LCOLOR(ITERM)) CALL PGSCI(1)
                    END DO
                  END DO
c..............................................................................
                END IF
              END IF
            END DO
          END IF
        ELSE
          WRITE(*,101)'ERROR: no indices available in current '//
     +     'wavelength range!'
        END IF
C
100     FORMAT(A,$)
101     FORMAT(A)
        END
C
C******************************************************************************
C
C Dibuja las bandas fotometricas
        SUBROUTINE SUBPBANDS(STWV,DISP,RCVEL1,LCOLOR)
        IMPLICIT NONE
C
        INCLUDE 'redlib.inc'
        REAL RCVEL1
        LOGICAL LCOLOR
C
        INTEGER NPBAND,I
        REAL WV(NPBANDMAX),RES(NPBANDMAX)
        REAL WLMIN,WLMAX
        REAL XX(NPBANDMAX)
        REAL XMIN,XMAX,YMIN,YMAX
C
        COMMON/BLKDATA1/NSCAN,NCHAN
        COMMON/BLKLIMITS/XMIN,XMAX,YMIN,YMAX
C------------------------------------------------------------------------------
        WLMIN=STWV-DISP/2.
        WLMAX=STWV+REAL(NCHAN-1)*DISP+DISP/2.0
C
        CALL PGSLS(4)
        CALL PGSCH(0.8)
C
C Banda U
        CALL SELBANDS('U',NPBAND,WV,RES)
        IF(LCOLOR) CALL PGSCI(4)
        DO I=1,NPBAND
          XX(I)=(WV(I)*RCVEL1-STWV)/DISP+1.
          RES(I)=RES(I)*(YMAX*0.30)
        END DO
        CALL PGLINE(NPBAND,XX,RES)
        DO I=1,NPBAND
          IF((WV(I).GE.WLMIN).AND.(WV(I).LE.WLMAX))THEN
            CALL PGPTEXT(XX(I),RES(I),0.,0.5,'U')
          END IF
        END DO
C Banda B
        CALL SELBANDS('B',NPBAND,WV,RES)
        IF(LCOLOR) CALL PGSCI(3)
        DO I=1,NPBAND
          XX(I)=(WV(I)*RCVEL1-STWV)/DISP+1.
          RES(I)=RES(I)*(YMAX*0.30)
        END DO
        CALL PGLINE(NPBAND,XX,RES)
        DO I=1,NPBAND
          IF((WV(I).GE.WLMIN).AND.(WV(I).LE.WLMAX))THEN
            CALL PGPTEXT(XX(I),RES(I),0.,0.5,'B')
          END IF
        END DO
C Banda V
        CALL SELBANDS('V',NPBAND,WV,RES)
        IF(LCOLOR) CALL PGSCI(2)
        DO I=1,NPBAND
          XX(I)=(WV(I)*RCVEL1-STWV)/DISP+1.
          RES(I)=RES(I)*(YMAX*0.30)
        END DO
        CALL PGLINE(NPBAND,XX,RES)
        DO I=1,NPBAND
          IF((WV(I).GE.WLMIN).AND.(WV(I).LE.WLMAX))THEN
            CALL PGPTEXT(XX(I),RES(I),0.,0.5,'V')
          END IF
        END DO
C
        IF(LCOLOR) CALL PGSCI(1)
        CALL PGSCH(1.0)
        CALL PGSLS(1)
C
        END
C
C******************************************************************************
C
C Dibuja las bandas de los indices seleccionados
        SUBROUTINE SUBPLINES(STWV,DISP,RCVEL1,LCOLOR)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        REAL RCVEL1
        LOGICAL LCOLOR
C
        INTEGER NLINE,NLINEST
        REAL WLMIN,WLMAX
        REAL XMIN,XMAX,YMIN,YMAX
        REAL X0,X0MIN,X0MAX,Y0MIN,Y0MAX
        REAL WV(NLINMAX)
        REAL DY
        CHARACTER*8 CLABEL(NLINMAX)
C
        COMMON/BLKDATA1/NSCAN,NCHAN
        COMMON/BLKLIMITS/XMIN,XMAX,YMIN,YMAX
C------------------------------------------------------------------------------
        CALL PGQWIN(X0MIN,X0MAX,Y0MIN,Y0MAX)
        DY=YMAX-YMIN
C
        WLMIN=STWV-DISP/2.
        WLMAX=STWV+REAL(NCHAN-1)*DISP+DISP/2.0
C
        CALL SELLINES(0,NLINEST,WV,CLABEL)                     !serie de Balmer
        IF(LCOLOR) CALL PGSCI(5)
        DO NLINE=1,NLINEST
          WV(NLINE)=WV(NLINE)*RCVEL1
          IF((WV(NLINE).GE.WLMIN).AND.(WV(NLINE).LE.WLMAX))THEN
            X0=(WV(NLINE)-STWV)/DISP+1.
            IF((X0.GE.X0MIN).AND.(X0.LE.X0MAX))THEN
              CALL PGMOVE(X0,YMAX)
              CALL PGDRAW(X0,YMIN)
              CALL PGSCH(0.8)
              CALL PGPTEXT(X0,YMAX+DY/100.,0.,.5,CLABEL(NLINE))
              CALL PGSCH(1.0)
            END IF
          END IF
        END DO
        IF(LCOLOR) CALL PGSCI(6)
        CALL SELLINES(1,NLINEST,WV,CLABEL)                        !otras lineas
        DO NLINE=1,NLINEST
          WV(NLINE)=WV(NLINE)*RCVEL1
          IF((WV(NLINE).GE.WLMIN).AND.(WV(NLINE).LE.WLMAX))THEN
            X0=(WV(NLINE)-STWV)/DISP+1.
            IF((X0.GE.X0MIN).AND.(X0.LE.X0MAX))THEN
              CALL PGMOVE(X0,YMAX)
              CALL PGDRAW(X0,YMIN)
              CALL PGSCH(0.8)
              CALL PGPTEXT(X0,YMAX+DY/100.,0.,.5,CLABEL(NLINE))
              CALL PGSCH(1.0)
            END IF
          END IF
        END DO
        IF(LCOLOR) CALL PGSCI(7)
        CALL SELLINES(2,NLINEST,WV,CLABEL)            !lineas tipicas del cielo
        DO NLINE=1,NLINEST
C NOTA: las lineas de cielo no hay que corregirlas de redshift (obviamente)
ccc       WV(NLINE)=WV(NLINE)*RCVEL1
          IF((WV(NLINE).GE.WLMIN).AND.(WV(NLINE).LE.WLMAX))THEN
            X0=(WV(NLINE)-STWV)/DISP+1.
            IF((X0.GE.X0MIN).AND.(X0.LE.X0MAX))THEN
              CALL PGMOVE(X0,YMAX)
              CALL PGDRAW(X0,YMIN)
              CALL PGSCH(0.8)
              CALL PGPTEXT(X0,YMIN+DY/100.,0.,.5,CLABEL(NLINE))
              CALL PGSCH(1.0)
            END IF
          END IF
        END DO
        IF(LCOLOR) CALL PGSCI(8)
        CALL SELLINES(3,NLINEST,WV,CLABEL)         !lineas tipicas de absorción
        DO NLINE=1,NLINEST
          WV(NLINE)=WV(NLINE)*RCVEL1
          IF((WV(NLINE).GE.WLMIN).AND.(WV(NLINE).LE.WLMAX))THEN
            X0=(WV(NLINE)-STWV)/DISP+1.
            IF((X0.GE.X0MIN).AND.(X0.LE.X0MAX))THEN
              CALL PGMOVE(X0,YMAX)
              CALL PGDRAW(X0,YMIN)
              CALL PGSCH(0.8)
              CALL PGPTEXT(X0,YMIN+DY/100.,0.,.5,CLABEL(NLINE))
              CALL PGSCH(1.0)
            END IF
          END IF
        END DO
        IF(LCOLOR) CALL PGSCI(1)
C
        END
C
C******************************************************************************
C Chequea si las nuevas dimensiones son compatibles con las ya existentes
        SUBROUTINE CCSIZE(NSCAN,NCHAN,STWV,DISP,
     +   NSCAN2,NCHAN2,STWV2,DISP2,LOK)
        IMPLICIT NONE
        INTEGER NSCAN,NCHAN
        INTEGER NSCAN2,NCHAN2
        REAL STWV,DISP
        REAL STWV2,DISP2
        LOGICAL LOK
C------------------------------------------------------------------------------
        LOK=.TRUE.
        IF(NSCAN2.NE.NSCAN)THEN
          WRITE(*,100)'ERROR: invalid NSCAN: '
          WRITE(*,*)NSCAN2
          WRITE(*,100)'> expected NSCAN: '
          WRITE(*,*)NSCAN
          WRITE(*,101)'Last file will not be read'
          WRITE(*,*)
          WRITE(*,100)'Press <CR> to continue...'
          READ(*,*)
          LOK=.FALSE.
          RETURN
        END IF
        IF(NCHAN2.NE.NCHAN)THEN
          WRITE(*,100)'ERROR: invalid NCHAN: '
          WRITE(*,*)NCHAN2
          WRITE(*,100)'> expected NCHAN: '
          WRITE(*,*)NCHAN
          WRITE(*,101)'Last file will not be read'
          WRITE(*,*)
          WRITE(*,100)'Press <CR> to continue...'
          READ(*,*)
          LOK=.FALSE.
          RETURN
        END IF
        IF(STWV2.NE.STWV)THEN
          WRITE(*,100)'WARNING: new STWV: '
          WRITE(*,*)STWV2
          WRITE(*,100)'> expected STWV: '
          WRITE(*,*)STWV
          WRITE(*,101)'Last STWV value will not be employed'
          WRITE(*,*)
          WRITE(*,100)'Press <CR> to continue...'
          READ(*,*)
        END IF
        IF(DISP2.NE.DISP)THEN
          WRITE(*,100)'WARNING: new DISP: '
          WRITE(*,*)DISP2
          WRITE(*,100)'> expected DISP: '
          WRITE(*,*)DISP
          WRITE(*,101)'Last DISP value will not be employed'
          WRITE(*,*)
          WRITE(*,100)'Press <CR> to continue...'
          READ(*,*)
        END IF
100     FORMAT(A,$)
101     FORMAT(A)
C
        END
C
C******************************************************************************
C******************************************************************************
C
        SUBROUTINE SUBIMATH
        IMPLICIT NONE
C
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INTEGER TRUELEN
        REAL READF
C
        INTEGER NYBUTT
        INTEGER NBUFF,NBUFF0,NBUFF1,NB
        INTEGER I,J,L,K
        INTEGER NC1,NC2,NS1,NS2
        INTEGER NSCAN2,NCHAN2
        REAL A(NCMAX,NSMAX,NMAXBUFF)
        REAL SSX(NCMAX),SSY(NSMAX)
        REAL XCUT(NCMAX),YCUT(NSMAX)
        REAL X1B,X2B,Y1B,Y2B
        REAL XC,YC
        REAL RFACTOR,FMEANSIGMA,FMEDIAN,FMEDIAN1
        REAL STWV2,DISP2
        REAL FMEAN1,FMEAN2
        CHARACTER*1 CH
        CHARACTER*1 COPER,CUTIL,CMODE
        CHARACTER*50 CDUMMY
        CHARACTER*75 INFILEBUFF(NMAXBUFF),NEWFILE
        LOGICAL LDEFBUFF(NMAXBUFF),LUSEBUFF(NMAXBUFF)
        LOGICAL LOK,LDIVZERO
        LOGICAL IFSCAN(NSMAX),IFCHAN(NCMAX)
C
        COMMON/BLKDATA1/NSCAN,NCHAN
        COMMON/BLKDATA2/A
        COMMON/BLKBUFF1/LDEFBUFF,LUSEBUFF
        COMMON/BLKBUFF2/INFILEBUFF
C------------------------------------------------------------------------------
        OUTFILEX=OUTFILEX
C
        CALL BUTTQBR(X1B,X2B,Y1B,Y2B)             !initial button region limits
        CALL BUTTQYB(NYBUTT)                                 !initial MAX_YBUTT
C
        CALL BUTTSBR(X1B,X2B,0.05,Y2B)
        CALL BUTTSYB(8)
C------------------------------------------------------------------------------
5       CALL PGERAS
        WRITE(*,*)
        DO NB=1,48
          CALL BUTTSEX(NB,.FALSE.)
        END DO
C
        CALL BUTTON(30,'cancel',0)
        CALL BUTTON(30,'cancel',3)
        CALL BUTTON(36,'go',0)
        CALL BUTTON(36,'go',3)
        CALL BUTTON(42,'exit',0)
C------------------------------------------------------------------------------
C Mostramos todos los buffers y seleccionamos uno (NBUFF0)
        DO NBUFF=1,NMAXBUFF
          NB=REAL(NBUFF)*6+1
          IF(LDEFBUFF(NBUFF))THEN
            WRITE(CDUMMY,'(A8,I1)')'Buffer #',NBUFF
            CALL BUTTON(NB,CDUMMY(1:9),0)
          ELSE
            CALL BUTTON(NB,'undefined',0)
          END IF
        END DO
C
10      CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
        CALL IFBUTTON(XC,YC,NB)
C
        IF(NB.EQ.42)THEN
          GOTO 900
        ELSE
          IF(MOD(NB-7,6).EQ.0)THEN
            NBUFF0=(NB-7)/6+1
          ELSE
            GOTO 10
          END IF
        END IF
C
        DO NBUFF=1,NMAXBUFF
          NB=REAL(NBUFF)*6+1
          IF(LDEFBUFF(NBUFF))THEN
            WRITE(CDUMMY,'(A8,I1)')'Buffer #',NBUFF
            IF(NBUFF0.EQ.NBUFF)THEN
              CALL BUTTON(NB,CDUMMY(1:9),1)
            ELSE
              CALL BUTTON(NB,CDUMMY(1:9),3)
            END IF
          ELSE
            IF(NBUFF0.EQ.NBUFF)THEN
              WRITE(CDUMMY,'(A8,I1)')'Buffer #',NBUFF
              CALL BUTTON(NB,CDUMMY(1:9),0)
              CALL BUTTON(NB,CDUMMY(1:9),1)
              DO I=1,NSCAN               !inicializamos a cero el buffer NBUFF0
                DO J=1,NCHAN
                  A(J,I,NBUFF0)=0.0
                END DO
              END DO
            ELSE
              CALL BUTTON(NB,'undefined',3)
            END IF
          END IF
        END DO
C------------------------------------------------------------------------------
C Mostramos operaciones posibles con el buffer seleccionado y seleccionamos
C la operacion deseada (COPER). Si la opcion es limpiar el buffer seleccionado,
C dicha accion se realiza inmediatamente.
        CALL BUTTON(30,'cancel',2)

        IF(LDEFBUFF(NBUFF0)) CALL BUTTON(02,'clear buffer',0)
        CALL BUTTON(08,'+',0)
        CALL BUTTON(14,'-',0)
        CALL BUTTON(20,'*',0)
        CALL BUTTON(26,'/',0)
C
20      CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
        CALL IFBUTTON(XC,YC,NB)
C
        IF(NB.EQ.30)THEN
          GOTO 5
        ELSEIF(NB.EQ.42)THEN
          GOTO 900
        ELSEIF(NB.EQ.02)THEN
          COPER='0'
        ELSEIF(NB.EQ.08)THEN
          COPER='+'
        ELSEIF(NB.EQ.14)THEN
          COPER='-'
        ELSEIF(NB.EQ.20)THEN
          COPER='*'
        ELSEIF(NB.EQ.26)THEN
          COPER='/'
        ELSE
          GOTO 20
        END IF
C
        IF(LDEFBUFF(NBUFF0)) CALL BUTTON(02,'clear buffer',3)
        CALL BUTTON(08,'+',3)
        CALL BUTTON(14,'-',3)
        CALL BUTTON(20,'*',3)
        CALL BUTTON(26,'/',3)
C
        IF(COPER.EQ.'0')THEN
          CALL BUTTON(02,'clear buffer',1)
          CALL BUTTON(36,'go',2)
          CALL BUTTON(36,'go',-3)
22        CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
          CALL IFBUTTON(XC,YC,NB)
          IF(NB.EQ.30)THEN
            GOTO 5
          ELSEIF(NB.EQ.36)THEN
            DO I=1,NSCAN
              DO J=1,NCHAN
                A(J,I,NBUFF0)=0.0
              END DO
            END DO
            WRITE(*,'(A,I1,A)')'>>> Buffer #',NBUFF0,
     +       ' has been cleared.'
            INFILEBUFF(NBUFF0)='[cleared]'
            GOTO 5
          ELSEIF(NB.EQ.42)THEN
            GOTO 900
          ELSE
            GOTO 22
          END IF
        ELSEIF(COPER.EQ.'+')THEN
          CALL BUTTON(08,'+',1)
        ELSEIF(COPER.EQ.'-')THEN
          CALL BUTTON(14,'-',1)
        ELSEIF(COPER.EQ.'*')THEN
          CALL BUTTON(20,'*',1)
        ELSEIF(COPER.EQ.'/')THEN
          CALL BUTTON(26,'/',1)
        END IF
C------------------------------------------------------------------------------
C Determinamos si la accion se va a realizar sobre todo el buffer NBUFF0 o
C sobre una region. En cualquier caso, la region a utilizar queda definida por
C NS1,NS2,NC1,NC2.
        CALL BUTTON(32,'full frame',0)
        CALL BUTTON(38,'region',0)
C
25      CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
        CALL IFBUTTON(XC,YC,NB)
C
        IF(NB.EQ.30)THEN
          GOTO 5
        ELSEIF(NB.EQ.42)THEN
          GOTO 900
        ELSEIF(NB.EQ.32)THEN
          CALL BUTTON(32,'full frame',1)
          CALL BUTTON(38,'region',3)
          NS1=1
          NS2=NSCAN
          NC1=1
          NC2=NCHAN
        ELSEIF(NB.EQ.38)THEN
          CALL BUTTON(32,'full frame',3)
          CALL BUTTON(38,'region',1)
          WRITE(CDUMMY,'(A,I10)')'1,',NSCAN
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          WRITE(*,100)'1st and last scan '
          CALL READ2I(CDUMMY(1:L),NS1,NS2)
          WRITE(CDUMMY,'(A,I10)')'1,',NCHAN
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          WRITE(*,100)'1st and last channel '
          CALL READ2I(CDUMMY(1:L),NC1,NC2)
        ELSE
          GOTO 25
        END IF
C------------------------------------------------------------------------------
C Mostramos los posibles operandos a utilizar sobre NBUFF0 (CUTIL). Si la
C opcion es operar con una constante, la accion se realiza inmediatamente.
        CALL BUTTON(09,'constant',0)
        CALL BUTTON(15,'frame',0)
        CALL BUTTON(21,'x-cut',0)
        CALL BUTTON(27,'y-cut',0)
C
30      CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
        CALL IFBUTTON(XC,YC,NB)
C
        IF(NB.EQ.30)THEN
          GOTO 5
        ELSEIF(NB.EQ.42)THEN
          GOTO 900
        ELSEIF(NB.EQ.09)THEN
          CUTIL='c'
        ELSEIF(NB.EQ.15)THEN
          CUTIL='f'
        ELSEIF(NB.EQ.21)THEN
          CUTIL='x'
        ELSEIF(NB.EQ.27)THEN
          CUTIL='y'
        ELSE
          GOTO 30
        END IF
C
        CALL BUTTON(09,'constant',3)
        CALL BUTTON(15,'frame',3)
        CALL BUTTON(21,'x-cut',3)
        CALL BUTTON(27,'y-cut',3)
C
        IF(CUTIL.EQ.'c')THEN
          CALL BUTTON(09,'constant',1)
31        WRITE(*,100)'Constant factor'
          RFACTOR=READF('@')
          IF((COPER.EQ.'/').AND.(RFACTOR.EQ.0.0))THEN
            WRITE(*,101)'ERROR: division by zero. '
            WRITE(*,100)'Insert a new constant.'
            GOTO 31
          END IF
          CALL BUTTON(36,'go',2)
          CALL BUTTON(36,'go',-3)
32        CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
          CALL IFBUTTON(XC,YC,NB)
          IF(NB.EQ.30)THEN
            GOTO 5
          ELSEIF(NB.EQ.36)THEN
            IF(COPER.EQ.'+')THEN
              DO I=NS1,NS2
                DO J=NC1,NC2
                  A(J,I,NBUFF0)=A(J,I,NBUFF0)+RFACTOR
                END DO
              END DO
            ELSEIF(COPER.EQ.'-')THEN
              DO I=NS1,NS2
                DO J=NC1,NC2
                  A(J,I,NBUFF0)=A(J,I,NBUFF0)-RFACTOR
                END DO
              END DO
            ELSEIF(COPER.EQ.'*')THEN
              DO I=NS1,NS2
                DO J=NC1,NC2
                  A(J,I,NBUFF0)=A(J,I,NBUFF0)*RFACTOR
                END DO
              END DO
            ELSEIF(COPER.EQ.'/')THEN
              DO I=NS1,NS2
                DO J=NC1,NC2
                  A(J,I,NBUFF0)=A(J,I,NBUFF0)/RFACTOR
                END DO
              END DO
            END IF
            WRITE(*,'(A,I1,A)')'>>> Buffer #',NBUFF0,
     +       ' has been modified.'
            LDEFBUFF(NBUFF0)=.TRUE.
            LUSEBUFF(NBUFF0)=.TRUE.
            L=TRUELEN(INFILEBUFF(NBUFF0))
            IF(L.LT.75)INFILEBUFF(NBUFF0)(L+1:L+1)='*'
            GOTO 5
          ELSEIF(NB.EQ.42)THEN
            GOTO 900
          ELSE
            GOTO 32
          END IF
        ELSEIF(CUTIL.EQ.'f')THEN
          CALL BUTTON(15,'frame',1)
        ELSEIF(CUTIL.EQ.'x')THEN
          CALL BUTTON(21,'x-cut',1)
        ELSEIF(CUTIL.EQ.'y')THEN
          CALL BUTTON(27,'y-cut',1)
        END IF
C------------------------------------------------------------------------------
C Si el operando no es una constante, podemos escoger entre una nueva imagen
C o un buffer ya definido (NBUFF1).
        DO NBUFF=1,NMAXBUFF
          NB=REAL(NBUFF)*6+4
          IF(LDEFBUFF(NBUFF))THEN
            WRITE(CDUMMY,'(A8,I1)')'Buffer #',NBUFF
            CALL BUTTON(NB,CDUMMY(1:9),0)
          ELSE
            CALL BUTTON(NB,'undefined',0)
            CALL BUTTON(NB,'undefined',3)
          END IF
        END DO
        CALL BUTTON(4,'new file',0)
C
40      CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
        CALL IFBUTTON(XC,YC,NB)
C
        IF(NB.EQ.30)THEN
          GOTO 5
        ELSEIF(NB.EQ.42)THEN
          GOTO 900
        ELSE
          IF(NB.EQ.4)THEN
            NBUFF1=0
          ELSEIF(NB.EQ.10)THEN
            NBUFF1=1
          ELSEIF(NB.EQ.16)THEN
            NBUFF1=2
          ELSEIF(NB.EQ.22)THEN
            NBUFF1=3
          ELSEIF(NB.EQ.28)THEN
            NBUFF1=4
          ELSEIF(NB.EQ.34)THEN
            NBUFF1=5
          ELSEIF(NB.EQ.40)THEN
            NBUFF1=6
          ELSE
            GOTO 40
          END IF
        END IF
C
        IF(NBUFF1.EQ.0)THEN
          CALL BUTTON(4,'new file',1)
        ELSE
          CALL BUTTON(4,'new file',3)
        END IF
        DO NBUFF=1,NMAXBUFF
          NB=REAL(NBUFF)*6+4
          IF(LDEFBUFF(NBUFF))THEN
            WRITE(CDUMMY,'(A8,I1)')'Buffer #',NBUFF
            IF(NBUFF1.EQ.NBUFF)THEN
              CALL BUTTON(NB,CDUMMY(1:9),1)
            ELSE
              CALL BUTTON(NB,CDUMMY(1:9),3)
            END IF
          END IF
        END DO
C------------------------------------------------------------------------------
C Opcion operar con un frame: la accion se realiza directamente.
        IF(CUTIL.EQ.'f')THEN
C..............................................................................
          CALL BUTTON(36,'go',2)
          CALL BUTTON(36,'go',-3)
C
45        CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
          CALL IFBUTTON(XC,YC,NB)
C
          IF(NB.EQ.30)THEN
            GOTO 5
          ELSEIF(NB.EQ.42)THEN
            GOTO 900
          ELSE
            IF(NB.NE.36) GOTO 45
          END IF
C..............................................................................
          IF(NBUFF1.EQ.0)THEN                           !usamos fichero externo
            WRITE(*,100)'Name of the external file'
            NEWFILE=INFILEX(20,'@',NSCAN2,NCHAN2,STWV2,DISP2,1,.FALSE.)
            CALL CCSIZE(NSCAN,NCHAN,STWV,DISP,
     +       NSCAN2,NCHAN2,STWV2,DISP2,LOK)
            IF(.NOT.LOK)THEN
              CLOSE(20)
              GOTO 5
            END IF
            IF(NS1.GT.1)THEN        !saltamos los scans que no vamos a utilizar
              DO I=1,NS1-1
                READ(20)(XCUT(J),J=1,NCHAN)
              END DO
            END IF
            DO I=NS1,NS2                  !realizamos la operacion seleccionada
              READ(20)(XCUT(J),J=1,NCHAN)
              IF(COPER.EQ.'/')THEN
                LDIVZERO=.FALSE.
                DO J=NC1,NC2
                  IF(XCUT(J).EQ.0.0) LDIVZERO=.TRUE.
                END DO
                IF(LDIVZERO)THEN
                  WRITE(*,100)'ERROR: division by zero '
                  WRITE(*,101)'at scan: '
                  WRITE(*,*)I
                  WRITE(*,100)'Press <CR> to continue...'
                  READ(*,*)
                  CLOSE(20)
                  GOTO 5
                END IF
              END IF
              IF(COPER.EQ.'+')THEN
                DO J=NC1,NC2
                  A(J,I,NBUFF0)=A(J,I,NBUFF0)+XCUT(J)
                END DO
              ELSEIF(COPER.EQ.'-')THEN
                DO J=NC1,NC2
                  A(J,I,NBUFF0)=A(J,I,NBUFF0)-XCUT(J)
                END DO
              ELSEIF(COPER.EQ.'*')THEN
                DO J=NC1,NC2
                  A(J,I,NBUFF0)=A(J,I,NBUFF0)*XCUT(J)
                END DO
              ELSEIF(COPER.EQ.'/')THEN
                DO J=NC1,NC2
                  A(J,I,NBUFF0)=A(J,I,NBUFF0)/XCUT(J)
                END DO
              END IF
            END DO
            CLOSE(20)
            WRITE(*,'(A,I1,A)')'>>> Buffer #',NBUFF0,
     +       ' has been modified.'
            LDEFBUFF(NBUFF0)=.TRUE.      !el buffer se activa (si no lo estaba)
            LUSEBUFF(NBUFF0)=.TRUE.
            L=TRUELEN(INFILEBUFF(NBUFF0))
            IF(L.LT.75)INFILEBUFF(NBUFF0)(L+1:L+1)='*'
C..............................................................................
          ELSE                                       !usamos buffer ya definido
            DO I=NS1,NS2                  !realizamos la operacion seleccionada
              IF(COPER.EQ.'/')THEN
                LDIVZERO=.FALSE.
                DO J=NC1,NC2
                  IF(A(J,I,NBUFF1).EQ.0.0) LDIVZERO=.TRUE.
                END DO
                IF(LDIVZERO)THEN
                  WRITE(*,100)'ERROR: division by zero '
                  WRITE(*,101)'at scan: '
                  WRITE(*,*)I
                  WRITE(*,100)'Press <CR> to continue...'
                  READ(*,*)
                  GOTO 5
                END IF
              END IF
              IF(COPER.EQ.'+')THEN
                DO J=NC1,NC2
                  A(J,I,NBUFF0)=A(J,I,NBUFF0)+A(J,I,NBUFF1)
                END DO
              ELSEIF(COPER.EQ.'-')THEN
                DO J=NC1,NC2
                  A(J,I,NBUFF0)=A(J,I,NBUFF0)-A(J,I,NBUFF1)
                END DO
              ELSEIF(COPER.EQ.'*')THEN
                DO J=NC1,NC2
                  A(J,I,NBUFF0)=A(J,I,NBUFF0)*A(J,I,NBUFF1)
                END DO
              ELSEIF(COPER.EQ.'/')THEN
                DO J=NC1,NC2
                  A(J,I,NBUFF0)=A(J,I,NBUFF0)/A(J,I,NBUFF1)
                END DO
              END IF
            END DO
            CLOSE(20)
            WRITE(*,'(A,I1,A)')'>>> Buffer #',NBUFF0,
     +       ' has been modified.'
            LDEFBUFF(NBUFF0)=.TRUE.      !el buffer se activa (si no lo estaba)
            LUSEBUFF(NBUFF0)=.TRUE.
            L=TRUELEN(INFILEBUFF(NBUFF0))
            IF(L.LT.75)INFILEBUFF(NBUFF0)(L+1:L+1)='*'
C..............................................................................
          END IF
          GOTO 5
        END IF
C------------------------------------------------------------------------------
C Si lo que vamos a utilizar como operando es un corte unidimensional, y dicho
C corte puede calcularse a partir de un frame (mediante colapso de una
C dimension), tenemos que seleccionar la region a colapsar.
C En caso contrario, realizamos directamente la operacion con el corte
C seleccionado.
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C Si es un corte de un fichero externo, solo podemos calcular (sin dimensionar
C mas matrices) la media o la suma total. En el caso de un buffer ya existente,
C podemos calcular tambien la media (eliminando puntos a varios sigmas), la
C mediana y la moda.
        CALL BUTTON(05,'mean cut',0)
        CALL BUTTON(11,'added cut',0)
        IF(NBUFF1.NE.0)THEN
          CALL BUTTON(17,'mean (<3\\gs)',0)
          CALL BUTTON(23,'median',0)
        END IF
C
ccc50      CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
        CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
        CALL IFBUTTON(XC,YC,NB)
C
        IF(NB.EQ.30)THEN
          GOTO 5
        ELSEIF(NB.EQ.42)THEN
          GOTO 900
        ELSEIF(NB.EQ.05)THEN
          CMODE='1'
        ELSEIF(NB.EQ.11)THEN
          CMODE='2'
        ELSEIF(NB.EQ.17)THEN
          CMODE='3'
        ELSEIF(NB.EQ.23)THEN
          CMODE='4'
        ELSE
          GOTO 60
        END IF
C
        CALL BUTTON(05,'mean cut',3)
        CALL BUTTON(11,'added cut',3)
        IF(NBUFF1.NE.0)THEN
          CALL BUTTON(17,'mean (<3\\gs)',3)
          CALL BUTTON(23,'median',3)
        END IF
        IF(CMODE.EQ.'1')THEN
          CALL BUTTON(05,'mean cut',1)
        ELSEIF(CMODE.EQ.'2')THEN
          CALL BUTTON(11,'added cut',1)
        ELSEIF(CMODE.EQ.'3')THEN
          CALL BUTTON(17,'mean (<3\\gs)',1)
        ELSEIF(CMODE.EQ.'4')THEN
          CALL BUTTON(23,'median',1)
        END IF

C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
        IF(NBUFF1.EQ.0)THEN            !el corte proviene de un fichero externo
          WRITE(*,100)'Name of the external file'
          NEWFILE=INFILEX(20,'@',NSCAN2,NCHAN2,STWV2,DISP2,1,.FALSE.)
C
          IF(CUTIL.EQ.'x')THEN
            CALL CCSIZE(NSCAN,NCHAN,STWV,DISP,
     +       NSCAN,NCHAN2,STWV2,DISP2,LOK)
          ELSEIF(CUTIL.EQ.'y')THEN
            CALL CCSIZE(NSCAN,NCHAN,STWV,DISP,
     +       NSCAN2,NCHAN,STWV2,DISP2,LOK)
          END IF
C
          IF(.NOT.LOK)THEN
            CLOSE(20)
            GOTO 5
          END IF
C
          IF(CUTIL.EQ.'x')THEN
            IF(NSCAN2.EQ.1)THEN
              READ(20)(XCUT(J),J=1,NCHAN)
            ELSE
              CALL PIDELIMITS(1,NSCAN2,NSMAX,'y',IFSCAN,K)
              DO J=1,NCHAN
                XCUT(J)=0.0
              END DO
              K=0
              DO I=1,NSCAN2
                READ(20)(SSX(J),J=1,NCHAN)
                IF(IFSCAN(I))THEN
                  K=K+1
                  DO J=1,NCHAN
                    XCUT(J)=XCUT(J)+SSX(J)
                  END DO
                END IF
              END DO
              IF(CMODE.EQ.'1')THEN
                DO J=1,NCHAN
                  XCUT(J)=XCUT(J)/REAL(K)
                END DO
              END IF
            END IF
          ELSEIF(CUTIL.EQ.'y')THEN
            IF(NCHAN2.EQ.1)THEN
              DO I=1,NSCAN
                READ(20)YCUT(I)
              END DO
            ELSE
              CALL PIDELIMITS(1,NCHAN2,NCMAX,'x',IFCHAN,K)
              DO I=1,NSCAN
                YCUT(I)=0.0
              END DO
              K=0
              DO J=1,NCHAN2
                DO I=1,NSCAN
                  READ(20)SSY(I)
                END DO
                IF(IFCHAN(J))THEN
                  K=K+1
                  DO I=1,NSCAN
                    YCUT(I)=YCUT(I)+SSY(I)
                  END DO
                END IF
              END DO
              IF(CMODE.EQ.'1')THEN
                DO I=1,NSCAN
                  YCUT(I)=YCUT(I)/REAL(K)
                END DO
              END IF
            END IF
          END IF
C
          CLOSE(20)
        END IF
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
        IF(NBUFF1.NE.0)THEN                     !el corte proviene de un buffer
C..............................................................................
          IF(CUTIL.EQ.'x')THEN
            IF(NSCAN.EQ.1)THEN
              DO J=1,NCHAN
                XCUT(J)=A(J,1,NBUFF1)
              END DO
            ELSE
              CALL PIDELIMITS(1,NSCAN,NSMAX,'y',IFSCAN,K)
              IF((CMODE.EQ.'1').OR.(CMODE.EQ.'2'))THEN
                DO J=1,NCHAN
                  XCUT(J)=0.0
                END DO
                K=0
                DO I=1,NSCAN
                  IF(IFSCAN(I))THEN
                    K=K+1
                    DO J=1,NCHAN
                      XCUT(J)=XCUT(J)+A(J,I,NBUFF1)
                    END DO
                  END IF
                END DO
                IF(CMODE.EQ.'1')THEN
                  DO J=1,NCHAN
                    XCUT(J)=XCUT(J)/REAL(K)
                  END DO
                END IF
              ELSE
                DO J=1,NCHAN
                  K=0
                  DO I=1,NSCAN
                    IF(IFSCAN(I))THEN
                      K=K+1
                      YCUT(K)=A(J,I,NBUFF1)
                    END IF
                  END DO
                  IF(K.GT.2)THEN
                    IF(CMODE.EQ.'3')THEN
                      FMEANSIGMA=FMEAN2(K,YCUT,3.0)
                      XCUT(J)=FMEANSIGMA
                    ELSEIF(CMODE.EQ.'4')THEN
                      FMEDIAN=FMEDIAN1(K,YCUT)
                      XCUT(J)=FMEDIAN
                    END IF
                  ELSE
                    XCUT(J)=FMEAN1(K,YCUT)
                  END IF
                END DO
              END IF
            END IF
C..............................................................................
          ELSEIF(CUTIL.EQ.'y')THEN
            IF(NCHAN.EQ.1)THEN
              DO I=1,NSCAN
                YCUT(I)=A(1,I,NBUFF1)
              END DO
            ELSE
              CALL PIDELIMITS(1,NCHAN,NCMAX,'x',IFCHAN,K)
              IF((CMODE.EQ.'1').OR.(CMODE.EQ.'2'))THEN
                DO I=1,NSCAN
                  YCUT(I)=0.0
                END DO
                K=0
                DO J=1,NCHAN
                  IF(IFCHAN(J))THEN
                    K=K+1
                    DO I=1,NSCAN
                      YCUT(I)=YCUT(I)+A(J,I,NBUFF1)
                    END DO
                  END IF
                END DO
                IF(CMODE.EQ.'1')THEN
                  DO I=1,NSCAN
                    YCUT(I)=YCUT(I)/REAL(K)
                  END DO
                END IF
              ELSE
                DO I=1,NSCAN
                  K=0
                  DO J=1,NCHAN
                    IF(IFCHAN(J))THEN
                      K=K+1
                      XCUT(K)=A(J,I,NBUFF1)
                    END IF
                  END DO
                  IF(K.GT.2)THEN
                    IF(CMODE.EQ.'3')THEN
                      FMEANSIGMA=FMEAN2(K,XCUT,3.0)
                      YCUT(I)=FMEANSIGMA
                    ELSEIF(CMODE.EQ.'4')THEN
                      FMEDIAN=FMEDIAN1(K,XCUT)
                      YCUT(I)=FMEDIAN
                    END IF
                  ELSE
                    YCUT(I)=FMEAN1(K,XCUT)
                  END IF
                END DO
              END IF
            END IF
C..............................................................................
          END IF
        END IF
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
        CALL BUTTON(36,'go',2)
        CALL BUTTON(36,'go',-3)
C
60      CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
        CALL IFBUTTON(XC,YC,NB)
C
        IF(NB.EQ.30)THEN
          GOTO 5
        ELSEIF(NB.EQ.42)THEN
          GOTO 900
        ELSE
          IF(NB.NE.36) GOTO 60
        END IF
C
        IF(CUTIL.EQ.'x')THEN
          IF(COPER.EQ.'/')THEN
            LDIVZERO=.FALSE.
            DO J=NC1,NC2
              IF(XCUT(J).EQ.0.0) LDIVZERO=.TRUE.
            END DO
            IF(LDIVZERO)THEN
              WRITE(*,101)'ERROR: division by zero '
              WRITE(*,100)'Press <CR> to continue...'
              READ(*,*)
              GOTO 5
            END IF
          END IF
          IF(COPER.EQ.'+')THEN
            DO I=NS1,NS2
              DO J=NC1,NC2
                A(J,I,NBUFF0)=A(J,I,NBUFF0)+XCUT(J)
              END DO
            END DO
          ELSEIF(COPER.EQ.'-')THEN
            DO I=NS1,NS2
              DO J=NC1,NC2
                A(J,I,NBUFF0)=A(J,I,NBUFF0)-XCUT(J)
              END DO
            END DO
          ELSEIF(COPER.EQ.'*')THEN
            DO I=NS1,NS2
              DO J=NC1,NC2
                A(J,I,NBUFF0)=A(J,I,NBUFF0)*XCUT(J)
              END DO
            END DO
          ELSEIF(COPER.EQ.'/')THEN
            DO I=NS1,NS2
              DO J=NC1,NC2
                A(J,I,NBUFF0)=A(J,I,NBUFF0)/XCUT(J)
              END DO
            END DO
          END IF
          WRITE(*,'(A,I1,A)')'>>> Buffer #',NBUFF0,
     +     ' has been modified.'
          LDEFBUFF(NBUFF0)=.TRUE.        !el buffer se activa (si no lo estaba)
          LUSEBUFF(NBUFF0)=.TRUE.
          L=TRUELEN(INFILEBUFF(NBUFF0))
          IF(L.LT.75)INFILEBUFF(NBUFF0)(L+1:L+1)='*'
          GOTO 5
        END IF
C------------------------------------------------------------------------------
        IF(CUTIL.EQ.'y')THEN
          IF(COPER.EQ.'/')THEN
            LDIVZERO=.FALSE.
            DO I=NS1,NS2
              IF(YCUT(I).EQ.0.0) LDIVZERO=.TRUE.
            END DO
            IF(LDIVZERO)THEN
              WRITE(*,101)'ERROR: division by zero '
              WRITE(*,100)'Press <CR> to continue...'
              READ(*,*)
              GOTO 5
            END IF
          END IF
          IF(COPER.EQ.'+')THEN
            DO I=NS1,NS2
              DO J=NC1,NC2
                A(J,I,NBUFF0)=A(J,I,NBUFF0)+YCUT(I)
              END DO
            END DO
          ELSEIF(COPER.EQ.'-')THEN
            DO I=NS1,NS2
              DO J=NC1,NC2
                A(J,I,NBUFF0)=A(J,I,NBUFF0)-YCUT(I)
              END DO
            END DO
          ELSEIF(COPER.EQ.'*')THEN
            DO I=NS1,NS2
              DO J=NC1,NC2
                A(J,I,NBUFF0)=A(J,I,NBUFF0)*YCUT(I)
              END DO
            END DO
          ELSEIF(COPER.EQ.'/')THEN
            DO I=NS1,NS2
              DO J=NC1,NC2
                A(J,I,NBUFF0)=A(J,I,NBUFF0)/YCUT(I)
              END DO
            END DO
          END IF
          WRITE(*,'(A,I1,A)')'>>> Buffer #',NBUFF0,
     +     ' has been modified.'
          LDEFBUFF(NBUFF0)=.TRUE.        !el buffer se activa (si no lo estaba)
          LUSEBUFF(NBUFF0)=.TRUE.
          L=TRUELEN(INFILEBUFF(NBUFF0))
          IF(L.LT.75)INFILEBUFF(NBUFF0)(L+1:L+1)='*'
          GOTO 5
        END IF
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
900     CALL PGERAS
        CALL BUTTSBR(X1B,X2B,Y1B,Y2B)
        RETURN
C
100     FORMAT(A,$)
101     FORMAT(A)
        END
C
C******************************************************************************
C
        SUBROUTINE PIDELIMITS(N1,N2,NMAX,CAXIS,IFF,K)
        IMPLICIT NONE
C
        INTEGER N1,N2,NMAX
        CHARACTER*1 CAXIS
        LOGICAL IFF(NMAX)
C
        INTEGER I,I1,I2
        INTEGER K
C------------------------------------------------------------------------------
        DO I=N1,N2
          IFF(I)=.FALSE.
        END DO
C
10      IF(CAXIS.EQ.'x')THEN
          WRITE(*,100)'Channel region (0,0=EXIT) '
        ELSE
          WRITE(*,100)'Scan region (0,0=EXIT) '
        END IF
        CALL READ2I('0,0',I1,I2)
        IF((I1.EQ.0).AND.(I2.EQ.0)) GOTO 20
        IF((I1.LT.N1).OR.(I2.GT.N2).OR.(I1.GT.I2))THEN
          WRITE(*,101)'ERROR: invalid entry. Try again.'
          GOTO 10
        END IF
        DO I=I1,I2
          IFF(I)=.TRUE.
        END DO
        GOTO 10
C
20      K=0
        DO I=N1,N2
          IF(IFF(I)) K=K+1
        END DO
        IF(K.EQ.0)THEN
          WRITE(*,100)'ERROR: number of '
          IF(CAXIS.EQ.'x')THEN
            WRITE(*,100)' channels=0!'
          ELSE
            WRITE(*,100)' scans=0!'
          END IF
          WRITE(*,101)' Try again.'
          WRITE(*,*)
          GOTO 10
        ELSE
          WRITE(*,100)'> Number of '
          IF(CAXIS.EQ.'x')THEN
            WRITE(*,100)' channels'
          ELSE
            WRITE(*,100)' scans'
          END IF
          WRITE(*,100)' to be employed: '
          WRITE(*,*)K
        END IF
C
100     FORMAT(A,$)
101     FORMAT(A)
        END
C
C******************************************************************************
C Determina si hemos indicado un subconjunto de scans despues del nombre
C del fichero. En ese caso retornamos en NS1 y NS2 el intervalo
C correspondiente. Si hay algun problema con los numeros, en NS1 y NS2
C retornamos 0. En CSEP retorna el separador de los numeros (en caso de
C que haya mas de uno)
        SUBROUTINE HAYCOMA(INFILE,NSCAN,NS1,NS2,CSEP)
        IMPLICIT NONE
        CHARACTER*(*) INFILE
        INTEGER NSCAN
        INTEGER NS1,NS2
        CHARACTER*1 CSEP
C
        INTEGER LCOMA1,LCOMA2,LMAS
C------------------------------------------------------------------------------
        CSEP=' '
        LCOMA1=INDEX(INFILE,',')
        IF(LCOMA1.EQ.0)THEN             !no hay ningun numero detras del nombre
          NS1=1
          NS2=NSCAN
        ELSE           !hay un scan (o un intervalo de scans) detras del nombre
          LCOMA2=INDEX(INFILE(LCOMA1+1:),',')
          IF(LCOMA2.EQ.0)THEN                                 !no hay dos comas
            LMAS=INDEX(INFILE(LCOMA1+1:),'+')
            IF(LMAS.EQ.0)THEN                               !solo hay un numero
              READ(INFILE(LCOMA1+1:),*,ERR=900) NS1
              NS2=NS1
            ELSE                                    !hay que sumar un intervalo
              READ(INFILE(LCOMA1+1:LCOMA1+LMAS-1),*,ERR=900) NS1
              READ(INFILE(LCOMA1+LMAS+1:),*,ERR=900) NS2
              CSEP='+'
            END IF
          ELSE                          !hay dos numeros separados por una coma
            READ(INFILE(LCOMA1+1:),*,ERR=900) NS1,NS2
            CSEP=','
          END IF
          IF((NS1.LT.1).OR.(NS2.GT.NSCAN).OR.(NS1.GT.NS2))THEN
            WRITE(*,'(A)') 'ERROR: scan number(s) out of range'
            NS1=0
            NS2=0
          END IF
        END IF
        RETURN
C
900     WRITE(*,'(A)') 'ERROR: invalid scan number(s)'
        NS1=0
        NS2=0
C
        END
