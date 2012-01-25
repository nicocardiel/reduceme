C------------------------------------------------------------------------------
C Version 13-October-2007                                         file: index.f
C @ ncl & fjg
C------------------------------------------------------------------------------
Comment
C
C Program: index
C Classification: measurement
C Description: Measures indices in a fully calibrated image. The program
C requires an external file containing the index definitions. For this purpose,
C the program looks first for a file called 'myindex.dat' in the current
C directory. If this file does not exist, the program then looks for a file
C called 'index.dat' (located in the subdirectory 'files' of the distribution
C package). If this last file is also missing, the program stops.
C
Comment
C
C Programa para medir indices en un espectro
C
C El programa devuelve una columna con un valor STATUS=ABC con tres digitos, 
C cuyo significado es el siguiente:
C digito A: 0 --> el indice se midio normalmente
C           1 --> error al medir el indice
C digito B: 0 --> error en Vradial medido sin problemas
C           1 --> error en Vradial imposible de medir en un lado
C           2 --> error en Vradial imposible de medir a ambos lados
C digito C: 0 --> error en cal. flujo medido en todas las curvas disponibles
C           1 --> idem imposible de medir en 1 de las curvas disponibles
C           2 --> idem en 2 de las curvas disponibles
C           n --> idem en n de las curvas disponibles (max. 9). Si se producen
C                 mas de 9 medidas no posibles, el valor mostrado seguira
C                 siendo 9.
C 
        PROGRAM INDEX
        IMPLICIT NONE
C
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER TRUELEN
        INTEGER READI
        INTEGER READILIM
        REAL READF
C
        INTEGER NFMAX                        !numero maximo de curvas respuesta
        PARAMETER (NFMAX=101)
        REAL C                                      !velocidad de la luz (km/s)
        PARAMETER (C=2.9979246E+5)
C
        INTEGER I,J,K,L,NL,M
        INTEGER IALPHA,IWL
        INTEGER NCREST,NINDEXT
        INTEGER ITI                    !tipo de indice (ver subrutina SELINDEX)
        INTEGER NS1,NS2
        INTEGER NSCAN1,NSCAN2
        INTEGER NCHAN2
        INTEGER NINDEX,NCRES
        INTEGER NCONTI,NABSOR
        INTEGER MIDEIND,ID,ID0
        INTEGER NIND1,NIND2
        INTEGER NSIMUL,NSEED     !numero simulaciones y semilla num. aleatorios
        INTEGER NWINX,NWINY
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        REAL A(NCMAX,NSMAX),ERR(NCMAX,NSMAX)
        REAL FLUX(NCMAX,NFMAX)
        REAL FINDEX,EINDEX,FINDEXX,EINDEXRV,EINDEXLDO,EJJGG,ESIMU
        REAL EINDEXX,EJJGGX,ESIMUX,SNX
        REAL STWV2,DISP2
        REAL FRACTION
        REAL WLMIN,WLMAX,WLMINZ,WLMAXZ,WLMIN0,WLMAX0
        REAL RVEL(NSMAX),RCVEL,RCVEL1(NSMAX),ERVEL(NSMAX)
        REAL RVELERR0,RVELERR1(NSMAX),RVELERR2(NSMAX)
        REAL RVELMIN,RVELMAX,RCVEL1MIN,RCVEL1MAX
        REAL WV(NWVMAX),FWV(NWVMAX/4)             !limites l.d.o. de las bandas
        REAL SN    !senhal/ruido por Angs. promedio en la region de los indices
        REAL ALPHAPLAW                   !valor de alfa en una ley de potencias
        REAL RESCHAN(NCMAX)    !curva de la banda para fijar la fraccion de luz
        REAL EBVPLAW            !E(B-V) para la correccion de extincion interna
        CHARACTER*1 CERR,CFCAL,COUT,CRVTABLE,CPLOT,CDUM
        CHARACTER*1 CRCS,CMIDE,CASK,CFHEAD,CFHEAD2,CPAUSA,CPLAW
        CHARACTER*1 CNEWGO,CDEFAULT,CFF,CSIMUL
        CHARACTER*1 CBANDPLAW !banda para calcular fraccion de luz en power law
        CHARACTER*8 CLABEL                               !nombre de los indices
        CHARACTER*15 SCANNAME(NSMAX)  !nombre de cada scan (del fichero RVFILE)
        CHARACTER*50 CDUMMY
        CHARACTER*75 INFILE,ERRFILE,FLUXFILE,RVFILE,FILEHEAD,OUTFILE
        CHARACTER*75 FILESIMUL
        LOGICAL LANY                                !algun indice puede medirse
        LOGICAL LINDOK(NINDMAX)   !indice medible en rango espectral disponible
        LOGICAL FIRSTPLOT
        LOGICAL IFSCAN(NSMAX)
        LOGICAL LPLOT,LPERRI,LPLINES,LPBANDS
        LOGICAL LPLOT_OLD
        LOGICAL LYMIN                            !si .TRUE. YMIN=0 en los plots
        LOGICAL LOGFILE
        LOGICAL LCOLOR(MAX_ID_RED)
        LOGICAL LRUN,LMANUAL,LHTML
C variables con informacion de cabecera en el caso de que la imagen a medir
C se haya generado como union de espectros de imagenes diferentes
        INTEGER NSCANHEAD(NSMAX),NCHANHEAD(NSMAX)
        REAL STWVHEAD(NSMAX),DISPHEAD(NSMAX)
        REAL AIRMASSHEAD(NSMAX),TIMEXPOSHEAD(NSMAX)
        CHARACTER*20 OBJECTHEAD(NSMAX),FITSFILEHEAD(NSMAX)
C variable con informacion de binning en el caso de que la imagen se haya
C creado con adnsc.f
        CHARACTER*15 CSCANBINNING(NSMAX)
C
        COMMON/BLKINDEX0/A,ERR
        COMMON/BLKINDEX1/FLUX
        COMMON/BLKINDEX2/WLMIN
        COMMON/BLKINDEX3/NSCAN,NCHAN,STWV,DISP
        COMMON/BLKINDEX4/NCREST
        COMMON/BLKINDEX6/INFILE
        COMMON/BLKINDEX8A/LPLOT,LPERRI           !dibujar o no, con/sin errores
        COMMON/BLKINDEX8B/LPLINES,LPBANDS           !dibujar o no lineas/bandas
        COMMON/BLKINDEX9/NSIMUL,NSEED
        COMMON/BLKINDEX10A/ALPHAPLAW,EBVPLAW
        COMMON/BLKINDEX10B/CBANDPLAW
        COMMON/BLKINDEX10C/RESCHAN
        COMMON/BLKINDEX11A/OBJECTHEAD,CSCANBINNING
        COMMON/BLKINDEX11B/CFHEAD,CFHEAD2,CRVTABLE
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C NOTA: si en un futuro el programa se cambiara y se permitiera volver al
C principio del programa con un GOTO, hay que considerar el hecho de que la
C semilla de los numeros aleatorios ya estaria calculada.
        THISPROGRAM='index'
        CALL WELCOME('13-October-2007')
C
        CALL CHEQUEA_FILEINDEX
C
        CALL PIDEGTER(NTERM,IDN,LCOLOR)
C
        FIRSTPLOT=.TRUE.                       !el siguiente plot es el primero
        CPLOT='y'
        CPLAW='n'
        COUT='n'
        OUTFILE=' '
        CSIMUL='n'
        ALPHAPLAW=1.0
        CBANDPLAW='V'
        EBVPLAW=0.0
C
        WRITE(*,100)'Work with error images (y/n) '
        CERR(1:1)=READC('y','yn')
C
        WRITE(*,100)'Input file name......'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        FILEHEAD=INFILE
        IF(CERR.EQ.'y') CALL GUESSEF(INFILE,ERRFILE)
        IF(TRUELEN(OBJECT).GT.0)THEN
          INFILE=INFILE(1:TRUELEN(INFILE))//' ['//
     +     OBJECT(1:TRUELEN(OBJECT))//']'
        END IF
        DO I=1,NSCAN
          READ(20) (A(J,I),J=1,NCHAN)
        END DO
        CLOSE(20)
C
        IF(CERR.EQ.'y')THEN
          WRITE(*,100)'Input error file name '
          ERRFILE=INFILEX(22,ERRFILE,NSCAN,NCHAN,STWV,DISP,21,.TRUE.)  !match
          DO I=1,NSCAN
            READ(22) (ERR(J,I),J=1,NCHAN)
          END DO
          CLOSE(22)
        ELSE
          ERRFILE='[none]'
        END IF
C
        IF(CERR.EQ.'y')THEN
          NSEED=-1
        END IF
C------------------------------------------------------------------------------
C velocidad radial y error
        WRITE(*,*)
        WRITE(*,100)'Radial velocity from table file'//
     +   '............. (y/n) '
        CRVTABLE(1:1)=READC('n','yn')
        IF(CRVTABLE.EQ.'y')THEN
          IF(CERR.EQ.'n')THEN
            WRITE(*,101)'* NOTE: File format must be the following:'//
     +       '  FORMAT(F8.1,1X,A15)'
            WRITE(*,*)
            WRITE(*,101)'      Rad.Vel. =====Object===='
            WRITE(*,101)'      ------------------------'
            WRITE(*,101)'      xxxxxx.x NameDescription'
            WRITE(*,101)'      123456789012345678901234'
            WRITE(*,101)'      000000000111111111122222'
            WRITE(*,*)
            WRITE(*,100)'Table file name'
            RVFILE=INFILEX(25,'@',0,0,0.,0.,3,.FALSE.)
            DO I=1,NSCAN
              READ(25,'(F8.1,1X,A15)') RVEL(I),SCANNAME(I)
              RCVEL=RVEL(I)/C
              RCVEL1(I)=1.+RCVEL                        !(1+z)
              RCVEL1(I)=RCVEL1(I)/SQRT(1.-RCVEL*RCVEL)  !correccion relativista
            END DO
            CLOSE(25)
            WRITE(*,*)
          ELSE
            WRITE(*,101)'* NOTE: File format must be the following:'//
     +       '  FORMAT(A15,1X,2F8.1)'
            WRITE(*,*)
            WRITE(*,101)'      =====Object====Rad.Vel.=Err.Vel='
            WRITE(*,101)'      --------------------------------'
            WRITE(*,101)'      NameDescription xxxxxx.xyyyyyy.y'
            WRITE(*,101)'      01234567890123456789012345678901'
            WRITE(*,*)
            WRITE(*,100)'Table file name'
            RVFILE=INFILEX(25,'@',0,0,0.,0.,3,.FALSE.)
            DO I=1,NSCAN
              READ(25,'(A15,1X,2F8.1)') SCANNAME(I),RVEL(I),ERVEL(i)
              RCVEL=RVEL(I)/C
              RCVEL1(I)=1.+RCVEL                        !(1+z)
              RCVEL1(I)=RCVEL1(I)/SQRT(1.-RCVEL*RCVEL)  !correccion relativista
              RCVEL=(RVEL(I)-ERVEL(i))/C
              RVELERR1(I)=1.+RCVEL
              RVELERR1(I)=RVELERR1(I)/SQRT(1.-RCVEL*RCVEL)   !corr. relativista
              RCVEL=(RVEL(I)+ERVEL(i))/C
              RVELERR2(I)=1.+RCVEL
              RVELERR2(I)=RVELERR2(I)/SQRT(1.-RCVEL*RCVEL)   !corr. relativista
            END DO
            RVELERR0=0.0
            CLOSE(25)
            WRITE(*,*)
          END IF
        ELSE
          RVFILE='[none]'
          WRITE(*,100)'Radial velocity (km/sec)'//
     +     '...............................'
          RVEL(1)=READF('@')
          RCVEL=RVEL(1)/C
          DO I=1,NSCAN
            RVEL(I)=RVEL(1)
            RCVEL1(I)=1.+RCVEL                          !(1+z)
            RCVEL1(I)=RCVEL1(I)/SQRT(1.-RCVEL*RCVEL)    !correccion relativista
          END DO
          IF(CERR.EQ.'y')THEN
            WRITE(*,100)'Radial velocity shift '//
     +       'to estimate error (km/sec) '
            RVELERR0=READF('0.0')
            DO I=1,NSCAN
              ERVEL(I)=RVELERR0
              RCVEL=(RVEL(I)-RVELERR0)/C
              RVELERR1(I)=1.+RCVEL
              RVELERR1(I)=RVELERR1(I)/SQRT(1.-RCVEL*RCVEL)   !corr. relativista
              RCVEL=(RVEL(I)+RVELERR0)/C
              RVELERR2(I)=1.+RCVEL
              RVELERR2(I)=RVELERR2(I)/SQRT(1.-RCVEL*RCVEL)   !corr. relativista
            END DO
          END IF
        END IF
C------------------------------------------------------------------------------
        WRITE(*,100)'Header information from gluesc.f (y/n) '
        CFHEAD(1:1)=READC('n','yn')
        IF(CFHEAD.EQ.'y')THEN
          WRITE(*,100)'Header information file name'
          FILEHEAD=INFILEX(27,'@',0,0,.0,.0,3,.FALSE.)
          DO I=1,NSCAN
            READ(27,'(A20,1X,A20,2(1X,I6),1X,F8.2,1X,F8.3,1X,
     +       F8.6,1X,F8.1)')
     +       OBJECTHEAD(I),FITSFILEHEAD(I),NSCANHEAD(I),NCHANHEAD(I),
     +       STWVHEAD(I),DISPHEAD(I),AIRMASSHEAD(I),TIMEXPOSHEAD(I)
          END DO
          CLOSE(27)
        ELSE
          IF(CRVTABLE.EQ.'y')THEN
            DO I=1,NSCAN
              OBJECTHEAD(I)=SCANNAME(I)
            END DO
          END IF
        END IF
C
        IF((CFHEAD.EQ.'n').AND.(CRVTABLE.EQ.'n'))THEN
          WRITE(*,*)
          WRITE(*,101)'(note: # uses the own scan number)'
          WRITE(*,100)'Header information from adnsc.f  (y/#/n) '
          CFHEAD2(1:1)=READC('n','y#n')
          IF(CFHEAD2.EQ.'y')THEN
            WRITE(*,100)'Header information file name '
            FILEHEAD=INFILEX(27,FILEHEAD(1:TRUELEN(FILEHEAD))//'.log',
     +       0,0,.0,.0,3,.FALSE.)
            DO I=1,NSCAN
              READ(27,'(A15)') CSCANBINNING(I)
            END DO
            CLOSE(27)
          ELSEIF(CFHEAD2.EQ.'#')THEN
            DO I=1,NSCAN
              WRITE(CDUMMY,'(I7,A1,I7)')I,',',I
              CALL RMBLANK(CDUMMY,CDUMMY,L)
              WRITE(CSCANBINNING(I),'(A15)') CDUMMY(1:L)
            END DO
            CFHEAD2='y'
          END IF
        ELSE
          CFHEAD2='n'
        END IF
C------------------------------------------------------------------------------
C l.d.o. del extremo izquierdo del primer pixel (STWV se refiere al centro del
C primer pixel) y del extremo derecho del ultimo pixel a Vrad=0 y a Vrad
C correspondiente a la velocidad radial introducida (si hay diferentes 
C velocidades radiales para cada scan, se determina la velocidad radial maxima
C y minima, de modo que los limites en l.d.o. se calculan en el caso mas
C exigente (el valor minimo es el mayor de todos los valores minimos y el
C valor maximo el menor de todos los valores maximos)
        WLMIN=STWV-DISP/2.0
        WLMAX=STWV+REAL(NCHAN-1)*DISP+DISP/2.0
        WRITE(*,*)
        WRITE(CDUMMY,*)WLMIN
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,100)'>>> Limits: from '//CDUMMY(1:L)
        WRITE(CDUMMY,*)WLMAX
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,101)' to '//CDUMMY(1:L)//
     +   ' angstroms (at v=0 km/sec)'
C
        IF(CRVTABLE.EQ.'y')THEN
          RCVEL1MIN=RCVEL1(1)
          RCVEL1MAX=RCVEL1(1)
          DO I=2,NSCAN
            IF(RCVEL1(I).LT.RCVEL1MIN) RCVEL1MIN=RCVEL1(I)
            IF(RCVEL1(I).GT.RCVEL1MAX) RCVEL1MAX=RCVEL1(I)
          END DO
          RVELMIN=RVEL(1)
          RVELMAX=RVEL(1)
          DO I=2,NSCAN
            IF(RVEL(I).LT.RVELMIN) RVELMIN=RVEL(I)
            IF(RVEL(I).GT.RVELMAX) RVELMAX=RVEL(I)
          END DO
        ELSE
          RCVEL1MIN=RCVEL1(1)
          RCVEL1MAX=RCVEL1(1)
          RVELMIN=RVEL(1)
          RVELMAX=RVEL(1)
        END IF
C
        WLMINZ=WLMIN/RCVEL1MIN                             !l.d.o. minima mayor
        WLMAXZ=WLMAX/RCVEL1MAX                             !l.d.o. maxima menor
        IF(WLMINZ.GE.WLMAXZ)THEN
          WRITE(*,101)'FATAL ERROR: WLMINZ.eq.WLMAXZ'
          STOP
        END IF
        WRITE(CDUMMY,*)WLMINZ
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,100)'>>> Limits: from '//CDUMMY(1:L)
        WRITE(CDUMMY,*)WLMAXZ
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,101)' to '//CDUMMY(1:L)//
     +   ' angstroms (at v=Vrad km/sec)'
C------------------------------------------------------------------------------
C determinamos que indices pueden ser medidos (ver subrutina SELINDEX)
        LANY=.FALSE.
        CALL SELINDEX(0,WV,FWV,NINDEXT,CLABEL)    !numero de indices disponible
        DO K=1,NINDEXT
          CALL SELINDEX(K,WV,FWV,ITI,CLABEL)
          IF((ITI.GT.-100).AND.(ITI.LT.-1))THEN
            NCONTI=-ITI
            NABSOR=0
            WLMIN0=WV(1)
            DO IWL=2,2*NCONTI
              WLMIN0=AMIN1(WLMIN0,WV(IWL))
            END DO
            WLMAX0=WV(1)
            DO IWL=2,2*NCONTI
              WLMAX0=AMAX1(WLMAX0,WV(IWL))
            END DO
            LINDOK(K)=((WLMIN0.GE.WLMINZ).AND.(WLMAX0.LE.WLMAXZ))
          ELSEIF((ITI.EQ.1).OR.(ITI.EQ.2))THEN
            WLMIN0=WV(1)
            DO IWL=2,6
              WLMIN0=AMIN1(WLMIN0,WV(IWL))
            END DO
            WLMAX0=WV(1)
            DO IWL=2,6
              WLMAX0=AMAX1(WLMAX0,WV(IWL))
            END DO
            LINDOK(K)=((WLMIN0.GE.WLMINZ).AND.(WLMAX0.LE.WLMAXZ))
          ELSEIF((ITI.EQ.3).OR.(ITI.EQ.4).OR.(ITI.EQ.5))THEN
            LINDOK(K)=((WV(1).GE.WLMINZ).AND.(WV(4).LE.WLMAXZ))
          ELSEIF((ITI.GE.101).AND.(ITI.LE.9999))THEN
            NCONTI=(ITI/100)                    !numero de regiones de continuo
            NABSOR=ITI-NCONTI*100           !numero de regiones con absorciones
            WLMIN0=WV(1)
            DO IWL=2,2*(NCONTI+NABSOR)
              WLMIN0=AMIN1(WLMIN0,WV(IWL))
            END DO
            WLMAX0=WV(1)
            DO IWL=2,2*(NCONTI+NABSOR)
              WLMAX0=AMAX1(WLMAX0,WV(IWL))
            END DO
            LINDOK(K)=((WLMIN0.GE.WLMINZ).AND.(WLMAX0.LE.WLMAXZ))
          ELSE
            WRITE(*,100) 'ITI='
            WRITE(*,*) ITI
            WRITE(*,101) 'FATAL ERROR: invalid ITI value.'
          END IF
          IF(LINDOK(K))THEN
            LANY=.TRUE.
          END IF
        END DO
        IF(.NOT.LANY)THEN
          WRITE(*,*)
          WRITE(*,101)'WARNING: No. of indices which can be measured=0'
          STOP
        END IF
C------------------------------------------------------------------------------
C curva respuesta
        WRITE(*,*)
        WRITE(*,100)'Flux calibration.................. (y/n) '
        CFCAL(1:1)=READC('y','yn')
        IF(CFCAL.EQ.'y')THEN
          WRITE(*,100)'Response curve file name '
          FLUXFILE=INFILEX(23,'../stand/curvresf',NCREST,NCHAN2,
     +     STWV2,DISP2,1,.FALSE.)
          IF(NCHAN2.NE.NCHAN)THEN
            WRITE(*,101)'FATAL ERROR: NCHAN in this file is different.'
            CLOSE(23)
            STOP
          END IF
          IF(STWV2.NE.STWV)THEN
            WRITE(*,101)'FATAL ERROR: STWV in this file is different.'
            CLOSE(23)
            STOP
          END IF
          IF(DISP2.NE.DISP)THEN
            WRITE(*,101)'FATAL ERROR: DISP in this file is different.'
            CLOSE(23)
            STOP
          END IF
          IF(NCREST.GT.1)THEN
            WRITE(*,*)
            WRITE(*,101)'NOTE: this file contains more than 1 spectrum.'
            WRITE(*,101)'- The 1st spectrum will be employed as '//
     +       'response curve.'
            WRITE(*,100)'- Are you using the rest of the spectra to '//
     +       'estimate flux errors (y/n) '
            CRCS(1:1)=READC('y','yn')
            IF(CRCS.EQ.'n') NCREST=1
          END IF
          IF(NCREST.GT.NFMAX)THEN
            WRITE(*,101)'FATAL ERROR: number of response curves '//
     +       'too large. Redim NFMAX.'
            CLOSE(23)
            STOP
          END IF
          DO I=1,NCREST
            READ(23) (FLUX(J,I),J=1,NCHAN)
          END DO
          CLOSE(23)
        ELSE
          FLUXFILE='[none]'
          NCREST=1              !es igual a 1 para saltar proteccion en MIDEIND
        END IF
C------------------------------------------------------------------------------
C dividimos por curva respuesta
        IF(CFCAL.EQ.'y')THEN
          DO I=1,NSCAN
            DO J=1,NCHAN
              A(J,I)=A(J,I)/FLUX(J,1)
            END DO
          END DO
          IF(CERR.EQ.'y')THEN
            DO I=1,NSCAN
              DO J=1,NCHAN
                ERR(J,I)=ERR(J,I)/FLUX(J,1)
              END DO
            END DO
          END IF
        END IF
C------------------------------------------------------------------------------
C elegimos los indices a medir
10      WRITE(*,*)
        CALL SHINDEX(LINDOK,0)
12      WRITE(*,100)'Index number '
        NINDEX=READILIM('@',-1,NINDEXT)
        IF(NINDEX.EQ.-1) THEN
          IF(.NOT.FIRSTPLOT) CALL PGEND
          STOP
        END IF
        IF(NINDEX.EQ.0)THEN
          NIND1=1
          NIND2=NINDEXT
        ELSE
          IF(.NOT.LINDOK(NINDEX))THEN
            WRITE(*,101)'ERROR: this index is not available. Try again.'
            GOTO 12
          END IF
          NIND1=NINDEX
          NIND2=NINDEX
        END IF
C------------------------------------------------------------------------------
C decidimos pintar o no los indices mientras se miden
        WRITE(*,100)'Plot indices....................... (y/n) '
        CPLOT(1:1)=READC(CPLOT,'yn')
C
        IF(CPLOT.EQ.'y')THEN
          IF(.NOT.FIRSTPLOT)THEN
            WRITE(*,100)'Change graphic options............. (y/n) '
            CNEWGO(1:1)=READC('n','yn')
            IF(CNEWGO.EQ.'y') FIRSTPLOT=.TRUE.
          END IF
C
          IF(FIRSTPLOT)THEN
            WRITE(*,101)'Default options:'
            WRITE(*,101)'* Plot indices in all the spectra.... [y]'
            IF(CERR.EQ.'y')THEN
              WRITE(*,101)'* Plot error bars in spectrum data... [y]'
            END IF
            WRITE(*,101)'* Pause between plots................ [y]'
            WRITE(*,101)'* YMIN=0............................. [n]'
            WRITE(*,101)'* No. windows in X .................. [1]'
            WRITE(*,101)'* No. windows in Y .................. [1]'
            WRITE(*,101)'* Plot location of typical lines..... [n]'
            WRITE(*,101)'* Plot photometric bands............. [n]'
C
            WRITE(*,100)'Accept default options (y/n) '
            CDEFAULT(1:1)=READC('y','ny')
            IF(CDEFAULT.EQ.'y')THEN
              DO I=1,NSCAN
                IFSCAN(I)=.TRUE.
              END DO
              IF(CERR.EQ.'y')THEN
                LPERRI=.TRUE.
              ELSE
                LPERRI=.FALSE.
              END IF
              CASK='y'
              LYMIN=.FALSE.
              NWINX=1
              NWINY=1
              LPLINES=.FALSE.
              LPBANDS=.FALSE.
            ELSE
              DO I=1,NSCAN
                IFSCAN(I)=.FALSE.                    !scans que vamos a dibujar
              END DO
              WRITE(*,100)'Plot indices in all the spectra.... (y/n) '
              CMIDE(1:1)=READC('y','yn')
              IF(CMIDE.EQ.'y')THEN
                DO I=1,NSCAN
                  IFSCAN(I)=.TRUE.
                END DO
              ELSE
20              WRITE(*,100)'Scan region to be plotted (0,0=EXIT) '
                CALL READ2I('0,0',NS1,NS2)
                IF((NS1.EQ.0).AND.(NS2.EQ.0)) GOTO 21
                IF((NS1.LT.1).OR.(NS1.GT.NSCAN).OR.(NS1.GT.NS2))THEN
                  WRITE(*,101)'ERROR: invalid entry. Try again.'
                  GOTO 20
                END IF
                DO I=NS1,NS2
                  IFSCAN(I)=.TRUE.
                END DO
                GOTO 20
21              CONTINUE
              END IF
              IF(CERR.EQ.'y')THEN
                WRITE(*,100)'Plot error bars in spectrum data... (y/n) '
                CDUM(1:1)=READC('y','yn')
                LPERRI=(CDUM.EQ.'y')
              ELSE
                LPERRI=.FALSE.
              END IF
              WRITE(*,100)'Pause between plots................ (y/n) '
              CASK(1:1)=READC('y','yn')
              WRITE(*,100)'YMIN=0............................. (y/n) '
              CDUM(1:1)=READC('n','yn')
              LYMIN=(CDUM.EQ.'y')
              WRITE(*,100)'No. windows in X '
              NWINX=READI('1')
              WRITE(*,100)'No. windows in Y '
              NWINY=READI('1')
              WRITE(*,100)'Plot location of typical lines..... (y/n) '
              CDUM(1:1)=READC('n','yn')
              LPLINES=(CDUM.EQ.'y')
              WRITE(*,100)'Plot photometric bands............. (y/n) '
              CDUM(1:1)=READC('n','yn')
              LPBANDS=(CDUM.EQ.'y')
            END IF
          END IF
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGSUBP(NWINX,NWINY)
          END DO
          FIRSTPLOT=.FALSE.
        END IF
C------------------------------------------------------------------------------
C Numero de simulaciones para el calculo de errores (si procede)
        IF(CERR.EQ.'y')THEN
          WRITE(*,100)'No. of simulations for error estimation '
          NSIMUL=READI('0')
ccc       NSIMUL=0
        ELSE
          NSIMUL=0
        END IF
C------------------------------------------------------------------------------
C podemos calcular el efecto en los indices al introducir una ley de potencias
        WRITE(*,100)'Add variable power law to the data. (y/n) '
        CPLAW(1:1)=READC(CPLAW,'yn')
        IF(CPLAW.EQ.'y')THEN
          WRITE(*,*)
          WRITE(*,101)'Power law:'
          WRITE(*,101)'     Flux(lambda) = k * (lambda)**(alpha-2.0)'
          WRITE(*,100)'alpha '
          WRITE(CDUMMY,*) ALPHAPLAW
          ALPHAPLAW=READF(CDUMMY)
          WRITE(*,*)
          CALL BANDFRACTION('U',0.5*(RCVEL1MIN+RCVEL1MAX),FRACTION)
          WRITE(*,100)'> U band sampled (%): '
          WRITE(*,'(F5.1)') FRACTION*100.
          CALL BANDFRACTION('B',0.5*(RCVEL1MIN+RCVEL1MAX),FRACTION)
          WRITE(*,100)'> B band sampled (%): '
          WRITE(*,'(F5.1)') FRACTION*100.
          CALL BANDFRACTION('V',0.5*(RCVEL1MIN+RCVEL1MAX),FRACTION)
          WRITE(*,100)'> V band sampled (%): '
          WRITE(*,'(F5.1)') FRACTION*100.
          WRITE(*,*)
          FRACTION=0.0
          DO WHILE(FRACTION.LE.0.)
            WRITE(*,100)'Band to fix fraction of the power law '//
     +       'light (u/b/v) '
            CBANDPLAW(1:1)=READC(CBANDPLAW,'uUbBvV')
            CALL CHUPPER(CBANDPLAW)
C volvemos a llamar a la subrutina para calcular de nuevo la curva respuesta
C de la banda RESCHAN que sera utilizada para medir la fraccion de luz
            CALL BANDFRACTION(CBANDPLAW,0.5*(RCVEL1MIN+RCVEL1MAX),
     +       FRACTION)
            IF(FRACTION.EQ.0.)THEN
              WRITE(*,100)'ERROR: this band is outside spectrum limits!'
            END IF
          END DO
          WRITE(*,100)'E(B-V) internal extinction '
          WRITE(CDUMMY,*) EBVPLAW
          EBVPLAW=READF(CDUMMY)
        END IF
C------------------------------------------------------------------------------
C scans a medir
        IF(NSCAN.GT.1)THEN
          WRITE(*,100)'First scan to be measured '
          NSCAN1=READILIM('1',1,NSCAN)
          WRITE(*,100)'Last scan to be measured '
          WRITE(CDUMMY,*)NSCAN
          NSCAN2=READILIM(CDUMMY,NSCAN1,NSCAN)
        ELSE
          NSCAN1=1
          NSCAN2=1
        END IF
C------------------------------------------------------------------------------
C crear fichero de salida e informacion general
        WRITE(*,100)'Create output file name............ (y/n) '
        COUT(1:1)=READC(COUT,'yn')
        IF(COUT.EQ.'y')THEN
          IF(TRUELEN(OUTFILE).GT.0)THEN
            WRITE(*,100)'Output file name '
            OUTFILE(1:75)=READC(OUTFILE,'@')
          ELSE
            WRITE(*,100)'Output file name'
            OUTFILE(1:75)=READC('@','@')
          END IF
          INQUIRE(FILE=OUTFILE,EXIST=LOGFILE)
          IF(LOGFILE)THEN
            WRITE(*,*)
            WRITE(*,101)'WARNING: This file already exist.'
            WRITE(*,100)'Insert <FF> (y/n) '
            CFF(1:1)=READC('n','yn')
            WRITE(*,*)
            OPEN(30,FILE=OUTFILE,STATUS='OLD',FORM='FORMATTED')
777         READ(30,*,END=778)
            GOTO 777
778         IF(CFF.EQ.'n')THEN
              WRITE(30,*)
              WRITE(30,*)
              WRITE(30,*)
            ELSE
              WRITE(30,100) CHAR(12)
            END IF
          ELSE
            OPEN(30,FILE=OUTFILE,STATUS='NEW',FORM='FORMATTED')
          END IF
          WRITE(30,101)'--------------------------------------------'//
     +     '-----------------------------------'
          WRITE(30,101)'This file is..........................: '//
     +     OUTFILE(1:TRUELEN(OUTFILE))
          WRITE(30,101)'Input file name.......................: '//
     +     INFILE(1:TRUELEN(INFILE))
          WRITE(30,101)'Input error file name.................: '//
     +     ERRFILE(1:TRUELEN(ERRFILE))
          WRITE(30,101)'Radial velocity file name.............: '//
     +     RVFILE(1:TRUELEN(RVFILE))
          IF(CRVTABLE.EQ.'y')THEN
            WRITE(30,100)'Minimum radial velocity value (km/sec): '
            WRITE(30,'(F7.0)')RVELMIN
            WRITE(30,100)'Maximum radial velocity value (km/sec): '
            WRITE(30,'(F7.0)')RVELMAX
          ELSE
            WRITE(30,100)'Radial velocity value........ (km/sec): '
            WRITE(30,'(F7.0)')RVELMIN
          END IF
          IF(RVELERR0.NE.0.)THEN
            WRITE(30,100)'Radial velocity shift........ (km/sec): '
            WRITE(30,'(F7.0)')RVELERR0
          END IF
          WRITE(30,101)'Response curve file name..............: '//
     +     FLUXFILE(1:TRUELEN(FLUXFILE))
          IF(NCREST.GT.1)THEN
            WRITE(30,100)'No. of response curves (flux error)...: '
            WRITE(30,*)NCREST-1
          END IF
          IF(NSIMUL.GT.0)THEN
            WRITE(30,100)'No. of simulations....................: '
            WRITE(30,*)NSIMUL
          END IF
          IF(CFHEAD.EQ.'y')THEN
            WRITE(30,101)'Header information file name (gluesc).: '//
     +       FILEHEAD(1:TRUELEN(FILEHEAD))
          ELSEIF(CFHEAD2.EQ.'y')THEN
            WRITE(30,101)'Header information file name (adnsc)..: '//
     +       FILEHEAD(1:TRUELEN(FILEHEAD))
          END IF
          IF(CPLAW.EQ.'y')THEN
            WRITE(30,100)'Power law applied with an alpha value.: '
            WRITE(30,*)ALPHAPLAW
            WRITE(30,100)'E(B-V) internal extinct. (power law)..: '
            WRITE(30,*)EBVPLAW
          END IF
          WRITE(30,101)'--------------------------------------------'//
     +     '-----------------------------------'
        END IF
C------------------------------------------------------------------------------
C informacion general
        WRITE(*,*)
        WRITE(*,101)'--------------------------------------------'//
     +   '-----------------------------------'
        WRITE(*,101)'Input file name.......................: '//
     +   INFILE(1:TRUELEN(INFILE))
        WRITE(*,101)'Input error file name.................: '//
     +   ERRFILE(1:TRUELEN(ERRFILE))
        WRITE(*,101)'Radial velocity file name.............: '//
     +   RVFILE(1:TRUELEN(RVFILE))
        IF(CRVTABLE.EQ.'y')THEN
          WRITE(*,100)'Minimum radial velocity value (km/sec): '
          WRITE(*,'(F7.0)')RVELMIN
          WRITE(*,100)'Maximum radial velocity value (km/sec): '
          WRITE(*,'(F7.0)')RVELMAX
        ELSE
          WRITE(*,100)'Radial velocity value........ (km/sec): '
          WRITE(*,'(F7.0)')RVELMIN
        END IF
        IF(RVELERR0.NE.0.)THEN
          WRITE(*,100)'Radial velocity shift........ (km/sec): '
          WRITE(*,'(F7.0)')RVELERR0
        END IF
        WRITE(*,101)'Response curve file name..............: '//
     +   FLUXFILE(1:TRUELEN(FLUXFILE))
        IF(NCREST.GT.1)THEN
          WRITE(*,100)'No. of response curves (flux error)...: '
          WRITE(*,*)NCREST-1
        END IF
        IF(NSIMUL.GT.0)THEN
          WRITE(*,100)'No. of simulations....................: '
          WRITE(*,*)NSIMUL
        END IF
        IF(CFHEAD.EQ.'y')THEN
          WRITE(*,101)'Header information file name (gluesc).: '//
     +     FILEHEAD(1:TRUELEN(FILEHEAD))
        ELSEIF(CFHEAD2.EQ.'y')THEN
          WRITE(*,101)'Header information file name (adnsc)..: '//
     +     FILEHEAD(1:TRUELEN(FILEHEAD))
        END IF
        IF(CPLAW.EQ.'y')THEN
          WRITE(*,100)'Power law applied with an alpha value.: '
          WRITE(*,*)ALPHAPLAW
          WRITE(*,100)'E(B-V) internal extinct. (power law)..: '
          WRITE(*,*)EBVPLAW
        END IF
        WRITE(*,101)'--------------------------------------------'//
     +   '-----------------------------------'
C------------------------------------------------------------------------------
C crear fichero de salida e informacion general
        IF(NSIMUL.GT.0)THEN
          WRITE(*,100)'Create output file name with simul. (y/n) '
          CSIMUL(1:1)=READC(CSIMUL,'yn')
          IF(CSIMUL.EQ.'y')THEN
            WRITE(*,100)'Output file name'
            FILESIMUL=OUTFILEX(44,'@',0,0,0.,0.,3,.FALSE.)
          END IF
        ELSE
          CSIMUL='n'
        END IF
C------------------------------------------------------------------------------
C bucle en indices y numero de scans
        DO K=NIND1,NIND2                                      !bucle en indices
          IF(LINDOK(K))THEN                    !solo si el indice puede medirse
C..............................................................................
            CALL SELINDEX(K,WV,FWV,ITI,CLABEL)
            WRITE(*,*)
            WRITE(*,'(A11,I2,A10,$)') '>>> Index #',K,': '//CLABEL
            IF((ITI.GT.-100).AND.(ITI.LT.-1))THEN
              WRITE(*,100) ' (slope)     '
            ELSEIF(ITI.EQ.1)THEN
              WRITE(*,100) ' (molecular) '
            ELSEIF(ITI.EQ.2)THEN
              WRITE(*,100) ' (atomic)    '
            ELSEIF(ITI.EQ.3)THEN
              WRITE(*,100) '             '
            ELSEIF(ITI.EQ.4)THEN
              WRITE(*,100) '             '
            ELSEIF((ITI.GE.101).AND.(ITI.LE.9999))THEN
              WRITE(*,100) ' (generic)   '
            END IF
            IF(COUT.EQ.'y')THEN
              WRITE(30,*)
              WRITE(30,'(A11,I2,A10,$)') '>>> Index #',K,': '//CLABEL
              IF((ITI.GT.-100).AND.(ITI.LT.-1))THEN
                WRITE(30,100) ' (slope)     '
              ELSEIF(ITI.EQ.1)THEN
                WRITE(30,100) ' (molecular) '
              ELSEIF(ITI.EQ.2)THEN
                WRITE(30,100) ' (atomic)    '
              ELSEIF(ITI.EQ.3)THEN
                WRITE(30,100) '             '
              ELSEIF(ITI.EQ.4)THEN
                WRITE(30,100) '             '
              ELSEIF((ITI.GE.101).AND.(ITI.LE.9999))THEN
                WRITE(30,100) ' (generic)   '
              END IF
            END IF
            IF(CSIMUL.EQ.'y')THEN
              WRITE(44,*)
              WRITE(44,'(A11,I2,A10,$)') '>>> Index #',K,': '//CLABEL
              IF((ITI.GT.-100).AND.(ITI.LT.-1))THEN
                WRITE(44,100) ' (slope)     '
              ELSEIF(ITI.EQ.1)THEN
                WRITE(44,100) ' (molecular) '
              ELSEIF(ITI.EQ.2)THEN
                WRITE(44,100) ' (atomic)    '
              ELSEIF(ITI.EQ.3)THEN
                WRITE(44,100) '             '
              ELSEIF(ITI.EQ.4)THEN
                WRITE(44,100) '             '
              ELSEIF((ITI.GE.101).AND.(ITI.LE.9999))THEN
                WRITE(44,100) ' (generic)   '
              END IF
            END IF
            L=TRUELEN(INFILE)
            IF(L+42.LT.79)THEN
              DO M=1,37-L
                WRITE(*,100)' '
                IF(COUT.EQ.'y') WRITE(30,100)' '
                IF(CSIMUL.EQ.'y') WRITE(44,100)' '
              END DO
            ELSE
              WRITE(*,*)
              IF(COUT.EQ.'y') WRITE(30,*)
              IF(CSIMUL.EQ.'y') WRITE(44,*)
            END IF
            WRITE(*,101)'File: '//INFILE(1:L)
            WRITE(*,101)'--------------------------------------------'//
     +       '-----------------------------------'
            WRITE(*,100)'Scan#   INDEX    Efoto    Ervel  '//
     +       '  Eflux  (S/N)/A   Rvel ###'
            IF(CRVTABLE.EQ.'y')THEN
              WRITE(*,100)' Object'
              IF(CFHEAD.EQ.'y')THEN
                WRITE(*,101)'          Object              '//
     +           ' Fitsfile            '//
     +           '  NSCAN  NCHAN   STWV      DISP'//
     +           '  AIRMASS  TIMEXPOS'
              ELSE
                WRITE(*,101)
              END IF
            ELSE
              IF(CFHEAD.EQ.'y')THEN
                WRITE(*,101)'                 Object              '//
     +           ' Fitsfile            '//
     +           '  NSCAN  NCHAN   STWV      DISP'//
     +           '  AIRMASS  TIMEXPOS'
              ELSEIF(CFHEAD2.EQ.'y')THEN
                WRITE(*,101)' Binning (Scans)'
              ELSE
                WRITE(*,101)
              END IF
            END IF
            IF(COUT.EQ.'y')THEN
              WRITE(30,101)'File: '//INFILE(1:L)
              WRITE(30,101)'----------------------------------------'//
     +         '---------------------------------------'
              WRITE(30,100)'Scan#   INDEX    Efoto    Ervel  '//
     +         '  Eflux  (S/N)/A   Rvel ###'
              IF(CRVTABLE.EQ.'y')THEN
                WRITE(30,100)' Object'
                IF(CFHEAD.EQ.'y')THEN
                  WRITE(30,101)'          Object              '//
     +             ' Fitsfile            '//
     +             '  NSCAN  NCHAN   STWV      DISP'//
     +             '  AIRMASS  TIMEXPOS'
                ELSE
                  WRITE(30,101)
                END IF
              ELSE
                IF(CFHEAD.EQ.'y')THEN
                  WRITE(30,101)'                 Object              '//
     +             ' Fitsfile            '//
     +             '  NSCAN  NCHAN   STWV      DISP'//
     +             '  AIRMASS  TIMEXPOS'
                ELSEIF(CFHEAD2.EQ.'y')THEN
                  WRITE(30,101)' Binning (Scans)'
                ELSE
                  WRITE(30,101)
                END IF
              END IF
            END IF
            IF(CSIMUL.EQ.'y')THEN
              WRITE(44,101)'File: '//INFILE(1:L)
              WRITE(44,101)'-----------------------------------------'//
     +         '--------------------------------------'
              WRITE(44,101)'#1) Scan, #2) Index, #3) Efoto, '//
     +         '#4) Ejjgg, #5) Esimu, #6) S/N-ratio/A'
            END IF
C..............................................................................
            DO I=NSCAN1,NSCAN2                         !bucle en numero de scan
              IF(CPLOT.EQ.'y')THEN
                LPLOT=IFSCAN(I)
              ELSE
                LPLOT=.FALSE.
              END IF
C.............
              RCVEL=RCVEL1(I)
              NCRES=1           !usamos la curva respuesta promedio (o ninguna)
              IF(LPLOT)THEN
                DO ITERM=NTERM,1,-1
                  CALL PGSLCT(IDN(ITERM))
                  CALL PLOTINDEX(K,I,RCVEL,LYMIN,LCOLOR(ITERM))
                  IF(LCOLOR(ITERM)) CALL PGSCI(6)
                END DO
              END IF
              ID0=MIDEIND(I,I,ITI,WV,FWV,CERR,RCVEL,NCRES,FINDEX,EINDEX,
     +         EJJGG,ESIMU,SN,-1,.FALSE.)
              IF(LPLOT)THEN
                DO ITERM=NTERM,1,-1
                  CALL PGSLCT(IDN(ITERM))
                  IF(LCOLOR(ITERM)) CALL PGSCI(1)
                END DO
              END IF
              IF(CSIMUL.EQ.'y')THEN
                WRITE(44,*)I,FINDEX,EINDEX,EJJGG,ESIMU,SN
              END IF
C.............
              IF(CERR.EQ.'y')THEN
                EINDEXRV=0.          !error debido al shift en velocidad radial
                IF(ERVEL(I).NE.0.)THEN
                  NL=0
                  IF(LPLOT)THEN
                    DO ITERM=NTERM,1,-1
                      CALL PGSLCT(IDN(ITERM))
                      CALL PGSLS(2)
                      IF(LCOLOR(ITERM)) CALL PGSCI(4)
                    END DO
                  END IF
                  ID=MIDEIND(I,I,ITI,WV,FWV,CERR,RVELERR1(I),NCRES,
     +             FINDEXX,EINDEXX,EJJGGX,ESIMUX,SNX,-1,.FALSE.)
                  IF(ID.EQ.0)THEN
                    EINDEXRV=EINDEXRV+(FINDEX-FINDEXX)*(FINDEX-FINDEXX)
                    NL=NL+1
                  ELSE
                    ID0=ID0+10
                  END IF
                  IF(LPLOT)THEN
                    DO ITERM=NTERM,1,-1
                      CALL PGSLCT(IDN(ITERM))
                      IF(LCOLOR(ITERM)) CALL PGSCI(2)
                    END DO
                  END IF
                  ID=MIDEIND(I,I,ITI,WV,FWV,CERR,RVELERR2(I),NCRES,
     +             FINDEXX,EINDEXX,EJJGGX,ESIMUX,SNX,-1,.FALSE.)
                  IF(ID.EQ.0)THEN
                    EINDEXRV=EINDEXRV+(FINDEX-FINDEXX)*(FINDEX-FINDEXX)
                    NL=NL+1
                  ELSE
                    ID0=ID0+10
                  END IF
                  IF(NL.EQ.1)THEN
                    EINDEXRV=SQRT(EINDEXRV)
                  ELSEIF(NL.EQ.2)THEN
                    EINDEXRV=SQRT(EINDEXRV/REAL(NL-1))
                  ELSE
                    EINDEXRV=0.     !redundante, pero lo dejamos (por claridad)
                  END IF
                  IF(LPLOT)THEN
                    DO ITERM=NTERM,1,-1
                      CALL PGSLCT(IDN(ITERM))
                      CALL PGSLS(1)
                      IF(LCOLOR(ITERM)) CALL PGSCI(1)
                    END DO
                  END IF
                END IF
              ELSE
                EINDEXRV=0.
              END IF
C.............
              IF(NCREST.GT.1)THEN
                EINDEXLDO=0.                        !error calibracion en flujo
                NL=0
                DO NCRES=2,NCREST
                  ID=MIDEIND(I,I,ITI,WV,FWV,CERR,RCVEL,NCRES,
     +             FINDEXX,EINDEXX,EJJGGX,ESIMUX,SNX,-1,.FALSE.)
                  IF(ID.EQ.0)THEN
                    EINDEXLDO=EINDEXLDO+
     +               (FINDEX-FINDEXX)*(FINDEX-FINDEXX)
                    NL=NL+1
                  ELSE
                    IF(ID0.LT.900) ID0=ID0+100
                  END IF
                END DO
                IF(NL.GT.1)THEN
                  EINDEXLDO=SQRT(EINDEXLDO/REAL(NL-1))
                ELSEIF(NL.EQ.1)THEN
                  EINDEXLDO=SQRT(EINDEXLDO)
                ELSE
                  EINDEXLDO=0.       !redundante pero lo dejamos (por claridad)
                END IF
              ELSE
                EINDEXLDO=0.
              END IF
C.............
              IF(CFHEAD2.EQ.'y')THEN
                WRITE(*,'(I4,4(1X,F8.3),1X,F6.1,1X,F8.1,1X,I3.3,
     +           1X,A15,$)')I,FINDEX,EINDEX,EINDEXRV,
     +           EINDEXLDO,SN,RVEL(I),ID0,CSCANBINNING(I)
                IF(COUT.EQ.'y')THEN
                  WRITE(30,'(I4,4(1X,F8.3),1X,F6.1,1X,F8.1,1X,I3.3,
     +             1X,A15,$)')I,FINDEX,EINDEX,EINDEXRV,
     +             EINDEXLDO,SN,RVEL(I),ID0,CSCANBINNING(I)
                END IF
              ELSE
                IF(CRVTABLE.EQ.'n') SCANNAME(I)='               '
                WRITE(*,'(I4,4(1X,F8.3),1X,F6.1,1X,F8.1,1X,I3.3,
     +           1X,A15,$)')I,FINDEX,EINDEX,EINDEXRV,
     +           EINDEXLDO,SN,RVEL(I),ID0,SCANNAME(I)
                IF(COUT.EQ.'y')THEN
                  WRITE(30,'(I4,4(1X,F8.3),1X,F6.1,1X,F8.1,1X,I3.3,
     +             1X,A15,$)')I,FINDEX,EINDEX,EINDEXRV,
     +             EINDEXLDO,SN,RVEL(I),ID0,SCANNAME(I)
                END IF
              END IF
C.............
              IF(CFHEAD.EQ.'y')THEN
                WRITE(*,'(1X,A20,1X,A20,2(1X,I6),1X,F8.2,1X,F8.3,1X,
     +           F8.6,1X,F8.1,$)') OBJECTHEAD(I),FITSFILEHEAD(I),
     +           NSCANHEAD(I),NCHANHEAD(I),STWVHEAD(I),DISPHEAD(I),
     +           AIRMASSHEAD(I),TIMEXPOSHEAD(I)
                IF(COUT.EQ.'y')THEN
                  WRITE(30,'(1X,A20,1X,A20,2(1X,I6),1X,F8.2,1X,F8.3,1X,
     +             F8.6,1X,F8.1)') OBJECTHEAD(I),FITSFILEHEAD(I),
     +             NSCANHEAD(I),NCHANHEAD(I),STWVHEAD(I),DISPHEAD(I),
     +             AIRMASSHEAD(I),TIMEXPOSHEAD(I)
                END IF
              ELSE
                IF(COUT.EQ.'y') WRITE(30,*)
              END IF
C.............
              IF(CPLAW.EQ.'y')THEN
                LPLOT_OLD=LPLOT
                LPLOT=.FALSE.
                NCRES=1         !usamos la curva respuesta promedio (o ninguna)
                WRITE(*,*)
                IF(COUT.EQ.'y') WRITE(30,*)
                DO IALPHA=0,10        !estimamos efecto de una ley de potencias
                  ID0=MIDEIND(I,I,ITI,WV,FWV,CERR,RCVEL,NCRES,FINDEX,
     +             EINDEX,EJJGG,ESIMU,SN,IALPHA,.FALSE.)
                  WRITE(*,'(A1,F3.1,2(1X,F8.3))')
     +             '>',REAL(IALPHA)/10.,FINDEX,EINDEX
                  IF(COUT.EQ.'y')THEN
                    WRITE(30,'(A1,F3.1,2(1X,F8.3))')
     +               '>',REAL(IALPHA)/10.,FINDEX,EINDEX
                  END IF
                END DO
                CALL INDEXPLAW(K,ALPHAPLAW,FINDEX)     !indice de una Power-Law
                WRITE(*,'(A5,F8.3,$)')'PLaw>',FINDEX
                IF(COUT.EQ.'y')THEN
                  WRITE(30,'(A4,F8.3,$)')'PLaw>',FINDEX
                END IF
                LPLOT=LPLOT_OLD
              END IF
C.............
              IF((LPLOT).AND.(CASK.EQ.'y'))THEN
                READ(*,'(A)')CPAUSA
                CALL LRUNX(LRUN,LMANUAL,LHTML)
                IF((LRUN).OR.(LMANUAL).OR.(LHTML)) WRITE(*,*)
                IF(CPAUSA.EQ.'q')THEN
                  IF(COUT.EQ.'y') CLOSE(30)
                  GOTO 10
                END IF
              ELSE
                WRITE(*,*)
              END IF
C..............................................................................
            END DO                                              !numero de scan
          END IF                                                     !LINDOK(K)
        END DO                                                !bucle en indices
C------------------------------------------------------------------------------
C regreso a menu principal
        IF(COUT.EQ.'y') CLOSE(30)
        IF(CSIMUL.EQ.'y') CLOSE(44)
        GOTO 10
C------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
        END
C
C******************************************************************************
C******************************************************************************
C                               SUBROUTINE PLOTINDEX(NI,NS,RCVEL1,LYMIN,LCOLOR)
C******************************************************************************
C
C Dibuja la region alrededor del indice NI, en el scan NS, tomando como 
C velocidad radial RCVEL1 (en realidad RCVEL1 = 1+z )
C
        SUBROUTINE PLOTINDEX(NI,NS,RCVEL1,LYMIN,LCOLOR)
        IMPLICIT NONE
        INTEGER NI      !numero de indice a dibujar (NI=0 --> todo el espectro)
        INTEGER NS                                              !numero de scan
        REAL RCVEL1                                           !velocidad radial
        LOGICAL LYMIN                            !si .TRUE. YMIN=0 en los plots
        LOGICAL LCOLOR                                !si .TRUE. usamos colores
C
        INCLUDE 'redlib.inc'
        INTEGER TRUEBEG,TRUELEN
C
        INTEGER ITI
        INTEGER J,L,NINDEXT,LL1,LL2
        INTEGER NLINE,NLINEST
        INTEGER NPB,NPBAND
        INTEGER NCONTI,NABSOR
        REAL WV(NWVMAX),FWV(NWVMAX/4),WVLIN(NLINMAX)
        REAL WVBAND(NPBANDMAX),RESBAND(NPBANDMAX)
        REAL WVMIN,WVMAX    !limites reales por si las bandas no estan en orden
        REAL A(NCMAX,NSMAX),ERR(NCMAX,NSMAX)
        REAL WLMIN
        REAL XP(NCMAX),YP(NCMAX),EYP(NCMAX)
        REAL XMIN,XMAX,YMIN,YMAX,DX,DY,XMINL,XMAXL,XMINLR,XMAXLR
        REAL RL1,RL2,DR
        CHARACTER*1 CFHEAD,CFHEAD2,CRVTABLE
        CHARACTER*8 CLABEL,CLABELLIN(NLINMAX)
        CHARACTER*15 CSCANBINNING(NSMAX)
        CHARACTER*20 OBJECTHEAD(NSMAX)
        CHARACTER*50 CDUMMY
        CHARACTER*75 INFILE
        LOGICAL LPLOT,LPERRI,LPLINES,LPBANDS
ccc
        INTEGER NSUBX
        REAL XLEN0,YLEN0
        REAL XINT,XINT2
        REAL VPX1,VPX2,VPY1,VPY2
        REAL PGRND
        integer idum,ixint,ixint2
ccc
C
        COMMON/BLKINDEX0/A,ERR
        COMMON/BLKINDEX2/WLMIN
        COMMON/BLKINDEX3/NSCAN,NCHAN,STWV,DISP
        COMMON/BLKINDEX6/INFILE
        COMMON/BLKINDEX8A/LPLOT,LPERRI           !dibujar o no, con/sin errores
        COMMON/BLKINDEX8B/LPLINES,LPBANDS           !dibujar o no lineas/bandas
        COMMON/BLKINDEX11A/OBJECTHEAD,CSCANBINNING
        COMMON/BLKINDEX11B/CFHEAD,CFHEAD2,CRVTABLE
C------------------------------------------------------------------------------
        DO J=1,NCHAN                                        !espectro a dibujar
          XP(J)=REAL(J)
          YP(J)=A(J,NS)
          EYP(J)=ERR(J,NS)
        END DO
C
        IF(NI.EQ.0)THEN
          XMIN=1
          XMAX=REAL(NCHAN)
          DX=XMAX-XMIN
          XMIN=XMIN-DX/50.
          XMAX=XMAX+DX/50.
          CALL SELINDEX(0,WV,FWV,NINDEXT,CLABEL)     !numero de indices totales
        ELSE
          CALL SELINDEX(NI,WV,FWV,ITI,CLABEL)      !buscamos limites del indice
          IF((ITI.GT.-100).AND.(ITI.LT.-1))THEN
            NCONTI=-ITI                         !numero de regiones de continuo
            WVMIN=WV(1)
            WVMAX=WV(2)
            DO J=2,NCONTI
              IF(WVMIN.GT.WV(2*J-1)) WVMIN=WV(2*J-1)
              IF(WVMAX.LT.WV(2*J)) WVMAX=WV(2*J)
            END DO
            RL1=(WVMIN*RCVEL1-WLMIN)/DISP+1.             !limite azul (channel)
            RL2=(WVMAX*RCVEL1-WLMIN)/DISP                !limite rojo (channel)
          ELSEIF((ITI.EQ.1).OR.(ITI.EQ.2))THEN
            WVMIN=WV(1)
            WVMAX=WV(2)
            DO J=2,3
              IF(WVMIN.GT.WV(2*J-1)) WVMIN=WV(2*J-1)
              IF(WVMAX.LT.WV(2*J)) WVMAX=WV(2*J)
            END DO
            RL1=(WVMIN*RCVEL1-WLMIN)/DISP+1.             !limite azul (channel)
            RL2=(WVMAX*RCVEL1-WLMIN)/DISP                !limite rojo (channel)
          ELSEIF((ITI.EQ.3).OR.(ITI.EQ.4).OR.(ITI.EQ.5))THEN
            RL1=(WV(1)*RCVEL1-WLMIN)/DISP+1.             !limite azul (channel)
            RL2=(WV(4)*RCVEL1-WLMIN)/DISP                !limite rojo (channel)
          ELSEIF((ITI.GE.101).AND.(ITI.LE.9999))THEN
            NCONTI=(ITI/100)                    !numero de regiones de continuo
            NABSOR=ITI-NCONTI*100           !numero de regiones con absorciones
            WVMIN=WV(1)
            WVMAX=WV(2)
            DO J=2,NCONTI+NABSOR
              IF(WVMIN.GT.WV(2*J-1)) WVMIN=WV(2*J-1)
              IF(WVMAX.LT.WV(2*J)) WVMAX=WV(2*J)
            END DO
            RL1=(WVMIN*RCVEL1-WLMIN)/DISP+1.             !limite azul (channel)
            RL2=(WVMAX*RCVEL1-WLMIN)/DISP                !limite rojo (channel)
          ELSE
            STOP 'FATAL ERROR: in subroutine PLOTINDEX'
          END IF
          DR=RL2-RL1
          XMIN=RL1-DR/5.
          XMAX=RL2+DR/5.
        END IF
        YMIN=1.E20
        YMAX=-1.E20
C si hay que dibujar barras de error en los datos exigimos un rango en Y un
C poco mayor (teniendo en cuenta los errores). En caso contrario ajustamos 
C solo los datos (sin los errores)
        IF(LPERRI)THEN
          DO J=1,NCHAN          !buscamos YMAX,YMIN en el rango en X que usamos
            IF((XP(J).GE.XMIN).AND.(XP(J).LE.XMAX))THEN
              IF(YP(J)-EYP(J).LT.YMIN) YMIN=YP(J)-EYP(J)
              IF(YP(J)+EYP(J).GT.YMAX) YMAX=YP(J)+EYP(J)
            END IF
          END DO
        ELSE
          DO J=1,NCHAN          !buscamos YMAX,YMIN en el rango en X que usamos
            IF((XP(J).GE.XMIN).AND.(XP(J).LE.XMAX))THEN
              IF(YP(J).LT.YMIN) YMIN=YP(J)
              IF(YP(J).GT.YMAX) YMAX=YP(J)
            END IF
          END DO
        END IF
        DY=YMAX-YMIN
        YMIN=YMIN-DY/20.
        YMAX=YMAX+DY/20.
        IF(LYMIN) YMIN=0.
        XMINL=(XMIN-1.)*DISP+STWV                            !limites en l.d.o.
        XMAXL=(XMAX-1.)*DISP+STWV
        XMINLR=XMINL/RCVEL1                !limites en l.d.o. corregido de Vrad
        XMAXLR=XMAXL/RCVEL1
C
        CALL PGPAGE
        CALL PGIDEN_RED
C
        CALL PGVPORT(0.05,0.93,0.83,0.89)
        CALL PGSWIN(XMINL,XMAXL,0.,1.)   !escala en l.d.o. sin corregir de Vrad
        CALL PGSLW(3)
        CALL PGBOX('ICMTS',0.,0,'BC',0.,0)
        CALL PGPTEXT(XMAXL+(XMAXL-XMINL)/100.,1.30,0.,0.,
     +   '\\gl\\dobs')
        CALL PGPTEXT(XMAXL+(XMAXL-XMINL)/100.,0.50,0.,0.,
     +   '\\gl\\d(z=0)')
c..............................................................................
c dibujamos correctamente los ticks de la escala corregida de Vrad
        CALL PGLEN(0,'0',XLEN0,YLEN0)  !determinamos dimensiones del caracter 0
        CALL PGQVP(0,VPX1,VPX2,VPY1,VPY2) !dimension del Viewport
        XINT= MAX(0.05, MIN(7.0*XLEN0/(VPX2-VPX1), 0.20))
     +        *(XMAXLR-XMINLR)
        XINT = PGRND(XINT,NSUBX)
        XINT2 = XINT/NSUBX
c afinamos a la centesima de Angstrom
        idum=nint(xminlr*100.)
        ixint=nint(xint*100.)
        if(mod(idum,ixint).ne.0.) idum=(idum/ixint)*ixint
        do while(idum.lt.nint(xmaxlr*100.))
          idum=idum+ixint
          if(idum.lt.nint(xmaxlr*100.))then
            call pgnumb(idum,-2,0,cdummy,l)
            call pgMOVE(real(idum)/100.*rcvel1,1.)
            call pgDRAW(real(idum)/100.*rcvel1,.7)
            call pgptxt(real(idum)/100.*rcvel1,.4,0.,.5,cdummy(1:l))
          end if
        end do
        idum=nint(xminlr*100.)
        ixint2=nint(xint2*100.)
        if(mod(idum,ixint2).ne.0.) idum=(idum/ixint2)*ixint2
        do while(idum.lt.nint(xmaxlr*100.))
          idum=idum+ixint2
          if(idum.lt.nint(xmaxlr*100.))then
            call pgMOVE(real(idum)/100.*rcvel1,1.)
            call pgDRAW(real(idum)/100.*rcvel1,.85)
          end if
        end do
c..............................................................................
        CALL PGVPORT(0.05,0.93,0.10,0.83)
        CALL PGSWIN(XMIN,XMAX,YMIN,YMAX)                     !escala en canales
        CALL PGBOX('IBNTS',0.,0,'BCNST',0.,0)
        CALL PGBIN(NCHAN,XP,YP,.TRUE.)                      !dibujamos espectro
        CALL PGSLW(1)
        IF(LPERRI)THEN                   !hay que dibujar errores en cada pixel
          DO J=1,NCHAN
            IF((REAL(J).GE.XMIN).AND.(REAL(J).LE.XMAX))THEN
              CALL PGMOVE(REAL(J),YP(J)-EYP(J))                 !barra de error
              CALL PGDRAW(REAL(J),YP(J)+EYP(J))
              CALL PGMOVE(REAL(J)-.1,YP(J)-EYP(J))                !tip inferior
              CALL PGDRAW(REAL(J)+.1,YP(J)-EYP(J))
              CALL PGMOVE(REAL(J)-.1,YP(J)+EYP(J))                !tip superior
              CALL PGDRAW(REAL(J)+.1,YP(J)+EYP(J))
            END IF
          END DO
        END IF
        CALL PGSLW(3)
        CALL PGLABEL('channel','No. counts/(response curve)',CHAR(32))
C
        CALL PGMTEXT('T',4.7,0.,0.,'file: '//INFILE)
        IF(NI.EQ.0)THEN
          CALL PGMTEXT('T',4.7,.5,.5,'(all indices)')
        ELSE
          CALL PGMTEXT('T',4.7,.5,.5,'index: '//CLABEL)
        END IF
        WRITE(CDUMMY,*) NS
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        IF((CFHEAD.EQ.'y').OR.(CRVTABLE.EQ.'y'))THEN
          LL1=TRUEBEG(OBJECTHEAD(NS))
          LL2=TRUELEN(OBJECTHEAD(NS))
          CALL PGMTEXT('T',4.7,1.,1.,'scan #'//CDUMMY(1:L)//
     +     ': '//OBJECTHEAD(NS)(LL1:LL2))
        ELSEIF(CFHEAD2.NE.'n')THEN
          LL1=TRUEBEG(CSCANBINNING(NS))
          LL2=TRUELEN(CSCANBINNING(NS))
          CALL PGMTEXT('T',4.7,1.,1.,'scan #'//CDUMMY(1:L)//
     +     ': '//CSCANBINNING(NS)(LL1:LL2))
        ELSE
          CALL PGMTEXT('T',4.7,1.,1.,'scan #'//CDUMMY(1:L))
        END IF
        CALL PGSLW(1)
C
        CALL PGSWIN(XMINL,XMAXL,YMIN,YMAX)  !escala l.d.o. sin corregir de Vrad
C
C NOTA: como en la ultima escala del plot tenemos l.d.o. que no ha sido 
C corregida de Vrad, a la hora de dibujar las bandas necesitamos corregir los 
C limites multiplicando por (1+z) [=RCVEL1]
C Es importante no realizar los calculos en una escala de l.d.o. corregida de
C Vrad (donde dicha escala se aproxima a una escala lineal)
C------------------------------------------------------------------------------
        IF(NI.EQ.0) L=1
10      IF(NI.EQ.0) CALL SELINDEX(L,WV,FWV,ITI,CLABEL)
c..............................................................................
        IF((ITI.GE.-100).AND.(ITI.LT.-1))THEN
          NCONTI=-ITI                           !numero de regiones de continuo
          IF(LCOLOR) CALL PGSCI(5)             
          DO IDUM=1,NCONTI
            CALL PGMOVE(WV(2*IDUM-1)*RCVEL1,YMIN)
            CALL PGDRAW(WV(2*IDUM-1)*RCVEL1,YMAX)
            CALL PGMOVE(WV(2*IDUM)*RCVEL1,YMIN)
            CALL PGDRAW(WV(2*IDUM)*RCVEL1,YMAX)
            CALL PGRECT(WV(2*IDUM-1)*RCVEL1,WV(2*IDUM)*RCVEL1,
     +       YMIN,YMIN+(YMAX-YMIN)/160.)
            CALL PGRECT(WV(2*IDUM-1)*RCVEL1,WV(2*IDUM)*RCVEL1,
     +       YMAX,YMAX-(YMAX-YMIN)/160.)
          END DO
c..............................................................................
        ELSEIF((ITI.GE.1).AND.(ITI.LE.5))THEN
          IF(LCOLOR) CALL PGSCI(4)             !banda azul de todos los indices
          CALL PGMOVE(WV(1)*RCVEL1,YMIN)
          CALL PGDRAW(WV(1)*RCVEL1,YMAX)
          CALL PGMOVE(WV(2)*RCVEL1,YMIN)
          CALL PGDRAW(WV(2)*RCVEL1,YMAX)
          CALL PGRECT(WV(1)*RCVEL1,WV(2)*RCVEL1,
     +     YMIN,YMIN+(YMAX-YMIN)/160.)
          CALL PGRECT(WV(1)*RCVEL1,WV(2)*RCVEL1,
     +     YMAX,YMAX-(YMAX-YMIN)/160.)
          IF(LCOLOR)THEN
            IF(ITI.EQ.3)THEN                              !banda roja del D4000
              CALL PGSCI(2)
            ELSEIF(ITI.EQ.4)THEN                          !banda roja del B4000
              CALL PGSCI(2)
            ELSE                            !banda central de los demas indices
              CALL PGSCI(3)
            END IF
          END IF
          CALL PGMOVE(WV(3)*RCVEL1,YMIN)
          CALL PGDRAW(WV(3)*RCVEL1,YMAX)
          CALL PGMOVE(WV(4)*RCVEL1,YMIN)
          CALL PGDRAW(WV(4)*RCVEL1,YMAX)
          CALL PGRECT(WV(3)*RCVEL1,WV(4)*RCVEL1,
     +     YMIN,YMIN+(YMAX-YMIN)/160.)
          CALL PGRECT(WV(3)*RCVEL1,WV(4)*RCVEL1,
     +     YMAX,YMAX-(YMAX-YMIN)/160.)
          IF((ITI.EQ.1).OR.(ITI.EQ.2))THEN
            IF(LCOLOR) CALL PGSCI(2)                                !banda roja
            CALL PGMOVE(WV(5)*RCVEL1,YMIN)
            CALL PGDRAW(WV(5)*RCVEL1,YMAX)
            CALL PGMOVE(WV(6)*RCVEL1,YMIN)
            CALL PGDRAW(WV(6)*RCVEL1,YMAX)
            CALL PGRECT(WV(5)*RCVEL1,WV(6)*RCVEL1,YMIN,
     +       YMIN+(YMAX-YMIN)/160.)
            CALL PGRECT(WV(5)*RCVEL1,WV(6)*RCVEL1,YMAX,
     +       YMAX-(YMAX-YMIN)/160.)
          END IF
c..............................................................................
        ELSEIF((ITI.GE.101).AND.(ITI.LE.9999))THEN
          NCONTI=(ITI/100)                      !numero de regiones de continuo
          NABSOR=ITI-NCONTI*100             !numero de regiones con absorciones
          IF(LCOLOR) CALL PGSCI(5)             
          DO IDUM=1,NCONTI
            CALL PGMOVE(WV(2*IDUM-1)*RCVEL1,YMIN)
            CALL PGDRAW(WV(2*IDUM-1)*RCVEL1,YMAX)
            CALL PGMOVE(WV(2*IDUM)*RCVEL1,YMIN)
            CALL PGDRAW(WV(2*IDUM)*RCVEL1,YMAX)
            CALL PGRECT(WV(2*IDUM-1)*RCVEL1,WV(2*IDUM)*RCVEL1,
     +       YMIN,YMIN+(YMAX-YMIN)/160.)
            CALL PGRECT(WV(2*IDUM-1)*RCVEL1,WV(2*IDUM)*RCVEL1,
     +       YMAX,YMAX-(YMAX-YMIN)/160.)
          END DO
          IF(LCOLOR) CALL PGSCI(3)             
          DO IDUM=1,NABSOR
            CALL PGMOVE(WV(2*NCONTI+2*IDUM-1)*RCVEL1,YMIN)
            CALL PGDRAW(WV(2*NCONTI+2*IDUM-1)*RCVEL1,
     +       YMAX-2.*(YMAX-YMIN)/160.)
            CALL PGMOVE(WV(2*NCONTI+2*IDUM)*RCVEL1,YMIN)
            CALL PGDRAW(WV(2*NCONTI+2*IDUM)*RCVEL1,
     +       YMAX-2.*(YMAX-YMIN)/160.)
            CALL PGRECT(WV(2*NCONTI+2*IDUM-1)*RCVEL1,
     +       WV(2*NCONTI+2*IDUM)*RCVEL1,YMIN,YMIN+(YMAX-YMIN)/160.)
            CALL PGRECT(WV(2*NCONTI+2*IDUM-1)*RCVEL1,
     +       WV(2*NCONTI+2*IDUM)*RCVEL1,
     +       YMAX-2.*(YMAX-YMIN)/160.,
     +       YMAX-3.*(YMAX-YMIN)/160.)
          END DO
c..............................................................................
        ELSE
          WRITE(*,100) 'ITI='
          WRITE(*,*) ITI
          WRITE(*,101) 'FATAL ERROR: in subroutine PLOTINDEX.'
          WRITE(*,101) 'Invalid ITI value.'
          STOP
        END IF
C------------------------------------------------------------------------------
C si lo hemos indicado, dibujamos lineas tipicas en la region dibujada
        IF(LPLINES)THEN
          CALL SELLINES(0,NLINEST,WVLIN,CLABELLIN)             !serie de Balmer
          IF(LCOLOR) CALL PGSCI(5)
          DO NLINE=1,NLINEST
            WVLIN(NLINE)=WVLIN(NLINE)*RCVEL1
            IF((WVLIN(NLINE).GE.XMINL).AND.(WVLIN(NLINE).LE.XMAXL))THEN
              CALL PGSLS(2)
              CALL PGMOVE(WVLIN(NLINE),YMAX)
              CALL PGDRAW(WVLIN(NLINE),YMIN)
              CALL PGSLS(1)
              CALL PGSCH(0.8)
              CALL PGPTEXT(WVLIN(NLINE),YMAX+DY/100.,0.,.5,
     +         CLABELLIN(NLINE))
              CALL PGSCH(1.0)
            END IF
          END DO
          CALL SELLINES(1,NLINEST,WVLIN,CLABELLIN)      !lineas tipicas emision
          IF(LCOLOR) CALL PGSCI(6)
          DO NLINE=1,NLINEST
            WVLIN(NLINE)=WVLIN(NLINE)*RCVEL1
            IF((WVLIN(NLINE).GE.XMINL).AND.(WVLIN(NLINE).LE.XMAXL))THEN
              CALL PGSLS(2)
              CALL PGMOVE(WVLIN(NLINE),YMAX)
              CALL PGDRAW(WVLIN(NLINE),YMIN)
              CALL PGSLS(1)
              CALL PGSCH(0.8)
              CALL PGPTEXT(WVLIN(NLINE),YMAX+DY/100.,0.,.5,
     +         CLABELLIN(NLINE))
              CALL PGSCH(1.0)
            END IF
          END DO
          CALL SELLINES(2,NLINEST,WVLIN,CLABELLIN)        !lineas tipicas cielo
          IF(LCOLOR) CALL PGSCI(7)
          DO NLINE=1,NLINEST
C NOTA: las lineas de cielo no hay que corregirlas de velocidad
ccc            WVLIN(NLINE)=WVLIN(NLINE)*RCVEL1
            IF((WVLIN(NLINE).GE.XMINL).AND.(WVLIN(NLINE).LE.XMAXL))THEN
              CALL PGSLS(2)
              CALL PGMOVE(WVLIN(NLINE),YMAX)
              CALL PGDRAW(WVLIN(NLINE),YMIN)
              CALL PGSLS(1)
              CALL PGSCH(0.8)
              CALL PGPTEXT(WVLIN(NLINE),YMIN+DY/100.,0.,.5,
     +         CLABELLIN(NLINE))
              CALL PGSCH(1.0)
            END IF
          END DO
        END IF
C------------------------------------------------------------------------------
C Si lo hemos indicado, dibujamos bandas fotometricas en la region dibujada
        IF(LPBANDS)THEN
          CALL SELBANDS('U',NPBAND,WVBAND,RESBAND)                     !banda U
          IF(LCOLOR) CALL PGSCI(4)
          DO NPB=1,NPBAND
            WVBAND(NPB)=WVBAND(NPB)*RCVEL1
            RESBAND(NPB)=YMIN+RESBAND(NPB)*(YMAX-YMIN)*0.30
          END DO
          CALL PGLINE(NPBAND,WVBAND,RESBAND)
          DO NPB=1,NPBAND
            IF((WVBAND(NPB).GE.XMINL).AND.(WVBAND(NPB).LE.XMAXL))THEN
              CALL PGPTEXT(WVBAND(NPB),RESBAND(NPB),0.,0.5,'U')
            END IF
          END DO
          CALL SELBANDS('B',NPBAND,WVBAND,RESBAND)                     !banda B
          IF(LCOLOR) CALL PGSCI(3)
          DO NPB=1,NPBAND
            WVBAND(NPB)=WVBAND(NPB)*RCVEL1
            RESBAND(NPB)=YMIN+RESBAND(NPB)*(YMAX-YMIN)*0.30
          END DO
          CALL PGLINE(NPBAND,WVBAND,RESBAND)
          DO NPB=1,NPBAND
            IF((WVBAND(NPB).GE.XMINL).AND.(WVBAND(NPB).LE.XMAXL))THEN
              CALL PGPTEXT(WVBAND(NPB),RESBAND(NPB),0.,0.5,'B')
            END IF
          END DO
          CALL SELBANDS('V',NPBAND,WVBAND,RESBAND)                     !banda V
          IF(LCOLOR) CALL PGSCI(2)
          DO NPB=1,NPBAND
            WVBAND(NPB)=WVBAND(NPB)*RCVEL1
            RESBAND(NPB)=YMIN+RESBAND(NPB)*(YMAX-YMIN)*0.30
          END DO
          CALL PGLINE(NPBAND,WVBAND,RESBAND)
          DO NPB=1,NPBAND
            IF((WVBAND(NPB).GE.XMINL).AND.(WVBAND(NPB).LE.XMAXL))THEN
              CALL PGPTEXT(WVBAND(NPB),RESBAND(NPB),0.,0.5,'V')
            END IF
          END DO
        END IF
C------------------------------------------------------------------------------
        IF(LCOLOR) CALL PGSCI(1)
        IF(NI.NE.0)THEN
          RETURN
        ELSE
          L=L+1
          IF(L.GT.NINDEXT) RETURN
          GOTO 10
        END IF
100     FORMAT(A,$)
101     FORMAT(A)
        END
C
C******************************************************************************
C Determina la fraccion de banda fotometrica CBAND que cae dentro del espectro
        SUBROUTINE BANDFRACTION(CBAND,RCVEL1,FRACTION)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        CHARACTER*1 CBAND
        REAL RCVEL1
        REAL FRACTION
C
        INTEGER J,J1,J2
        INTEGER NPBAND,NPB
        INTEGER IMIN(NCMAX),IMAX(NCMAX)
        REAL X1,X2,WL0
        REAL WVBAND(NPBANDMAX),RESBAND(NPBANDMAX)
        REAL RESCHAN(NCMAX)
        REAL WLMIN,WLMAX
        REAL SUMTOT
        LOGICAL LOK
C
        COMMON/BLKINDEX3/NSCAN,NCHAN,STWV,DISP
        COMMON/BLKINDEX10C/RESCHAN
C------------------------------------------------------------------------------
C obtenemos los datos de la banda CBAND
        CALL SELBANDS(CBAND,NPBAND,WVBAND,RESBAND)
        DO NPB=1,NPBAND
          WVBAND(NPB)=WVBAND(NPB)*RCVEL1
        END DO
        WLMIN=STWV-DISP/2.
        WLMAX=WLMIN+REAL(NCHAN)*DISP
C caso trivial: banda fuera del espectro
        IF((WVBAND(1).GT.WLMAX).OR.(WVBAND(NPBAND).LT.WLMIN))THEN
          FRACTION=0.0
          RETURN
        END IF
C caso no trivial: banda en espectro (parcial o totalmente)
C calculamos la respueta de la banda en cada canal
        DO J=1,NCHAN
          IMIN(J)=0
          IMAX(J)=0
        END DO
        DO NPB=1,NPBAND-1
          X1=(WVBAND(NPB)-STWV)/DISP+1.
          X2=(WVBAND(NPB+1)-STWV)/DISP+1.
          J1=NINT(X1)
          J2=NINT(X2)
          LOK=.TRUE.
          IF(J1.LT.1) J1=MIN(1,J2)
          IF(J2.GT.NCHAN) J2=MAX(NCHAN,J1)
          IF(J2.LT.1) LOK=.FALSE.
          IF(J1.GT.NCHAN) LOK=.FALSE.
          IF(LOK)THEN
            DO J=J1,J2
              IMIN(J)=NPB
              IMAX(J)=NPB+1
            END DO
          END IF
        END DO
        DO J=1,NCHAN
          WL0=STWV+REAL(J-1)*DISP
          IF((IMIN(J).NE.0).AND.(IMAX(J).NE.0))THEN
            RESCHAN(J)=RESBAND(IMIN(J))+
     +       (RESBAND(IMAX(J))-RESBAND(IMIN(J)))*
     +       (WL0-WVBAND(IMIN(J)))/(WVBAND(IMAX(J))-WVBAND(IMIN(J)))
          ELSE
            RESCHAN(J)=0.
          END IF
        END DO
C calculamos la fraccion de banda
        SUMTOT=0.                        !integral de toda la banda fotometrica
        DO NPB=1,NPBAND-1
          SUMTOT=SUMTOT+0.5*(RESBAND(NPB)+RESBAND(NPB+1))*
     +     (WVBAND(NPB+1)-WVBAND(NPB))
        END DO
        FRACTION=0.     !integral de la banda en la region de espectro cubierta
        DO J=1,NCHAN
          FRACTION=FRACTION+RESCHAN(J)*DISP
        END DO
        FRACTION=FRACTION/SUMTOT                             !fraccion de banda
C
        END
C
C******************************************************************************
C Calcula el indice de una ley de potencias de exponente ALPHA
        SUBROUTINE INDEXPLAW(NINDEX,ALPHA,FINDEX)
        IMPLICIT NONE
C
        INCLUDE 'redlib.inc'
C
        INTEGER NINDEX
        REAL ALPHA
        REAL FINDEX
C
        INTEGER ITI
        REAL WV(NWVMAX),FWV(NWVMAX/4)
        REAL LA,LC,LR,FA,FC,FR,FCONT
        CHARACTER*8 CLABEL
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        CALL SELINDEX(NINDEX,WV,FWV,ITI,CLABEL)
        IF((ITI.GT.-100).AND.(ITI.LT.-1))THEN
          WRITE(*,*)
          WRITE(*,101) 'ERROR: Power law has not been implemented for '
          WRITE(*,101) '       this type of index yet. Sorry.'
          WRITE(*,*)
          RETURN
        ELSEIF((ITI.GE.101).AND.(ITI.LE.9999))THEN
          WRITE(*,*)
          WRITE(*,101) 'ERROR: Power law has not been implemented for '
          WRITE(*,101) '       this type of index yet. Sorry.'
          WRITE(*,*)
          RETURN
        ELSEIF(ITI.EQ.3)THEN                                          !D4000
          FINDEX=(WV(4)**(ALPHA+1.)-WV(3)**(ALPHA+1.))
          FINDEX=FINDEX/(WV(2)**(ALPHA+1.)-WV(1)**(ALPHA+1.))
        ELSEIF(ITI.EQ.4)THEN                                          !B4000
          FINDEX=(WV(4)**(ALPHA-1.)-WV(3)**(ALPHA-1.))
          FINDEX=FINDEX/(WV(2)**(ALPHA-1.)-WV(1)**(ALPHA-1.))
        ELSE                                    !indices atomicos y moleculares
          LA=(WV(1)+WV(2))/2.
          LC=(WV(3)+WV(4))/2.
          LR=(WV(5)+WV(6))/2.
          IF(ALPHA.EQ.1.0)THEN
            FA=ALOG(WV(2))-ALOG(WV(1))
            FC=ALOG(WV(4))-ALOG(WV(3))
            FR=ALOG(WV(6))-ALOG(WV(5))
          ELSE
            FA=WV(2)**(ALPHA-1.)-WV(1)**(ALPHA-1.)
            FC=WV(4)**(ALPHA-1.)-WV(3)**(ALPHA-1.)
            FR=WV(6)**(ALPHA-1.)-WV(5)**(ALPHA-1.)
          END IF
          FA=FA/(WV(2)-WV(1))
          FC=FC/(WV(4)-WV(3))
          FR=FR/(WV(6)-WV(5))
          FCONT=0.5*(FA+FR)+0.5*(FR-FA)*(2.*LC-LA-LR)/(LR-LA)
          IF(ITI.EQ.1)THEN                                    !indice molecular
            FINDEX=-2.5*ALOG10(FC/FCONT)
          ELSE                                                  !indice atomico
            FINDEX=(WV(4)-WV(3))*(1.-FC/FCONT)
          END IF
        END IF
C
101     FORMAT(A)
        END
