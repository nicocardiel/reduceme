C------------------------------------------------------------------------------
C Version 13-October-2007                                       file: snratio.f
C @ ncl & fjg
C------------------------------------------------------------------------------
Comment
C
C Program: snratio
C Classification: error handling
C Description: Determines the S/N ratio as a function of binning in the spatial
C direction.
C
Comment
C
C Calcula la relacion senhal/ruido en una imagen como funcion de un binning
C cualquiera. Tambien permite estudiar como cambian los indices y errores en
C funcion del binning realizado.
C Se ha incluido la posibilidad de realizar plots bidimensionales (en
C funcion del scan inicial y final sumados), que representan la relacion
C senhal/ruido, los indices medidos y los errores.
C
        PROGRAM SNRATIO
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER TRUELEN
        INTEGER READI
        INTEGER READILIM
        REAL READF
C
C parametros 
C (NOTA: sin se cambian revisar subrutinas)
        INTEGER NCOLMAX                 !numero maximo de columnas en la salida
        PARAMETER (NCOLMAX=30) 
        INTEGER NFMAX                        !numero maximo de curvas respuesta
        PARAMETER (NFMAX=101)
        REAL C                                             !velocidad de la luz
        PARAMETER (C=2.9979246E+5)
C
        INTEGER I,J,K,L,M,N
        INTEGER I1,I2,I1MAX,I2MAX
        INTEGER IMAX
        INTEGER NCONTI,NABSOR,IWL
        INTEGER NS,NS0,ITI
        INTEGER NINDEXT,NINDEX,NIND1,NIND2
        INTEGER NBIN1,NBIN,NCOL
        INTEGER NS1,NS2,NFS1,NFS2,NSPLOT
        INTEGER IL1,IL2
        INTEGER ID0,MIDEIND
        INTEGER NCREST,NCHAN2
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        INTEGER NBIN_OPT,NS1_OPT(NSMAX),NS2_OPT(NSMAX)
        REAL A(NCMAX,NSMAX),ERR(NCMAX,NSMAX),FLUX(NCMAX,NFMAX)
        REAL PLOTNS(NSMAX,NSMAX)
        REAL SNRAT(NSMAX,NCOLMAX),SN,SNSIM(NSMAX)
        REAL FINDEXSIM(NSMAX),EINDEXSIM(NSMAX)
        REAL SNRATMIN,SNRATMAX
        REAL FINDEXMAX,FINDEXMIN,XDUM
        REAL XT(NCMAX),YT(NCMAX),XPLOT(NSMAX),YPLOT(NSMAX)
        REAL XMIN,XMAX,YMIN,YMAX,DX,DY,YMINDUM,YMAXDUM,DYDUM
        REAL YMAXF
        REAL RVEL,RCVEL,RCVEL1
        REAL WLMIN,WLMAX,WLMINZ,WLMAXZ,WLMIN0,WLMAX0
        REAL WV(NWVMAX),FWV(NWVMAX/4)
        REAL FINDEX,EINDEX,EJJGG,ESIMU
        REAL STWV2,DISP2
        REAL TR(6),FG,BG,XC,YC
        REAL SNMIN
        CHARACTER*1 CPLOT,CMODE,COUT,CFF,CNEW,CFCAL,CCONT,CMERR
        CHARACTER*1 CWHATP,CCHANGE,CH,CCONFIRM
        CHARACTER*1 COUT_BIN
        CHARACTER*8 CLABEL
        CHARACTER*50 CDUMMY
        CHARACTER*75 INFILE,ERRFILE,OUTFILE,FLUXFILE,CBINFILE
        LOGICAL LANY,LEXIT
        LOGICAL LINDOK(NINDMAX)
        LOGICAL LCOLOR(MAX_ID_RED)
C
        COMMON/BLKINDEX0/A,ERR
        COMMON/BLKINDEX1/FLUX
        COMMON/BLKINDEX2/WLMIN
        COMMON/BLKINDEX3/NSCAN,NCHAN,STWV,DISP
        COMMON/BLKINDEX4/NCREST
        COMMON/BLKSNF0/NBIN_OPT,NS1_OPT,NS2_OPT
        COMMON/BLKSNF1/NCOL,IL1,IL2
        COMMON/BLKSNF2/COUT
        COMMON/BLKSNF3/SNRAT
C------------------------------------------------------------------------------
        THISPROGRAM='snratio'
        CALL WELCOME('13-October-2007')
C
        CALL CHEQUEA_FILEINDEX
C------------------------------------------------------------------------------
        CMODE='a'
C
        TR(1)=0.
        TR(2)=1.
        TR(3)=0.
        TR(4)=0.
        TR(5)=0.
        TR(6)=1.
C
        WRITE(*,100)'Plots (y/n) '
        CPLOT(1:1)=READC('y','yn')
C abrimos la salida grafica
        IF(CPLOT.EQ.'y')THEN
          CALL PIDEGTER(NTERM,IDN,LCOLOR)
        END IF
C leemos el fichero de datos
5       WRITE(*,100)'Input file name......'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        CALL GUESSEF(INFILE,ERRFILE)         !predecimos nombre fichero errores
        IF(TRUELEN(OBJECT).GT.0)THEN
          INFILE=INFILE(1:TRUELEN(INFILE))//' ['//
     +     OBJECT(1:TRUELEN(OBJECT))//']'
        END IF
        DO I=1,NSCAN
          READ(20) (A(J,I),J=1,NCHAN)
        END DO
        CLOSE(20)
C leemos el fichero de errores
        WRITE(*,100)'Input error file name '
        ERRFILE=INFILEX(22,ERRFILE,NSCAN,NCHAN,STWV,DISP,21,.TRUE.) !.....match
        DO I=1,NSCAN
          READ(22) (ERR(J,I),J=1,NCHAN)
        END DO
        CLOSE(22)
        WRITE(*,*)
C
        NFS1=1
        NFS2=NSCAN
C leemos curva respuesta
        WRITE(*,100)'Flux calibration (y/n) '
        CFCAL(1:1)=READC('n','yn')
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
          END IF
          IF(NCREST.GT.NFMAX)THEN
            WRITE(*,101)'FATAL ERROR: number of response curves '//
     +       'too large. Redim NFMAX.'
            CLOSE(23)
            STOP
          END IF
          READ(23) (FLUX(J,1),J=1,NCHAN)
          CLOSE(23)
          DO I=1,NSCAN
            DO J=1,NCHAN
              A(J,I)=A(J,I)/FLUX(J,1)
            END DO
            DO J=1,NCHAN
             ERR(J,I)=ERR(J,I)/FLUX(J,1)
            END DO
          END DO
        ELSE
          FLUXFILE='[NONE]'
          NCREST=1              !es igual a 1 para saltar proteccion en MIDEIND
        END IF
C------------------------------------------------------------------------------
C inicializamos a cero
        DO J=1,NCHAN
          YT(J)=0.
        END DO
C dibujamos corte en la direccion espacial y buscamos el maximo
        WRITE(*,100)'Wait...'
        IMAX=0
        DO I=1,NSCAN
          XT(I)=REAL(I)
          DO J=1,NCHAN
            YT(I)=YT(I)+A(J,I)/REAL(NCHAN)
          END DO
        END DO
        IF(CPLOT.EQ.'y')THEN
          XMIN=XT(1)
          XMAX=XT(NSCAN)
          DX=XMAX-XMIN
          XMIN=XMIN-DX/50.
          XMAX=XMAX+DX/50.
        END IF
7       IF(CPLOT.EQ.'y')THEN
          YMIN=1.0E33
          YMAX=-YMIN
          DO J=1,NCHAN
            IF((XT(J).GE.XMIN).AND.(XT(J).LE.XMAX))THEN
              IF(YT(J).LT.YMIN) YMIN=YT(J)
              IF(YT(J).GT.YMAX) YMAX=YT(J)
            END IF
          END DO
          DY=YMAX-YMIN
          YMIN=YMIN-DY/50.
          YMAX=YMAX+DY/50.
        END IF
        YMAXF=YMIN
        DO I=1,NSCAN
          IF((XT(I).GE.XMIN).AND.(XT(I).LE.XMAX))THEN
            IF(YT(I).GT.YMAXF)THEN
              YMAXF=YT(I)
              IMAX=I
            END IF
          END IF
        END DO
        IF(CPLOT.EQ.'y')THEN
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,0)
            CALL PGIDEN_RED
            IF(LCOLOR(ITERM)) CALL PGSCI(2)
            CALL PGRECT(XT(IMAX)-.5,XT(IMAX)+.5,0.0,YMAXF)
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
            CALL PGBIN(NSCAN,XT,YT,.TRUE.)
            CALL PGLABEL('scan','averaged no. of counts','file: '//
     +       INFILE)
          END DO
        END IF
        WRITE(*,101)'OK!'
        WRITE(*,110)'>>> Maximum is located at scan: ',IMAX
C------------------------------------------------------------------------------
C Pedimos posicion definitiva del maximo
ccc10      IF(CPLOT.EQ.'y')THEN
        IF(CPLOT.EQ.'y')THEN
          WRITE(*,100)'Center scan (0=replot with zoom) '
          WRITE(CDUMMY,*) IMAX
          NS0=READILIM(CDUMMY,0,NSCAN)
          IF(NS0.EQ.0)THEN
            WRITE(CDUMMY,*)XMIN
            WRITE(*,100)'Xmin '
            XMIN=READF(CDUMMY)
            WRITE(CDUMMY,*)XMAX
            WRITE(*,100)'Xmax '
            XMAX=READF(CDUMMY)
            GOTO 7
          END IF
        ELSE
          WRITE(*,100)'Center scan '
          WRITE(CDUMMY,*) IMAX
          NS0=READILIM(CDUMMY,1,NSCAN)
        END IF
C------------------------------------------------------------------------------
C velocidad radial
        WRITE(*,*)
        WRITE(*,100)'Radial velocity (km/sec) '
        RVEL=READF('0.0')
        RCVEL=RVEL/C
        RCVEL1=1.+RCVEL                                 !(1+z)
        RCVEL1=RCVEL1/SQRT(1.-(RCVEL*RCVEL)/(C*C))      !correccion relativista
C------------------------------------------------------------------------------
C l.d.o. del extremo izquierdo del primer pixel (STWV se refiere al centro del
C primer pixel) y del extremo derecho del ultimo pixel a Vrad=0 y a Vrad
C correspondiente a la velocidad radial introducida
        WLMIN=STWV-DISP/2.0
        WLMAX=STWV+REAL(NCHAN-1)*DISP+DISP/2.0
        WRITE(*,*)
        WRITE(CDUMMY,*)WLMIN
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,100)'>>> Limits: from '//CDUMMY(1:L)
        WRITE(CDUMMY,*)WLMAX
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,101)' to '//CDUMMY(1:L)//
     +   ' angstroms (at 1+z = 1.0000)'
C
        WLMINZ=WLMIN/RCVEL1                                !l.d.o. minima mayor
        WLMAXZ=WLMAX/RCVEL1                                !l.d.o. maxima menor
        IF(WLMINZ.GE.WLMAXZ)THEN
          WRITE(*,101)'FATAL ERROR: WLMINZ.eq.WLMAXZ'
          IF(CPLOT.EQ.'y') CALL PGEND
          STOP
        END IF
        WRITE(CDUMMY,*)WLMINZ
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,100)'>>> Limits: from '//CDUMMY(1:L)
        WRITE(CDUMMY,*)WLMAXZ
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,100)' to '//CDUMMY(1:L)//
     +   ' angstroms (at 1+z = '
        WRITE(CDUMMY,'(F6.4)') RCVEL1
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,101)CDUMMY(1:L)//')'
C------------------------------------------------------------------------------
C determinamos que indices pueden ser medidos (ver subrutina SELINDEX)
        LANY=.FALSE.
        CALL SELINDEX(0,WV,FWV,NINDEXT,CLABEL)    !numero de indices disponible
        DO K=1,NINDEXT
          CALL SELINDEX(K,WV,FWV,ITI,CLABEL)
          IF((ITI.EQ.1).OR.(ITI.EQ.2))THEN
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
            WRITE(*,101)'FATAL ERROR: invalid ITI value.'
          END IF
          IF(LINDOK(K))THEN
            LANY=.TRUE.
          END IF
        END DO
        IF(.NOT.LANY)THEN
          WRITE(*,*)
          WRITE(*,101)'WARNING: No. of indices which can be '//
     +     'measured = 0'
          IF(CPLOT.EQ.'y') CALL PGEND
          STOP
        END IF
C------------------------------------------------------------------------------
C elegimos los indices a medir
30      WRITE(*,*)
        CALL SHINDEX(LINDOK,0)
32      WRITE(*,100)'Index number..................................'
        NINDEX=READI('-1')
        IF(NINDEX.EQ.-1) THEN
          WRITE(*,100)'[n]ew image or [q]uit (n/q) '
          CNEW(1:1)=READC('q','qn')
          IF(CNEW.EQ.'q')THEN
            IF(CPLOT.EQ.'y') CALL PGEND
            STOP
          ELSE
            WRITE(*,*)
            GOTO 5
          END IF
        END IF
        IF((NINDEX.LT.0).OR.(NINDEX.GT.NINDEXT))THEN
          WRITE(*,101)'ERROR: number out of range. Try again.'
          GOTO 32
        END IF
        IF(NINDEX.EQ.0)THEN
          NIND1=1
          NIND2=NINDEXT
        ELSE
          IF(.NOT.LINDOK(NINDEX))THEN
            WRITE(*,101)'ERROR: this index is not available. Try again.'
            GOTO 32
          END IF
          NIND1=NINDEX
          NIND2=NINDEX
        END IF
C------------------------------------------------------------------------------
C Elegimos metodo de medida de la relacion S/N
        IF(CPLOT.EQ.'y')THEN
          WRITE(*,100)'[a]utomatic, [m]anual or [p]lot scan '//
     +     'selection (a/m/p) '
          CMODE(1:1)=READC(CMODE,'amp')
        ELSE
          WRITE(*,100)'[a]utomatic, or [m]anual scan '//
     +     'selection (a/m) '
          CMODE(1:1)=READC(CMODE,'am')
        END IF
35      IF((CMODE.EQ.'m').OR.(CMODE.EQ.'p'))THEN
          WRITE(*,*)
36        WRITE(*,100)'1st and last scan (0,0=EXIT).............'
          WRITE(CDUMMY,'(I10,A1,I10)')NFS1,',',NFS2
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          CALL READ2I(CDUMMY(1:L),NFS1,NFS2)
          IF((NFS1.EQ.0).AND.(NFS2.EQ.0)) GOTO 30
          IF((NFS1.LT.1).OR.(NFS2.GT.NSCAN).OR.(NFS1.GT.NFS2))THEN
            WRITE(*,101)'ERROR: invalid entry. Try again.'
            GOTO 36
          END IF
          IF((NS0.LT.NFS1).OR.(NS0.GT.NFS2))THEN
            WRITE(*,101)'ERROR: center scan out of range. Try again.'
            GOTO 36
          END IF
          COUT='n'
          GOTO 60
        END IF
C------------------------------------------------------------------------------
C Fijamos region (scans) a utilizar
40      WRITE(*,100)'1st and last scan to be employed.........'
        WRITE(CDUMMY,'(I10,A1,I10)')NFS1,',',NFS2
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        CALL READ2I(CDUMMY(1:L),NFS1,NFS2)
        IF((NFS1.LT.1).OR.(NFS2.GT.NSCAN).OR.(NFS1.GT.NFS2))THEN
          WRITE(*,101)'ERROR: invalid entry. Try again.'
          GOTO 40
        END IF
        IF((NS0.LT.NFS1).OR.(NS0.GT.NFS2))THEN
          WRITE(*,101)'ERROR: center scan out of range. Try again.'
          GOTO 40
        END IF
C------------------------------------------------------------------------------
        WRITE(*,100)'Create output file name (y/n)...............'
        COUT(1:1)=READC('n','yn')
        IF(COUT.EQ.'y')THEN
          WRITE(*,100)'Output file name.......................'//
     +     '..........'
          OUTFILE=OUTFILEX(30,'@',0,0,0.,0.,3,.FALSE.)
          IF((NINDEX.EQ.0).AND.(NINDEXT.GT.1))THEN      !insercion de Form Feed
            WRITE(*,100)'Insert <FF> between indices (y/n)...........'
            CFF(1:1)=READC('y','yn')
          ELSE
            CFF='n'
          END IF
        END IF
C------------------------------------------------------------------------------
C bucle en indices
60      DO K=NIND1,NIND2                                      !bucle en indices
          IF(LINDOK(K))THEN                    !solo si el indice puede medirse
C..............................................................................
            CALL SELINDEX(K,WV,FWV,ITI,CLABEL)
C..............................................................................
C cabecera con indice y nombre del fichero
            WRITE(*,*)
            WRITE(*,101)'-----------------------------------------'//
     +       '--------------------------------------'
            WRITE(*,'(A11,I2,A10,$)')'>>> Index #',K,': '//CLABEL
            IF(ITI.EQ.1)THEN
              WRITE(*,100)' (molecular) '
            ELSEIF(ITI.EQ.2)THEN
              WRITE(*,100)' (atomic)    '
            ELSEIF(ITI.EQ.3)THEN
              WRITE(*,100)'             '
            ELSEIF(ITI.EQ.4)THEN
              WRITE(*,100)'             '
            ELSEIF(ITI.EQ.5)THEN
              WRITE(*,100)'             '
            ELSEIF((ITI.GE.101).AND.(ITI.LE.9999))THEN
              WRITE(*,100)' (generic)   '
            END IF
            L=TRUELEN(INFILE)+TRUELEN(ERRFILE)+2
            IF(L+42.LT.79)THEN
              DO M=1,37-L
                WRITE(*,100)' '
              END DO
            ELSE
              WRITE(*,*)
            END IF
            WRITE(*,101)'Files: '//INFILE(1:TRUELEN(INFILE))//
     +       ','//ERRFILE(1:TRUELEN(ERRFILE))
            WRITE(*,100)'>>> Radial velocity (km/sec): '
            WRITE(*,'(F7.0)')RVEL
            WRITE(*,101)'>>> Response curve file name: '//
     +       FLUXFILE(1:TRUELEN(FLUXFILE))
            WRITE(*,100)'>>> Center scan: '
            WRITE(*,*) NS0
C
            IF(COUT.EQ.'y')THEN
              WRITE(30,101)'#-------------------------------------'//
     +         '-----------------------------------------'
              WRITE(30,'(A11,I2,A10,$)')'#>>> Index #',K,': '//CLABEL
              IF(ITI.EQ.1)THEN
                WRITE(30,100)' (molecular) '
              ELSEIF(ITI.EQ.2)THEN
                WRITE(30,100)' (atomic)    '
              ELSEIF(ITI.EQ.3)THEN
                WRITE(30,100)'             '
              ELSEIF(ITI.EQ.4)THEN
                WRITE(30,100)'             '
              ELSEIF(ITI.EQ.5)THEN
                WRITE(30,100)'             '
              ELSEIF((ITI.GE.101).AND.(ITI.LE.9999))THEN
                WRITE(30,100)' (generic)   '
              END IF
              L=TRUELEN(INFILE)+TRUELEN(ERRFILE)+2
              IF(L+42.LT.79)THEN
                DO M=1,37-L
                  WRITE(30,100)' '
                END DO
              ELSE
                WRITE(30,*)
              END IF
              WRITE(30,101)'#Files: '//INFILE(1:TRUELEN(INFILE))//
     +         ','//ERRFILE(1:TRUELEN(ERRFILE))
              WRITE(30,100)'#>>> Radial velocity (km/sec): '
              WRITE(30,'(F7.0)')RVEL
              WRITE(30,101)'#>>> Response curve file name: '//
     +         FLUXFILE(1:TRUELEN(FLUXFILE))
              WRITE(30,100)'#>>> Center scan: '
              WRITE(30,*) NS0
            END IF
C..............................................................................
C binning simetrico
            DO I=1,NSCAN                                !inicializamos variable
              SNSIM(I)=0.
              FINDEXSIM(I)=0.
              EINDEXSIM(I)=0.
            END DO
            LEXIT=.FALSE.
            NBIN=-1
            DO WHILE(.NOT.LEXIT)
              NBIN=NBIN+1
              NS1=NS0-NBIN
              IF(NS1.LT.NFS1) NS1=NFS1
              NS2=NS0+NBIN
              IF(NS2.GT.NFS2) NS2=NFS2
              ID0=MIDEIND(NS1,NS2,ITI,WV,FWV,'y',RCVEL1,1,FINDEX,
     +         EINDEX,EJJGG,ESIMU,SN,-1,.FALSE.)
              LEXIT=((NS1.EQ.NFS1).AND.(NS2.EQ.NFS2))
              NS1=NS0-NBIN
              NS2=NS0+NBIN
              IF(NS1.GE.NFS1)THEN
                SNSIM(NS1) = SN
                FINDEXSIM(NS1) = FINDEX
                EINDEXSIM(NS1) = EINDEX
              END IF
              IF(NS2.LE.NFS2)THEN
                SNSIM(NS2) = SN
                FINDEXSIM(NS2) = FINDEX
                EINDEXSIM(NS2) = EINDEX
              END IF
            END DO
C
            IF(CPLOT.EQ.'y')THEN
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                CALL PGPAGE
                CALL PGVSTAND
                CALL PGWINDOW(REAL(NFS1),REAL(NFS2),YMIN,YMAX)
                CALL PGBOX('BNTS',0.,0,' ',0.,0)
                CALL PGBOX('CMTS',0.,0,' ',0.,0)
                CALL PGIDEN_RED
                CALL PGLABEL('scan',' ',' ')
                WRITE(CDUMMY,*)K
                CALL RMBLANK(CDUMMY,CDUMMY,L)
                CALL PGMTEXT('T',3.,1.,1.,
     +           'Index #'//CDUMMY(1:L)//': '//CLABEL)
                CALL PGMTEXT('T',3.,0.,0.,'file: '//INFILE)
                WRITE(CDUMMY,*)NS0
                CALL RMBLANK(CDUMMY,CDUMMY,L)
                CALL PGMTEXT('B',3.,0.,0.,'Symmetric binning around '//
     +           CDUMMY(1:L))
                IF(LCOLOR(ITERM)) CALL PGSCI(4)
                CALL PGBIN(NSCAN,XT,YT,.TRUE.)
              END DO
C
              FINDEXMIN=FINDEXSIM(NFS1)-EINDEXSIM(NFS1)
              FINDEXMAX=FINDEXSIM(NFS1)+EINDEXSIM(NFS1)
              DO I=NFS1+1,NFS2
                XDUM=FINDEXSIM(I)-EINDEXSIM(I)
                IF(XDUM.LT.FINDEXMIN) FINDEXMIN=XDUM
                XDUM=FINDEXSIM(I)+EINDEXSIM(I)
                IF(XDUM.GT.FINDEXMAX) FINDEXMAX=XDUM
              END DO
              DY=FINDEXMAX-FINDEXMIN
              FINDEXMIN=FINDEXMIN-DY/20.
              FINDEXMAX=FINDEXMAX+DY/20.
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                CALL PGWINDOW(REAL(NFS1),REAL(NFS2),FINDEXMIN,FINDEXMAX)
                IF(LCOLOR(ITERM)) CALL PGSCI(1)
                CALL PGBOX(' ',0.,0,'BNST',0.,0)
                CALL PGMTEXT('L',3.,.5,.5,'Index')
                CALL PGLINE(NSCAN,XT,FINDEXSIM)
                IF(LCOLOR(ITERM)) CALL PGSCI(7)
                DO I=NFS1,NFS2
                  CALL PGMOVE(XT(I),FINDEXSIM(I)-EINDEXSIM(I))
                  CALL PGDRAW(XT(I),FINDEXSIM(I)+EINDEXSIM(I))
                  CALL PGMOVE(XT(I)-.1,FINDEXSIM(I)-EINDEXSIM(I))
                  CALL PGDRAW(XT(I)+.1,FINDEXSIM(I)-EINDEXSIM(I))
                  CALL PGMOVE(XT(I)-.1,FINDEXSIM(I)+EINDEXSIM(I))
                  CALL PGDRAW(XT(I)+.1,FINDEXSIM(I)+EINDEXSIM(I))
                END DO
              END DO
C
              SNRATMAX=SNSIM(NFS1)
              SNRATMIN=SNSIM(NFS1)
              DO I=NFS1+1,NFS2
                IF(SNSIM(I).LT.SNRATMIN) SNRATMIN=SNSIM(I)
                IF(SNSIM(I).GT.SNRATMAX) SNRATMAX=SNSIM(I)
              END DO
              DY=SNRATMAX-SNRATMIN
              SNRATMIN=SNRATMIN-DY/50.
              SNRATMAX=SNRATMAX+DY/50.
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                CALL PGWINDOW(REAL(NFS1),REAL(NFS2),SNRATMIN,SNRATMAX)
                IF(LCOLOR(ITERM)) CALL PGSCI(2)
                CALL PGBOX(' ',0.,0,'CMST',0.,0)
                CALL PGMTEXT('R',3.,.5,.5,'<Signal/Noise>')
                CALL PGLINE(NSCAN,XT,SNSIM)
                IF(LCOLOR(ITERM)) CALL PGSCI(3)
                CALL PGSLS(2)
              END DO
            END IF
            IL1=NS0
            DO WHILE(((IL1-1).GE.1).AND.(SNSIM(IL1-1).GT.SNSIM(IL1)))
              IL1=IL1-1
            END DO
            IF(CPLOT.EQ.'y')THEN
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                CALL PGMOVE(REAL(IL1),SNRATMIN)
                CALL PGDRAW(REAL(IL1),SNRATMAX)
              END DO
            END IF
            WRITE(CDUMMY,*)IL1
            CALL RMBLANK(CDUMMY,CDUMMY,L)
            IF(CPLOT.EQ.'y')THEN
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                CALL PGPTEXT(REAL(IL1),SNRATMAX+(SNRATMAX-SNRATMIN)/15.,
     +           0.,.5,CDUMMY(1:L))
              END DO
            END IF
            WRITE(*,100)'>>> Optimum symmetric scan region: ['//
     +       CDUMMY(1:L)//','
            IF(COUT.EQ.'y')
     +       WRITE(30,100)'#>>> Optimum symmetric scan region: ['//
     +       CDUMMY(1:L)//','
            IL2=NS0
            DO WHILE(((IL2+1).LE.NSCAN).AND.
     +       (SNSIM(IL2+1).GT.SNSIM(IL2)))
              IL2=IL2+1
            END DO
            IF(CPLOT.EQ.'y')THEN
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                CALL PGMOVE(REAL(IL2),SNRATMIN)
                CALL PGDRAW(REAL(IL2),SNRATMAX)
              END DO
            END IF
            WRITE(CDUMMY,*)IL2
            CALL RMBLANK(CDUMMY,CDUMMY,L)
            IF(CPLOT.EQ.'y')THEN
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                CALL PGPTEXT(REAL(IL2),SNRATMAX+(SNRATMAX-SNRATMIN)/15.,
     +           0.,.5,CDUMMY(1:L))
              END DO
            END IF
            WRITE(*,100)CDUMMY(1:L)//'],  S/N = '
            IF(COUT.EQ.'y') WRITE(30,100)CDUMMY(1:L)//'],  S/N = '
            SNRATMAX=SNSIM(IL1)
            IF(IL2.GT.IL1)THEN
              DO I=IL1+1,IL2
                IF(SNSIM(I).GT.SNRATMAX) SNRATMAX=SNSIM(I)
              END DO
            END IF
            WRITE(CDUMMY,*)SNRATMAX
            CALL RMBLANK(CDUMMY,CDUMMY,L)
            WRITE(*,101)CDUMMY(1:L)
            IF(COUT.EQ.'y') WRITE(30,101)CDUMMY(1:L)
            IF(CPLOT.EQ.'y')THEN
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                CALL PGSLS(1)
                IF(LCOLOR(ITERM)) CALL PGSCI(1)
              END DO
            END IF
            WRITE(*,101)'-----------------------------------------'//
     +       '--------------------------------------'
            IF(COUT.EQ.'y')THEN
              WRITE(30,101)'#-------------------------------------'//
     +         '-----------------------------------------'
            END IF
C..............................................................................
            IF(CMODE.EQ.'a')THEN
C mostramos resultado detallado del binning simetrico
              WRITE(*,101)'scan     S/N     index    error'
              DO I=NFS1,NS0
                WRITE(*,'(I4,3X,$)') I
                IF((I.GE.IL1).AND.(I.LE.IL2))THEN
                  WRITE(*,'(1A,F8.2,2(1X,F8.3))') '*',SNSIM(I),
     +             FINDEXSIM(I),EINDEXSIM(I)
                ELSE
                  WRITE(*,'(1X,F8.2,2(1X,F8.3))') SNSIM(I),FINDEXSIM(I),
     +             EINDEXSIM(I)
                END IF
              END DO
              WRITE(*,101)'--------------------------------------'//
     +         '-----------------------------------------'
              DO I=NS0,NFS2
                WRITE(*,'(I4,3X,$)') I
                IF((I.GE.IL1).AND.(I.LE.IL2))THEN
                  WRITE(*,'(1A,F8.2,2(1X,F8.3))') '*',SNSIM(I),
     +             FINDEXSIM(I),EINDEXSIM(I)
                ELSE
                  WRITE(*,'(1X,F8.2,2(1X,F8.3))') SNSIM(I),FINDEXSIM(I),
     +             EINDEXSIM(I)
                END IF
              END DO
              WRITE(*,101)'--------------------------------------'//
     +         '-----------------------------------------'
C lo mismo en fichero
              IF(COUT.EQ.'y')THEN
                WRITE(*,100)'Save last results into output file (y/n) '
                CCONFIRM(1:1)=READC('n','yn')
                IF(CCONFIRM.EQ.'y')THEN
                  WRITE(30,101)'#scan    S/N     index    error'
                  DO I=NFS1,NS0
                    WRITE(30,'(I4,3X,$)') I
                    IF((I.GE.IL1).AND.(I.LE.IL2))THEN
                      WRITE(30,'(1A,F8.2,2(1X,F8.3))') '*',SNSIM(I),
     +                 FINDEXSIM(I),EINDEXSIM(I)
                    ELSE
                      WRITE(30,'(1X,F8.2,2(1X,F8.3))') SNSIM(I),
     +                 FINDEXSIM(I),EINDEXSIM(I)
                    END IF
                  END DO
                  WRITE(30,101)'#--------------------------------'//
     +             '----------------------------------------------'
                  DO I=NS0,NFS2
                    WRITE(30,'(I4,3X,$)') I
                    IF((I.GE.IL1).AND.(I.LE.IL2))THEN
                      WRITE(30,'(1A,F8.2,2(1X,F8.3))') '*',SNSIM(I),
     +                 FINDEXSIM(I),EINDEXSIM(I)
                    ELSE
                      WRITE(30,'(1X,F8.2,2(1X,F8.3))') SNSIM(I),
     +                 FINDEXSIM(I),EINDEXSIM(I)
                    END IF
                  END DO
                  WRITE(30,101)'#--------------------------------'//
     +             '----------------------------------------------'
                END IF
              END IF
C binning no simetrico
              WRITE(*,100)'Continue with asymmetric binning (y/n) '
              CCONT(1:1)=READC('y','yn')
              IF(CCONT.EQ.'y')THEN
C Binning incial y numero de columnas
                WRITE(*,100)'Initial bin width '
                NBIN1=READI('1')
                WRITE(*,100)'No. of columns '
                NCOL=READILIM('8',1,NCOLMAX)
                WRITE(*,100)'Minimum S/N to mark binning'
                SNMIN=READF('@')
                WRITE(*,101)'--------------------------------------'//
     +           '-----------------------------------------'
                WRITE(*,101)'Asymmetric binning...'
                WRITE(*,100)'Minimum S/N to mark binning: '
                WRITE(*,*) SNMIN
                WRITE(*,101)'--------------------------------------'//
     +           '-----------------------------------------'
                IF(COUT.EQ.'y')THEN
                  WRITE(30,101)'#----------------------------------'//
     +             '--------------------------------------------'
                  WRITE(30,101)'#Asymmetric binning...'
                  WRITE(30,100)'#Minimum S/N to mark binning: '
                  WRITE(30,*) SNMIN
                  WRITE(30,101)'#----------------------------------'//
     +             '--------------------------------------------'
                END IF
                WRITE(*,100)' scan'
                DO L=1,NCOL
                  IF(L.EQ.NCOL)THEN
                    WRITE(*,'(6X,I3)')NBIN1+L-1
                  ELSE
                    WRITE(*,'(6X,I3,$)')NBIN1+L-1
                  END IF
                END DO
                IF(COUT.EQ.'y')THEN
                  WRITE(30,100)'#scan'
                  DO L=1,NCOL
                    IF(L.EQ.NCOL)THEN
                      WRITE(30,'(6X,I3)')NBIN1+L-1
                    ELSE
                      WRITE(30,'(6X,I3,$)')NBIN1+L-1
                    END IF
                  END DO
                END IF
C
                NBIN_OPT=0 !de momento no hemos determinado ningun binning
C 1st HALF
                DO NBIN=NBIN1,NBIN1+NCOL
                  DO NS=NFS1,NS0
                    NS1=NS-NBIN+1
                    NS2=NS
                    IF(NS1.GE.NFS1)THEN
                      ID0=MIDEIND(NS1,NS2,ITI,WV,FWV,'y',RCVEL1,1,
     +                 FINDEX,EINDEX,EJJGG,ESIMU,SN,-1,.TRUE.)
                      SNRAT(NS,NBIN-NBIN1+1)=SN
                    ELSE
                      SNRAT(NS,NBIN-NBIN1+1)=0.
                    END IF
                  END DO
                END DO
                CALL SHOWSN(NFS1,NS0,NBIN1,1,SNMIN)
                WRITE(*,101)'--------------------------------------'//
     +           '-----------------------------------------'
                IF(COUT.EQ.'y')THEN
                  WRITE(30,101)'#----------------------------------'//
     +             '--------------------------------------------'
                END IF
C 2nd HALF
                DO NBIN=NBIN1,NBIN1+NCOL
                  DO NS=NS0,NFS2
                    NS1=NS
                    NS2=NS+NBIN-1
                    IF(NS2.LE.NFS2)THEN
                      ID0=MIDEIND(NS1,NS2,ITI,WV,FWV,'y',RCVEL1,1,
     +                 FINDEX,EINDEX,EJJGG,ESIMU,SN,-1,.TRUE.)
                      SNRAT(NS,NBIN-NBIN1+1)=SN
                    ELSE
                      SNRAT(NS,NBIN-NBIN1+1)=0.
                    END IF
                  END DO
                END DO
                CALL SHOWSN(NS0,NFS2,NBIN1,2,SNMIN)
                WRITE(*,101)'-----------------------------------'//
     +           '--------------------------------------------'
                IF(COUT.EQ.'y')THEN
                  WRITE(30,101)'#----------------------------------'//
     +             '--------------------------------------------'
                  IF(CFF.EQ.'y') WRITE(30,101) CHAR(12)
                END IF
              END IF
C podemos salvar un fichero con solo los scans del binning determinado
              WRITE(*,100)'Create output file name with computed '//
     +         'binning regions '
              COUT_BIN(1:1)=READC('n','yn')
              IF(COUT_BIN.EQ.'y')THEN
                CALL ORDENA2I(NBIN_OPT,NS1_OPT,NS2_OPT) !ordenamos
                WRITE(*,101)
                WRITE(*,101) 'WARNING: the central scan will appear '//
     +           'twice in the next file. In the case the'
                WRITE(*,101) 'central binning is 1, you '//
     +           'can easily remove the duplicated entry (and the'
                WRITE(*,101) 'comments: lines starting by "#") using '//
     +           'the command line:'
                WRITE(*,101) '% cat filename | grep -v "^#" | uniq '//
     +           '> new_filename'
                WRITE(*,100)'Output file name'
                CBINFILE=OUTFILEX(34,'@',0,0,0.,0.,3,.FALSE.)
                WRITE(34,101)'#-------------------------------------'//
     +           '-----------------------------------------'
                WRITE(34,'(A11,I2,A10,$)')'#>>> Index #',K,': '//CLABEL
                IF(ITI.EQ.1)THEN
                  WRITE(34,100)' (molecular) '
                ELSEIF(ITI.EQ.2)THEN
                  WRITE(34,100)' (atomic)    '
                ELSEIF(ITI.EQ.3)THEN
                  WRITE(34,100)'             '
                ELSEIF(ITI.EQ.4)THEN
                  WRITE(34,100)'             '
                ELSEIF(ITI.EQ.5)THEN
                  WRITE(34,100)'             '
                ELSEIF((ITI.GE.101).AND.(ITI.LE.9999))THEN
                  WRITE(34,100)' (generic)   '
                END IF
                L=TRUELEN(INFILE)+TRUELEN(ERRFILE)+2
                IF(L+42.LT.79)THEN
                  DO M=1,37-L
                    WRITE(34,100)' '
                  END DO
                ELSE
                  WRITE(34,*)
                END IF
                WRITE(34,101)'#Files: '//INFILE(1:TRUELEN(INFILE))//
     +           ','//ERRFILE(1:TRUELEN(ERRFILE))
                WRITE(34,100)'#>>> Radial velocity (km/sec): '
                WRITE(34,'(F7.0)')RVEL
                WRITE(34,101)'#>>> Response curve file name: '//
     +           FLUXFILE(1:TRUELEN(FLUXFILE))
                WRITE(34,100)'#>>> Center scan: '
                WRITE(34,*) NS0
                WRITE(34,101)'#----------------------------------'//
     +           '--------------------------------------------'
                WRITE(34,101)'#Asymmetric binning...'
                WRITE(34,100)'#Minimum S/N to mark binning: '
                WRITE(34,*) SNMIN
                WRITE(34,101)'#----------------------------------'//
     +           '--------------------------------------------'
                IF(NBIN_OPT.GT.0)THEN
                  DO N=1,NBIN_OPT
                    WRITE(34,'(I4,1X,I4)') NS1_OPT(N),NS2_OPT(N)
                  END DO
                END IF
                CLOSE(34)
              END IF
C..............................................................................
            ELSEIF(CMODE.EQ.'m')THEN
              WRITE(*,100)'Measure index & error (y/n) '
              CMERR(1:1)=READC('y','yn')
              IF(CMERR.EQ.'y')THEN
                ID0=MIDEIND(NS1,NS2,ITI,WV,FWV,'y',RCVEL1,1,FINDEX,
     +           EINDEX,EJJGG,ESIMU,SN,-1,.FALSE.)
              ELSE
                ID0=MIDEIND(NS1,NS2,ITI,WV,FWV,'y',RCVEL1,1,FINDEX,
     +           EINDEX,EJJGG,ESIMU,SN,-1,.TRUE.)
              END IF
              WRITE(*,101)'--------------------------------------'//
     +         '-----------------------------------------'
              WRITE(CDUMMY,'(I10,A1,I10)')NS1,',',NS2
              CALL RMBLANK(CDUMMY,CDUMMY,L)
              WRITE(*,100)'S/N with binning ['//CDUMMY(1:L)//'] = '
              WRITE(*,*) SN
              IF(CMERR.EQ.'y')THEN
                WRITE(*,100)'Index: '
                WRITE(*,*) FINDEX
                WRITE(*,100)'Error: '
                WRITE(*,*) EINDEX
              END IF
              WRITE(*,101)'--------------------------------------'//
     +         '-----------------------------------------'
C..............................................................................
C Realizamos representaciones bidimensionales
            ELSEIF(CMODE.EQ.'p')THEN
              DO I1=1,NSCAN        !inicializamos a cero la mitad no utilizable
                DO I2=1,NSCAN
                  PLOTNS(I2,I1)=-999.
                END DO
              END DO
              CWHATP=' '
              DO WHILE(CWHATP.NE.'x')
                WRITE(*,100)'Plot [s]/n ratio, [i]ndex, [e]rrors or '//
     +           'e[x]it (i/e/s/x) '
                CWHATP(1:1)=READC('x','iesx')
                IF(CWHATP.NE.'x')THEN
                  WRITE(*,100)'Please, wait...     '
                  IF(CWHATP.EQ.'s')THEN
                    CDUMMY='S/N ratio: '//CLABEL
                  ELSEIF(CWHATP.EQ.'i')THEN
                    CDUMMY='Index: '//CLABEL
                  ELSEIF(CWHATP.EQ.'e')THEN
                    CDUMMY='Errors: '//CLABEL
                  END IF
                  DO I1=NFS1,NFS2
                    WRITE(*,'(A,I4.4,$)')'\b\b\b\b',I1
                    DO I2=I1,NFS2
                      ID0=MIDEIND(I1,I2,ITI,WV,FWV,'y',RCVEL1,1,
     +                 FINDEX,EINDEX,EJJGG,ESIMU,SN,-1,.FALSE.)
                      IF(CWHATP.EQ.'s')THEN
                        PLOTNS(I2,I1)=SN
                      ELSEIF(CWHATP.EQ.'i')THEN
                        PLOTNS(I2,I1)=FINDEX
                      ELSEIF(CWHATP.EQ.'e')THEN
                        PLOTNS(I2,I1)=EINDEX
                      END IF
                    END DO
                  END DO
                  WRITE(*,101)'\b\b\b\b  ...OK!'
                  BG=PLOTNS(NFS1,NFS1)
                  FG=BG
                  I1MAX=NFS1
                  I2MAX=NFS2
                  DO I1=NFS1,NFS2
                    DO I2=I1,NFS2
                      IF(PLOTNS(I2,I1).LT.BG) BG=PLOTNS(I2,I1)
                      IF(PLOTNS(I2,I1).GT.FG)THEN
                        I1MAX=I1
                        I2MAX=I2
                        FG=PLOTNS(I2,I1)
                      END IF
                    END DO
                  END DO
                  WRITE(*,*)
                  WRITE(*,100)'BG: '
                  WRITE(*,*)BG
                  WRITE(*,100)'FG: '
                  WRITE(*,*)FG
                  WRITE(*,100)'FG (First scan) :'
                  WRITE(*,*)I1MAX
                  WRITE(*,100)'FG (Last  scan) :'
                  WRITE(*,*)I2MAX
                  ID0=MIDEIND(I1MAX,I2MAX,ITI,WV,FWV,'y',RCVEL1,1,
     +             FINDEX,EINDEX,EJJGG,ESIMU,SN,-1,.FALSE.)
                  WRITE(*,100)'Index...........: '
                  WRITE(*,*)FINDEX
                  WRITE(*,100)'Error...........: '
                  WRITE(*,*)EINDEX
                  WRITE(*,100)'S/N ratio.......: '
                  WRITE(*,*)SN
                  DO ITERM=NTERM,1,-1
                    CALL PGSLCT(IDN(ITERM))
                    CALL PGENV(REAL(NFS1-1),REAL(NFS2+1),REAL(NFS1-1),
     +               REAL(NFS2+1),1,0)
                    CALL PGGRAY(PLOTNS,NSMAX,NSMAX,1,NSCAN,1,NSCAN,
     +               FG,BG,TR)
                    CALL PGLABEL('last scan','first scan',CDUMMY)
                    IF(LCOLOR(ITERM))THEN
                      CALL PGSCI(2)
                      CALL PGMOVE(REAL(I2MAX)-.5,REAL(I1MAX)-.5)
                      CALL PGDRAW(REAL(I2MAX)+.5,REAL(I1MAX)+.5)
                      CALL PGMOVE(REAL(I2MAX)-.5,REAL(I1MAX)+.5)
                      CALL PGDRAW(REAL(I2MAX)+.5,REAL(I1MAX)-.5)
                      CALL PGSCI(1)
                    END IF
                  END DO
                  CCHANGE=' '
                  DO WHILE(CCHANGE.NE.'x')
                    WRITE(*,*)
                    WRITE(*,100)'Use [m]ouse, [s]et BG/FG, '//
     +               'xy [c]uts, or e[x]it (m/s/c/x) '
                    CCHANGE(1:1)=READC('x','mscx')
                    IF(CCHANGE.EQ.'x')THEN
                    ELSEIF(CCHANGE.EQ.'m')THEN
                      CH=' '
                      DO WHILE(CH.NE.'X')
                        CALL PGBAND(0,0,0.,0.,XC,YC,CH)
                        IF(CH.EQ.'q') CH='X'
                        IF(CH.EQ.'Q') CH='X'
                        IF(CH.NE.'X')THEN
                          WRITE(*,*)
                          WRITE(*,100)'First scan....: '
                          WRITE(*,*)NINT(YC)
                          WRITE(*,100)'Last scan.....: '
                          WRITE(*,*)NINT(XC)
                          WRITE(*,100)'Pixel value...: '
                          WRITE(*,*)PLOTNS(NINT(XC),NINT(YC))
                        END IF
                      END DO
                      WRITE(*,*)
                    ELSEIF(CCHANGE.EQ.'s')THEN
                      WRITE(CDUMMY,*) BG
                      WRITE(*,100)'New BG '
                      BG=READF(CDUMMY)
                      WRITE(CDUMMY,*) FG
                      WRITE(*,100)'New FG '
                      FG=READF(CDUMMY)
                      DO ITERM=NTERM,1,-1
                        CALL PGSLCT(IDN(ITERM))
                        CALL PGGRAY(PLOTNS,NSMAX,NSMAX,1,NSCAN,1,NSCAN,
     +                   FG,BG,TR)
                      END DO
                    ELSEIF(CCHANGE.EQ.'c')THEN
                      WRITE(*,100)'Scan '
                      NSPLOT=READILIM('@',NFS1,NFS2)
                      DO I2=NFS1,NFS2
                        XPLOT(I2-NFS1+1)=REAL(I2)
                        YPLOT(I2-NFS1+1)=PLOTNS(I2,NSPLOT)
                      END DO
                      DYDUM=FG-BG
                      YMINDUM=BG-DYDUM/50.
                      YMAXDUM=FG+DYDUM/50.
                      DO ITERM=NTERM,1,-1
                        CALL PGSLCT(IDN(ITERM))
                        CALL PGENV(REAL(NFS1-1),REAL(NFS2+1),
     +                   REAL(NFS1-1),REAL(NFS2+1),1,0)
                        CALL PGGRAY(PLOTNS,NSMAX,NSMAX,1,NSCAN,1,NSCAN,
     +                   FG,BG,TR)
                        CALL PGLABEL('last scan','first scan',' ')
                        IF(LCOLOR(ITERM)) CALL PGSCI(2)
                        CALL PGSLS(2)
                        CALL PGMOVE(REAL(NFS1-1),REAL(NSPLOT))
                        CALL PGDRAW(REAL(NFS2+1),REAL(NSPLOT))
                        IF(LCOLOR(ITERM)) CALL PGSCI(3)
                        CALL PGMOVE(REAL(NSPLOT),REAL(NFS1-1))
                        CALL PGDRAW(REAL(NSPLOT),REAL(NFS2+1))
                        CALL PGSLS(1)
                        IF(LCOLOR(ITERM)) CALL PGSCI(2)
                        CALL PGWINDOW(REAL(NFS1-1),REAL(NFS2+1),
     +                   YMINDUM,YMAXDUM)
                        CALL PGBOX(' ',0.0,0,'CMTS',0.0,0)
                        CALL PGLINE(NFS2-NFS1+1,XPLOT,YPLOT)
                        CALL PGMTEXT('R',2.5,0.5,0.5,CDUMMY)
                        IF(LCOLOR(ITERM)) CALL PGSCI(3)
                      END DO
                      DO I1=NFS1,NFS2
                        XPLOT(I1-NFS1+1)=PLOTNS(NSPLOT,I1)
                        YPLOT(I1-NFS1+1)=REAL(I1)
                      END DO
                      DO ITERM=NTERM,1,-1
                        CALL PGSLCT(IDN(ITERM))
                        CALL PGWINDOW(YMIN,YMAX,REAL(NFS1-1),
     +                   REAL(NFS2+1))
                        CALL PGBOX('CMTS',0.0,0,' ',0.0,0)
                        CALL PGLINE(NFS2-NFS1+1,XPLOT,YPLOT)
                        CALL PGMTEXT('T',2.,0.5,0.5,CDUMMY)
                        IF(LCOLOR(ITERM)) CALL PGSCI(1)
                        CALL PGWINDOW(REAL(NFS1-1),REAL(NFS2+1),
     +                   REAL(NFS1-1),REAL(NFS2+1))
                      END DO
                    END IF
                  END DO
                END IF
              END DO
C..............................................................................
            END IF
          END IF
          IF(NINDEX.EQ.0)THEN
            WRITE(*,100)'Press <RETURN> to continue with next index...'
            READ(*,*)
          END IF
        END DO                                                !bucle en indices
C------------------------------------------------------------------------------
C regreso a menu principal o nueva entrada manual
        IF(CMODE.EQ.'a')THEN
          IF(COUT.EQ.'y') CLOSE(30)
          GOTO 30
        ELSEIF(CMODE.EQ.'p')THEN
          GOTO 30
        ELSE
          GOTO 35
        END IF
C------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
110     FORMAT(A,I6)
        END
C
C******************************************************************************
C
C Muestra en varias columnas la informacion calculada
C NUPDOWN=1 se mueve desde NS2 hacia NS1 (arriba)
C NUPDOWN=2 se mueve desde NS1 hacia NS2 (abajo)
        SUBROUTINE SHOWSN(NS1,NS2,NBIN1,NUPDOWN,SNMIN)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
C
C parametros de la subrutina
        INTEGER NS1,NS2
        INTEGER NBIN1,NUPDOWN
        REAL SNMIN
C
C parametros globales del programa
C (NOTA: sin se cambia, revisar todas las definiciones)
        INTEGER NCOLMAX
        PARAMETER (NCOLMAX=30)          !numero maximo de columnas en la salida
C
C variables globales
        INTEGER NBIN_OPT,NS1_OPT(NSMAX),NS2_OPT(NSMAX)
        INTEGER NCOL,IL1,IL2
        REAL SNRAT(NSMAX,NCOLMAX)
        CHARACTER*1 COUT
        CHARACTER*1 CSNRAT(NSMAX,NCOLMAX)
C
C variables locales
        INTEGER I
        INTEGER NC
        INTEGER NS0,DNS0,NC0
        INTEGER NS0_TEMP
C
        COMMON/BLKSNF0/NBIN_OPT,NS1_OPT,NS2_OPT
        COMMON/BLKSNF1/NCOL,IL1,IL2
        COMMON/BLKSNF2/COUT
        COMMON/BLKSNF3/SNRAT
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
C determinamos el binning optimo que verifica SN.GE.SNMIN
        IF(NUPDOWN.EQ.1)THEN
          NS0=NS2
          DNS0=-1
        ELSE
          NS0=NS1
          DNS0=+1
        END IF
        NC0=1 !empezamos buscando sin binning
C
        DO I=NS1,NS2
          DO NC=1,NCOL
            CSNRAT(I,NC)=' '
          END DO
        END DO
C
10      IF((NS0.LT.NS1).OR.(NS0.GT.NS2).OR.(NC0.GT.NCOL)) GOTO 90
        IF(SNRAT(NS0,NC0).GT.SNMIN)THEN
          IF(NBIN1+(NC0-1).EQ.1)THEN
            CSNRAT(NS0,NC0)='*'
            NBIN_OPT=NBIN_OPT+1
            NS1_OPT(NBIN_OPT)=NS0
            NS2_OPT(NBIN_OPT)=NS0
            NS0=NS0+DNS0
          ELSE
            CSNRAT(NS0,NC0)='*'
            NS0_TEMP=NS0
            IF(NBIN1+(NC0-1).GT.2)THEN
              DO I=2,NBIN1+(NC0-1)-1
                NS0=NS0+DNS0
                CSNRAT(NS0,NC0)='|'
              END DO
            END IF
            NS0=NS0+DNS0
            CSNRAT(NS0,NC0)='<'
            NBIN_OPT=NBIN_OPT+1
            IF(NUPDOWN.EQ.1)THEN 
              NS1_OPT(NBIN_OPT)=NS0
              NS2_OPT(NBIN_OPT)=NS0_TEMP
            ELSE
              NS1_OPT(NBIN_OPT)=NS0_TEMP
              NS2_OPT(NBIN_OPT)=NS0
            END IF
            NS0=NS0+DNS0
          END IF
        ELSE
          NC0=NC0+1
        END IF
        GOTO 10
C------------------------------------------------------------------------------
C mostramos resultado
90      DO I=NS1,NS2
          WRITE(*,'(I5,A1,$)') I,' '
          IF(COUT.EQ.'y') WRITE(30,'(I5,A1,$)') I,' '
          DO NC=1,NCOL
            IF(NC.EQ.NCOL)THEN
              WRITE(*,'(F8.2,A1)') SNRAT(I,NC),CSNRAT(I,NC)
              IF(COUT.EQ.'y') WRITE(30,'(F8.2,A1)') SNRAT(I,NC),
     +         CSNRAT(I,NC)
            ELSE
              WRITE(*,'(F8.2,A1,$)') SNRAT(I,NC),CSNRAT(I,NC)
              IF(COUT.EQ.'y') WRITE(30,'(F8.2,A1,$)') SNRAT(I,NC),
     +         CSNRAT(I,NC)
            END IF
          END DO
        END DO
        END
