C------------------------------------------------------------------------------
C Version 23-July-1998                                       file: interlines.f
C @ ncl & fjg
C------------------------------------------------------------------------------
Comment
C
C Program: interlines
C Classification: miscellany
C Description: Interpolates lines in a spectrum by using data from a similar 
C spectrum.
C
Comment
C
C Interpola lineas en un espectro sustituyendo pixels por la senhal en un
C espectro similar. Tambien permite medir lineas de emision y ajusta
C gaussianas.
C
        PROGRAM INTERLINES
        IMPLICIT NONE
C Include NSMAX,NCMAX,NSCAN,NCHAN
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER TRUELEN
        INTEGER READILIM
        REAL READF
C
        REAL C                                             !velocidad de la luz
        PARAMETER (C=299792.46)
C
        INTEGER I,J,K,L,L1,L2
        INTEGER NS0,NN
        INTEGER NS1,NS2
        INTEGER NSCAN2,NCHAN2            !para la curva de calibracion en flujo
        INTEGER NB,NBLOCAL
        INTEGER NINDEX
        INTEGER NFIT,NTERMS
        INTEGER NDIM,MP,JP,NEVAL                       !variables para DOWNHILL
        INTEGER NTERMS_CONT
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        INTEGER NITER
        REAL X(NCMAX)                                          !numero de canal
        REAL S(NCMAX),SS(NCMAX)                            !espectro del objeto
        REAL ES(NCMAX),ESS(NCMAX)          !error asociado al espectro anterior
        REAL T(NCMAX),TT(NCMAX)                              !espectro template
        REAL ET(NCMAX),ETT(NCMAX)          !error asociado al espectro anterior
        REAL D(NCMAX),ED(NCMAX)     !diferencia entre objeto y template y error
        REAL P(NCMAX)                          !polinomio para igualar continuo
        REAL Q(NCMAX)                                         !ley de potencias
        REAL A(NCMAX,NSMAX),ERRA(NCMAX,NSMAX)      !imagen de trabajo y errores
        REAL FLUX(NCMAX)                         !curva de calibracion en flujo
        REAL TEMP(NCMAX),ERRTEMP(NCMAX)              !espectro template y error
        REAL XMIN,XMAX,YMIN,YMAX,DY,XMINL,XMAXL,YMIN0,YMAX0
        REAL TMEAN,SMEAN
        REAL STWV2,DISP2                 !para la curva de calibracion en flujo
        REAL RVGAL,RVTEMP,RVEFF
        REAL GAMMAGAL,GAMMATEMP,CTE
        REAL XC,YC
        REAL WLMIN,WLMAX,RCVEL1,WV0,WV0MFIT
        REAL XFIT(NCMAX),YFIT(NCMAX),CHISQR,COEFF(20),POL
        REAL POL1,POL2,POL3,POL4
        REAL YRMSTOL,XX0(7),DXX0(7),XXF(7),DXXF(7)     !variables para DOWNHILL
        EXTERNAL FUNKPL,FUNKMPFIT,FUNKPOLT
        REAL FUNKPL,FUNKMPFIT,FUNKPOLT                 !funciones para DOWNHILL
        REAL XXOLD(7),DXXOLD(7)
        REAL K1,K2,ALPHA                  !parametros de ajuste a una power law
        REAL AIRMASS0,TIMEXPOS0
        CHARACTER*1 CERR,CFLUX,CH,CSAVE,CREP,CTFIT
        CHARACTER*1 COPC,COK,CINTER,CINTER0,CSURE,CSAVETXT
        CHARACTER*1 CSAME,CASCII,CMEAS,CMEASLAST
        CHARACTER*8 CLABEL0
        CHARACTER*75 INFILE,ERRINFILE,FLUXFILE,FILEASCII,ERRFILEASCII
        CHARACTER*75 TEMPFILE,ERRTEMPFILE,OUTFILE,TEXTFILE
        CHARACTER*75 CDUMMY
        CHARACTER*80 COMENTARIO
        CHARACTER*255 OBJECT0,FITSFILE0,COMMENT0
        LOGICAL IFSCAN(NSMAX)
        LOGICAL LPINDEX,LPLINES,LPERROR,LPBANDS,LPREGUS,LPRESIDUALS
        LOGICAL LINDOK(NINDMAX)
        LOGICAL LREGIONS          !normaliza los espectros a regiones definidas
        LOGICAL IFCHAN(NCMAX)         !canales a usar para normalizar espectros
        LOGICAL LFIT           !indica si hemos ajustado un polinomio/power law
        LOGICAL LYMINMAX                            !valores de YMIN YMAX fijos
        LOGICAL IFCHAN_INTER(NCMAX),IFCHAN_POLYN(NCMAX)
        LOGICAL IFCHAN_CONT(NCMAX),IFCHAN_LINE(NCMAX)
        LOGICAL IFCHAN_MFIT(NCMAX)
        LOGICAL LBEXIST                                        !el boton existe
        LOGICAL LBUFFER                !=.TRUE. si ya hemos medido alguna linea
        LOGICAL LCOLOR(MAX_ID_RED)
C
        COMMON/BLKSHOW/LINDOK
        COMMON/BLKGEN1/STWV,DISP
        COMMON/BLKGEN2/RCVEL1,WLMIN,WLMAX
        COMMON/BLKGEN3/YMIN,YMAX
        COMMON/BLKFUNK1/IFCHAN
        COMMON/BLKFUNK2/S,T
        COMMON/BLKFUNK3/NCHAN
        COMMON/BLKFUNKMPFIT/IFCHAN_MFIT
        COMMON/BLKMWIDTH0/CSAVETXT
        COMMON/BLKMWIDTH1/LBUFFER
        COMMON/BLKMWIDTH2/IFCHAN_CONT,IFCHAN_LINE
        COMMON/BLKMWIDTH3/NTERMS_CONT
        COMMON/BLKMWIDTH4/NS0
        COMMON/BLKMWIDTH5/CMEASLAST
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
        THISPROGRAM='interlines'
        CALL WELCOME('13-April-1998')
C Inicializamos algunas variables
        LPINDEX=.FALSE.
        LPLINES=.FALSE.
        LPERROR=.FALSE.
        LPBANDS=.FALSE.
        LPREGUS=.FALSE.
        LREGIONS=.FALSE.
        LPRESIDUALS=.TRUE.
        LFIT=.FALSE.
        LYMINMAX=.FALSE.
        LBUFFER=.FALSE.
        CINTER='0'
        CINTER0='1'
        CMEASLAST='0'
        YRMSTOL=1.E-4
C evitamos algunos warnings de compilacion
        YMIN0=0.0
        YMAX0=0.0
C------------------------------------------------------------------------------
C abrimos salida grafica
        CALL RPGBEGIN(NTERM,IDN,LCOLOR)
        CALL BUTTSYB(3)
C------------------------------------------------------------------------------
        WRITE(*,100)'Work with error images (y/n) '
        CERR(1:1)=READC('y','yn')
        IF(CERR.EQ.'y')THEN
          NITER=100
        ELSE
          NITER=0
        END IF
C------------------------------------------------------------------------------
C Leemos fichero(s) de entrada
        WRITE(*,100)'Input file name'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        DO I=1,NSCAN
          READ(20) (A(J,I),J=1,NCHAN)
        END DO
        CLOSE(20)
        OBJECT0=OBJECT
        FITSFILE0=FITSFILE
        COMMENT0=COMMENT
        AIRMASS0=AIRMASS
        TIMEXPOS0=TIMEXPOS
        IF(CERR.EQ.'y')THEN
          WRITE(*,100)'Input error file name '
          CALL GUESSEF(INFILE,ERRINFILE)
          ERRINFILE=INFILEX(21,ERRINFILE,NSCAN,NCHAN,STWV,DISP,
     +     21,.TRUE.)                                                 !match
          DO I=1,NSCAN
            READ(21) (ERRA(J,I),J=1,NCHAN)
          END DO
          CLOSE(21)
        END IF
C
        WRITE(*,100)'Radial velocity (km/sec)'
        RVGAL=READF('@')
C------------------------------------------------------------------------------
C coordenada X para plots
        DO J=1,NCHAN
          X(J)=REAL(J)
        END DO
C------------------------------------------------------------------------------
        IF(NSCAN.GT.1)THEN
          WRITE(CDUMMY,*)NSCAN
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          WRITE(*,101)'Valid range: from 1 to '//CDUMMY(1:L)
          WRITE(*,100)'First scan to start with '
          NS0=READILIM('@',1,NSCAN)
        ELSE
          NS0=1
        END IF
C------------------------------------------------------------------------------
C espectro de estrella template
        WRITE(*,100)'Template spectrum file name'
        TEMPFILE=INFILEX(20,'@',NSCAN2,NCHAN2,STWV2,DISP2,1,.FALSE.)
        IF(NCHAN2.NE.NCHAN)THEN
          CLOSE(20)
          WRITE(*,100)'FATAL ERROR: NCHAN in last image is different!'
          STOP
        END IF
        IF(STWV2.NE.STWV)THEN
          CLOSE(20)
          WRITE(*,100)'FATAL ERROR: STWV in last image is different!'
          STOP
        END IF
        IF(DISP2.NE.DISP)THEN
          CLOSE(20)
          WRITE(*,100)'FATAL ERROR: DISP in last image is different!'
          STOP
        END IF
        IF(NSCAN2.GT.1)THEN
          WRITE(*,101)'WARNING: last image contain more than '//
     +     '1 spectrum!'
          WRITE(*,101)'Enter scan region to be employed: '
          DO I=1,NSCAN2
            IFSCAN(I)=.FALSE.
          END DO
3         WRITE(*,100)'Enter scan region (0,0=EXIT) '
          CALL READ2I('0,0',NS1,NS2)
          IF((NS1.EQ.0).AND.(NS2.EQ.0)) GOTO 4
          IF((NS1.LT.1).OR.(NS2.GT.NSCAN2).OR.(NS1.GT.NS2))THEN
            WRITE(*,101)'ERROR: limits out of range. Try again.'
            GOTO 3
          END IF
          DO I=NS1,NS2
            IFSCAN(I)=.TRUE.
          END DO
          GOTO 3
4         NS1=0
          DO I=1,NSCAN2
            IF(IFSCAN(I)) NS1=NS1+1
          END DO
          IF(NS1.EQ.0)THEN
            WRITE(*,101)'ERROR: total number of selected scans = 0'
            GOTO 3
          END IF
          DO J=1,NCHAN
            TEMP(J)=0.
          END DO
          DO I=1,NSCAN2
            READ(20) (T(J),J=1,NCHAN)
            IF(IFSCAN(I))THEN
              DO J=1,NCHAN
                TEMP(J)=TEMP(J)+T(J)
              END DO
            END IF
          END DO
          DO J=1,NCHAN
            TEMP(J)=TEMP(J)/REAL(NS1)
          END DO
        ELSE
          READ(20) (TEMP(J),J=1,NCHAN)
        END IF
        CLOSE(20)
C
        IF(CERR.EQ.'y')THEN
          WRITE(*,100)'Template error file name '
          CALL GUESSEF(TEMPFILE,ERRTEMPFILE)
          ERRTEMPFILE=INFILEX(21,ERRTEMPFILE,NSCAN2,NCHAN2,STWV2,DISP2,
     +     21,.TRUE.)                                                 !match
          IF(NSCAN2.GT.1)THEN
            DO J=1,NCHAN
              ERRTEMP(J)=0.
            END DO
            DO I=1,NSCAN2
              READ(21) (ET(J),J=1,NCHAN)
              IF(IFSCAN(I))THEN
                DO J=1,NCHAN
                  ERRTEMP(J)=ERRTEMP(J)+ET(J)*ET(J)
                END DO
              END IF
            END DO
            DO J=1,NCHAN
              ERRTEMP(J)=SQRT(ERRTEMP(J))/REAL(NS1)
            END DO
          ELSE
            READ(21) (ERRTEMP(J),J=1,NCHAN)
          END IF
          CLOSE(21)
        END IF
C
        WRITE(*,100)'Radial velocity (km/sec)'
        RVTEMP=READF('@')
C------------------------------------------------------------------------------
C calibracion en flujo
        WRITE(*,100)'Are you flux calibrating the data (y/n) '
        CFLUX(1:1)=READC('y','yn')
        IF(CFLUX.EQ.'y')THEN
          WRITE(*,100)'Response curve file name '
          FLUXFILE=INFILEX(20,'../stand/curvresf',
     +     NSCAN2,NCHAN2,STWV2,DISP2,1,.FALSE.)
          IF(NSCAN2.GT.1)THEN
            WRITE(*,101)'WARNING: this file contains more than 1 '//
     +       'spectrum.'
            WRITE(*,101)'The first spectrum will be used as '//
     +       'response curve.'
          END IF
          IF(NCHAN2.NE.NCHAN)THEN
            CLOSE(20)
            WRITE(*,100)'FATAL ERROR: NCHAN in last image is different!'
            STOP
          END IF
          IF(STWV2.NE.STWV)THEN
            CLOSE(20)
            WRITE(*,100)'FATAL ERROR: STWV in last image is different!'
            STOP
          END IF
          IF(DISP2.NE.DISP)THEN
            CLOSE(20)
            WRITE(*,100)'FATAL ERROR: DISP in last image is different!'
            STOP
          END IF
          READ(20) (FLUX(J),J=1,NCHAN)
          CLOSE(20)
        ELSE
          DO J=1,NCHAN
            FLUX(J)=1.
          END DO
        END IF
C------------------------------------------------------------------------------
        WRITE(*,100)'Create log file (y/n) '
        CSAVETXT(1:1)=READC('n','yn')
        IF(CSAVETXT.EQ.'y')THEN
          WRITE(*,100)'Log file name'
          TEXTFILE=OUTFILEX(40,'@',0,0,0.,0.,3,.FALSE.)
          WRITE(40,150)
          WRITE(40,101)'> This file is....: '//
     +     TEXTFILE(1:TRUELEN(TEXTFILE))
          WRITE(40,101)'> Input file......: '//
     +     INFILE(1:TRUELEN(INFILE))
          IF(CERR.EQ.'y')THEN
            WRITE(40,101)'> error file......: '//
     +       ERRINFILE(1:TRUELEN(ERRINFILE))
          ELSE
            WRITE(40,101)'> error file......: [none]'
          END IF
          WRITE(40,100)'> Radial velocity (km/sec): '
          WRITE(40,*) RVGAL
          WRITE(40,101)'> Template file...: '//
     +     TEMPFILE(1:TRUELEN(TEMPFILE))
          IF(CERR.EQ.'y')THEN
            WRITE(40,101)'> error file......: '//
     +       ERRTEMPFILE(1:TRUELEN(ERRTEMPFILE))
          ELSE
            WRITE(40,101)'> error file......: [none]'
          END IF
          WRITE(40,100)'> Radial velocity (km/sec): '
          WRITE(40,*) RVTEMP
          IF(CFLUX.EQ.'y')THEN
            WRITE(40,101)'> Flux file.......: '//
     +       FLUXFILE(1:TRUELEN(FLUXFILE))
          ELSE
            WRITE(40,101)'> (no flux file)'
          END IF
          WRITE(40,150)
        END IF
C------------------------------------------------------------------------------
C calibracion en flujo de la template
ccc5       DO J=1,NCHAN
        DO J=1,NCHAN
          T(J)=TEMP(J)/FLUX(J)
        END DO
        IF(CERR.EQ.'y')THEN
          DO J=1,NCHAN
            ET(J)=ERRTEMP(J)/FLUX(J)
          END DO
        END IF
C------------------------------------------------------------------------------
C colocamos la template al "redshift" de la galaxia
C NOTA: es muy importante tener en cuenta que la velocidad que hay que aplicar
C para colocar la template al mismo redshift que la galaxia NO se obtiene
C directamente a partir de la diferencia de velocidades. Es necesario corregir
C dicha diferencia de un factor que depende de la velocidad inicial de la
C template. Si no se realiza dicha correccion, los dos espectros se superponen
C de manera inexacta, aunque a simple vista el efecto puede pasar
C desapercibido. Consideramos tambien la correccion relativista.
C
        WRITE(*,100)'>>> Velocity difference between galaxy '//
     +   'and template: '
        WRITE(*,*) RVGAL-RVTEMP
        IF(CSAVETXT.EQ.'y')THEN
          WRITE(40,100)'>>> Velocity difference between galaxy '//
     +     'and template: '
          WRITE(40,*) RVGAL-RVTEMP
        END IF
C
        RVEFF=(RVGAL-RVTEMP)/(1.+RVTEMP/C)
        WRITE(*,100)'>>> Effective velocity difference'//
     +   ' (no relativ.)....: '
        WRITE(*,*)RVEFF
        IF(CSAVETXT.EQ.'y')THEN
          WRITE(40,100)'>>> Effective velocity difference'//
     +     ' (no relativ.)....: '
          WRITE(40,*)RVEFF
        END IF
C
        GAMMATEMP=SQRT(1.-(RVTEMP*RVTEMP)/(C*C))
        GAMMAGAL=SQRT(1.-(RVGAL*RVGAL)/(C*C))
        CTE=(1.+RVGAL/C)/(1.+RVTEMP/C)*GAMMATEMP/GAMMAGAL
        RVEFF=(CTE*CTE-1.)/(CTE*CTE+1.)*C
        WRITE(*,100)'>>> Effective velocity difference'//
     +   ' (relativ.).......: '
        WRITE(*,*)RVEFF
        IF(CSAVETXT.EQ.'y')THEN
          WRITE(40,100)'>>> Effective velocity difference'//
     +     ' (relativ.).......: '
          WRITE(40,*)RVEFF
        END IF
C
        CALL RVREBIN(RVEFF,NCHAN,T,TT,STWV,DISP)
        IF(CERR.EQ.'y') CALL RVREBIN(RVEFF,NCHAN,ET,ETT,STWV,DISP)
C------------------------------------------------------------------------------
C calculamos el intervalo disponible en longitud de onda
        WLMIN=STWV-DISP/2.
        WLMAX=STWV+REAL(NCHAN-1)*DISP+DISP/2.0
        WRITE(*,100)'>>> WLMIN..: '
        WRITE(*,*) WLMIN
        WRITE(*,100)'>>> WLMAX..: '
        WRITE(*,*) WLMAX
        IF(CSAVETXT.EQ.'y')THEN
          WRITE(40,100)'>>> WLMIN..: '
          WRITE(40,*) WLMIN
          WRITE(40,100)'>>> WLMAX..: '
          WRITE(40,*) WLMAX
        END IF
        RCVEL1=1.0+RVGAL/C
        RCVEL1=RCVEL1/SQRT(1.-(RVGAL*RVGAL)/(C*C))      !correccion relativista
C------------------------------------------------------------------------------
ccc7       CALL BUTTON(1,'[z]oom (m)',0)
        CALL BUTTON(1,'[z]oom (m)',0)
        CALL BUTTON(2,'zoom (k)',0)
        CALL BUTTON(3,'[w]hole',0)
        CALL BUTTON(4,'[j]ump',0)
        IF(NSCAN.EQ.1) CALL BUTTON(4,'[j]ump',3)
        CALL BUTTON(5,'[p]revious',0)
        IF(NSCAN.EQ.1) CALL BUTTON(5,'[p]revious',3)
        CALL BUTTON(6,'[n]ext',0)
        IF(NSCAN.EQ.1) CALL BUTTON(6,'[n]ext',3)
        CALL BUTTON(7,'[r]egions',0)
        CALL BUTTON(8,'[f]it',0)
        CALL BUTTON(9,'interpol[.]',0)
        CALL BUTTON(10,'[y] limits',0)
        CALL BUTTON(11,'[m]easure',0)
        CALL BUTTON(12,'resid[u]als',0)
        CALL BUTTON(13,'[i]ndices',0)
        CALL BUTTON(14,'[l]ines',0)
        CALL BUTTON(15,'[e]rrors',0)
        CALL BUTTON(16,'UBV [b]ands',0)
        CALL BUTTON(17,'[d]ef. reg.',0)
        CALL BUTTON(17,'[d]ef. reg.',3)
        IF(CERR.EQ.'n') CALL BUTTON(15,'[e]rrors',3)
        CALL BUTTON(18,'E[x]it',0)
C------------------------------------------------------------------------------
        XMIN=1.
        XMAX=REAL(NCHAN)
C..............................................................................
10      DO J=1,NCHAN
          SS(J)=A(J,NS0)/FLUX(J)
        END DO
        IF(CERR.EQ.'y')THEN
          DO J=1,NCHAN
            ESS(J)=ERRA(J,NS0)/FLUX(J)
          END DO
        END IF
20      SMEAN=0.
        TMEAN=0.
        NN=0
        IF(LREGIONS)THEN
          DO J=1,NCHAN
            IF(IFCHAN(J))THEN
              SMEAN=SMEAN+SS(J)
              TMEAN=TMEAN+TT(J)
              NN=NN+1
            END IF
          END DO
          WRITE(*,101)'* Normalization performed in previously '//
     +     'defined channel regions.'
        ELSE
          DO J=1,NCHAN
            IF((X(J).GE.XMIN).AND.(X(J).LE.XMAX))THEN
              SMEAN=SMEAN+SS(J)
              TMEAN=TMEAN+TT(J)
              NN=NN+1
            END IF
          END DO
          WRITE(*,101)'* Normalization performed in plot '//
     +     'region.'
        END IF
        SMEAN=SMEAN/REAL(NN)
        TMEAN=TMEAN/REAL(NN)
        IF(SMEAN.EQ.0.) SMEAN=1.0
        IF(TMEAN.EQ.0.) TMEAN=1.0
        DO J=1,NCHAN
          S(J)=SS(J)/SMEAN
        END DO
        IF(CERR.EQ.'y')THEN
          DO J=1,NCHAN
            ES(J)=ESS(J)/SMEAN
          END DO
        END IF
        IF(LFIT)THEN
          IF(CTFIT.EQ.'1')THEN              !polinomio ajustado a la diferencia
            DO J=1,NCHAN
              T(J)=TT(J)/TMEAN+P(J)
            END DO
          ELSEIF(CTFIT.EQ.'3')THEN
            DO J=1,NCHAN
              T(J)=K1*TT(J)/TMEAN+P(J)
            END DO
          ELSEIF(CTFIT.EQ.'2')THEN                                   !power law
            DO J=1,NCHAN
              T(J)=K1*TT(J)/TMEAN+Q(J)
            END DO
          END IF
        ELSE
          DO J=1,NCHAN
            T(J)=TT(J)/TMEAN
          END DO
        END IF
        IF(CERR.EQ.'y')THEN
          DO J=1,NCHAN
            ET(J)=ETT(J)/TMEAN
          END DO
        END IF
        DO J=1,NCHAN
          D(J)=S(J)-T(J)
        END DO
        IF(CERR.EQ.'y')THEN
          DO J=1,NCHAN
            ED(J)=SQRT(ES(J)*ES(J)+ET(J)*ET(J))
          END DO
        END IF
25      YMIN=1.E20
        YMAX=-1.E20
        IF(LPERROR)THEN
          DO J=1,NCHAN
            IF((X(J).GE.XMIN).AND.(X(J).LE.XMAX))THEN
              IF(S(J)-ES(J).LT.YMIN) YMIN=S(J)-ES(J)
              IF(S(J)+ES(J).GT.YMAX) YMAX=S(J)+ES(J)
              IF(T(J)-ET(J).LT.YMIN) YMIN=T(J)-ET(J)
              IF(T(J)+ET(J).GT.YMAX) YMAX=T(J)+ET(J)
              IF(LPRESIDUALS)THEN
                IF(D(J)-ED(J).LT.YMIN) YMIN=D(J)-ED(J)
                IF(D(J)+ED(J).GT.YMAX) YMAX=D(J)+ED(J)
              END IF
            END IF
          END DO
        ELSE
          DO J=1,NCHAN
            IF((X(J).GE.XMIN).AND.(X(J).LE.XMAX))THEN
              IF(S(J).LT.YMIN) YMIN=S(J)
              IF(S(J).GT.YMAX) YMAX=S(J)
              IF(T(J).LT.YMIN) YMIN=T(J)
              IF(T(J).GT.YMAX) YMAX=T(J)
              IF(LPRESIDUALS)THEN
                IF(D(J).LT.YMIN) YMIN=D(J)
                IF(D(J).GT.YMAX) YMAX=D(J)
              END IF
            END IF
          END DO
        END IF
        IF(LYMINMAX)THEN
          YMIN=YMIN0
          YMAX=YMAX0
        END IF
        DY=YMAX-YMIN
        YMIN=YMIN-DY/20.
        IF(.NOT.LYMINMAX)THEN
          YMIN=YMIN-DY/20.
          YMAX=YMAX+DY/20.
        END IF
30      XMINL=(XMIN-1.)*DISP+STWV
        XMAXL=(XMAX-1.)*DISP+STWV
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(ITERM.EQ.1)THEN
            CALL RPGENV(XMIN,XMAX,YMIN,YMAX,0,-2)
          ELSE
            CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,-2)
          END IF
          CALL PGIDEN_RED
          CALL PGSWIN(XMINL,XMAXL,YMIN,YMAX)
          CALL PGBOX('CMTSI',0.,0,'BCNTS',0.,0)
          CALL PGSWIN(XMIN,XMAX,YMIN,YMAX)
          CALL PGBOX('BNTSI',0.,0,' ',0.,0)
          CALL PGBIN(NCHAN,X,S,.TRUE.)
          IF(LPERROR)THEN
            IF(LCOLOR(ITERM)) CALL PGSCI(12)
            DO J=1,NCHAN
              IF((X(J).GE.XMIN).AND.(X(J).LE.XMAX))THEN
                CALL PGMOVE(X(J),S(J)-ES(J))
                CALL PGDRAW(X(J),S(J)+ES(J))
                CALL PGMOVE(X(J)-.25,S(J)-ES(J))
                CALL PGDRAW(X(J)+.25,S(J)-ES(J))
                CALL PGMOVE(X(J)-.25,S(J)+ES(J))
                CALL PGDRAW(X(J)+.25,S(J)+ES(J))
                CALL PGMOVE(X(J),T(J)-ET(J))
                CALL PGDRAW(X(J),T(J)+ET(J))
                CALL PGMOVE(X(J)-.25,T(J)-ET(J))
                CALL PGDRAW(X(J)+.25,T(J)-ET(J))
                CALL PGMOVE(X(J)-.25,T(J)+ET(J))
                CALL PGDRAW(X(J)+.25,T(J)+ET(J))
                IF(LPRESIDUALS)THEN
                  CALL PGMOVE(X(J),D(J)-ED(J))
                  CALL PGDRAW(X(J),D(J)+ED(J))
                  CALL PGMOVE(X(J)-.25,D(J)-ED(J))
                  CALL PGDRAW(X(J)+.25,D(J)-ED(J))
                  CALL PGMOVE(X(J)-.25,D(J)+ED(J))
                  CALL PGDRAW(X(J)+.25,D(J)+ED(J))
                END IF
              END IF
            END DO
          END IF
          IF(LCOLOR(ITERM)) CALL PGSCI(2)
          CALL PGBIN(NCHAN,X,T,.TRUE.)
          IF(LPRESIDUALS)THEN
            IF(LCOLOR(ITERM)) CALL PGSCI(3)
            CALL PGBIN(NCHAN,X,D,.TRUE.)
          END IF
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
          CALL PGLABEL('channel','Normalized no. counts',CHAR(32))
          WRITE(CDUMMY,*) NS0
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          CALL PGMTEXT('T',2.5,.5,.5,'Scan #'//CDUMMY(1:L))
          L1=TRUELEN(INFILE)
          L2=TRUELEN(OBJECT)
          CALL PGMTEXT('T',2.5,0.,0.,'file: '//INFILE(1:L1)//' '//
     +     OBJECT(1:L2))
          L=TRUELEN(TEMPFILE)
          CALL PGMTEXT('T',2.5,1.,1.,'template: '//TEMPFILE(1:L))
C
        END DO
        IF(LPINDEX) CALL SUBPINDEX(.FALSE.,NINDEX)
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LPLINES) CALL SUBPLINES(LCOLOR(ITERM))
          IF(LPBANDS) CALL SUBPBANDS(LCOLOR(ITERM))
          IF(LPREGUS) CALL PLOTREGUSED(NCHAN,X,IFCHAN,LCOLOR(ITERM))
        END DO
C..............................................................................
50      WRITE(*,101)'* REMEMBER: <#> output displayed data into files'
        CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
        CALL CHLOWER(CH)
        CALL IFBUTTON(XC,YC,NB)
C
        NBLOCAL=INDEX('zkwjpnrf.ymuilebdx',CH)
        IF((NBLOCAL.NE.0).AND.(CH.NE.' '))THEN
          CALL BUTTQEX(NBLOCAL,LBEXIST)
          IF(LBEXIST) NB=NBLOCAL
        END IF
C..............................................................................
        IF(NB.EQ.0)THEN
c.........
          IF(CH.EQ.'#')THEN
            WRITE(*,101)'Write displayed data into:'
            WRITE(*,101)'(1) ascii file'
            WRITE(*,101)'(2) REDUCEME file'
            WRITE(*,101)'(0) EXIT'
            WRITE(*,100)'Option (0..2) '
            CASCII(1:1)=READC('0','012')
            IF(CASCII.EQ.'1')THEN
              WRITE(*,100)'Output file name'
              FILEASCII=OUTFILEX(70,'@',0,0,0.,0.,3,.FALSE.)
              IF(CERR.EQ.'y')THEN
                DO J=1,NCHAN
                  WRITE(70,*)STWV+(REAL(J)-1.)*DISP,
     +             (STWV+(REAL(J)-1.)*DISP)/(RCVEL1),S(J),T(J),D(J),
     +             ES(J),ET(J),ED(J)
                END DO
              ELSE
                DO J=1,NCHAN
                  WRITE(70,*)STWV+(REAL(J)-1.)*DISP,
     +             (STWV+(REAL(J)-1.)*DISP)/(RCVEL1),S(J),T(J),D(J)
                END DO
              END IF
              CLOSE(70)
            ELSEIF(CASCII.EQ.'2')THEN
              WRITE(*,100)'Output file name'
              FILEASCII=OUTFILEX(70,'@',3,NCHAN,STWV,DISP,1,.FALSE.)
              WRITE(70) (S(J),J=1,NCHAN)
              WRITE(70) (T(J),J=1,NCHAN)
              WRITE(70) (D(J),J=1,NCHAN)
              CLOSE(70)
              IF(CERR.EQ.'y')THEN
                CALL GUESSEF(FILEASCII,ERRFILEASCII)
                WRITE(*,100)'Output error file name'
                FILEASCII=OUTFILEX(70,ERRFILEASCII,3,NCHAN,STWV,DISP,
     +           1,.TRUE.)
                WRITE(70) (ES(J),J=1,NCHAN)
                WRITE(70) (ET(J),J=1,NCHAN)
                WRITE(70) (ED(J),J=1,NCHAN)
                CLOSE(70)
              END IF
            END IF
c.........
          ELSE
            WRITE(*,100)'Cursor at: '
            WRITE(CDUMMY,*) XC
            CALL RMBLANK(CDUMMY,CDUMMY,L)
            WRITE(*,100)CDUMMY(1:L)//','
            WRITE(CDUMMY,*) YC
            CALL RMBLANK(CDUMMY,CDUMMY,L)
            WRITE(*,100)' '//CDUMMY(1:L)//',  wavelength: '
            WRITE(CDUMMY,*) STWV+(XC-1.)*DISP
            CALL RMBLANK(CDUMMY,CDUMMY,L)
            WRITE(*,100) CDUMMY(1:L)//'  (z=0: '
            WRITE(CDUMMY,*) (STWV+(XC-1.)*DISP)/(RCVEL1)
            CALL RMBLANK(CDUMMY,CDUMMY,L)
            WRITE(*,101) CDUMMY(1:L)//')'
            CALL FINDNEARLINE((STWV+(XC-1.)*DISP)/RCVEL1,CLABEL0,WV0)
            WRITE(*,100)'> Nearest line is: '
            WRITE(*,100)CLABEL0(1:TRUELEN(CLABEL0))
            WRITE(*,*)WV0
c.........
          END IF
C..............................................................................
        ELSEIF(NB.EQ.1)THEN
          CALL BUTTON(1,'[z]oom (m)',5)
          WRITE(*,100)'Press mouse (limit #1)...'
          CALL RPGBAND(6,0,0.,0.,XC,YC,CH)
          XMIN=XC
          IF(XMIN.LT.1.) XMIN=1.
          IF(XMIN.GT.REAL(NCHAN)) XMIN=REAL(NCHAN)
          WRITE(*,100)'   cursor at '
          WRITE(*,*)XMIN
          WRITE(*,100)'Press mouse (limit #2)...'
          CALL RPGBAND(4,0,XMIN,0.,XC,YC,CH)
          XMAX=XC
          IF(XMAX.LT.1.) XMAX=1.
          IF(XMAX.GT.REAL(NCHAN)) XMAX=REAL(NCHAN)
          WRITE(*,100)'   cursor at '
          WRITE(*,*)XMAX
          IF(XMAX.LT.XMIN)THEN
            XC=XMIN
            XMIN=XMAX
            XMAX=XC
          END IF
          CALL BUTTON(1,'[z]oom (m)',0)
          CALL RPGERASW(0.,1.,0.,0.80)
          GOTO 20
C..............................................................................
        ELSEIF(NB.EQ.2)THEN
          CALL BUTTON(2,'zoom (k)',5)
          WRITE(*,100)'Xmin '
          WRITE(CDUMMY,*) XMIN
          XMIN=READF(CDUMMY)
          WRITE(*,100)'Xmax '
          WRITE(CDUMMY,*) XMAX
          XMAX=READF(CDUMMY)
          CALL BUTTON(2,'zoom (k)',0)
          CALL RPGERASW(0.,1.,0.,0.80)
          GOTO 20
C..............................................................................
        ELSEIF(NB.EQ.3)THEN
          CALL BUTTON(3,'[w]hole',5)
          XMIN=1.
          XMAX=REAL(NCHAN)
          CALL BUTTON(3,'[w]hole',0)
          CALL RPGERASW(0.,1.,0.,0.80)
          GOTO 20
C..............................................................................
        ELSEIF(NB.EQ.4)THEN
          CALL BUTTON(4,'[j]ump',5)
          WRITE(*,100)'Next scan to work with '
          NS0=READILIM('@',1,NSCAN)
          IF(LFIT)THEN
            CALL BUTTON(8,'[f]it',0)
            LFIT=.FALSE.
          END IF
          CALL BUTTON(4,'[j]ump',0)
          CALL RPGERASW(0.,1.,0.,0.80)
          GOTO 10
C..............................................................................
        ELSEIF(NB.EQ.5)THEN
          CALL BUTTON(5,'[p]revious',5)
          IF(NS0.GT.1) NS0=NS0-1
          IF(LFIT)THEN
            CALL BUTTON(8,'[f]it',0)
            LFIT=.FALSE.
          END IF
          CALL BUTTON(5,'[p]revious',0)
          CALL RPGERASW(0.,1.,0.,0.80)
          GOTO 10
C..............................................................................
        ELSEIF(NB.EQ.6)THEN
          CALL BUTTON(6,'[n]ext',5)
          IF(NS0.LT.NSCAN) NS0=NS0+1
          IF(LFIT)THEN
            CALL BUTTON(8,'[f]it',0)
            LFIT=.FALSE.
          END IF
          CALL BUTTON(6,'[n]ext',0)
          CALL RPGERASW(0.,1.,0.,0.80)
          GOTO 10
C..............................................................................
        ELSEIF(NB.EQ.7)THEN
          IF(LREGIONS)THEN
            WRITE(*,100)'[r]emove regions or [a]dd new regions (r/a) '
            COPC(1:1)=READC('a','ra')
            IF(COPC.EQ.'r')THEN
              CALL BUTTON(7,'[r]egions',0)
              CALL BUTTON(17,'[d]ef. reg.',3)
              LPREGUS=.FALSE.
              LREGIONS=.FALSE.
            ELSE
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                CALL PLOTREGUSED(NCHAN,X,IFCHAN,LCOLOR(ITERM))
              END DO
              CALL SUBREGIONS(IFCHAN,NCHAN,.TRUE.)
            END IF
          ELSE
            CALL BUTTON(7,'[r]egions',1)
            CALL BUTTON(17,'[d]ef. reg.',0)
            LREGIONS=.TRUE.
            CALL SUBREGIONS(IFCHAN,NCHAN,.FALSE.)
          END IF
          IF(LFIT)THEN
            CALL BUTTON(8,'[f]it',0)
            LFIT=.FALSE.
          END IF
          CALL RPGERASW(0.,1.,0.,0.80)
          GOTO 20
C..............................................................................
        ELSEIF(NB.EQ.8)THEN
          IF(LFIT)THEN
            CALL BUTTON(8,'[f]it',0)
            LFIT=.FALSE.
          ELSE
            CALL BUTTON(8,'[f]it',1)
            LFIT=.TRUE.
          END IF
C
          IF(.NOT.LFIT)THEN
            CALL RPGERASW(0.,1.,0.,0.80)
            GOTO 20
          END IF
C
          IF(LREGIONS)THEN
            WRITE(*,101)'* Fit will be performed in previously '//
     +       'defined channel regions.'
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              CALL PLOTREGUSED(NCHAN,X,IFCHAN,LCOLOR(ITERM))
            END DO
          ELSE
            WRITE(*,101)'* Fit will be performed in plot '//
     +       'region.'
          END IF
C
          WRITE(*,101)'(1) fit (polynomial) to residuals'
          WRITE(*,101)'(2) fit (power law)+(k*template) to object'
          WRITE(*,101)'(3) fit (polynomial)+(k*template) to object'
          WRITE(*,100)'Option (1/2/3) '
          CTFIT(1:1)=READC('1','123')
C-------->
          IF((CTFIT.EQ.'1').OR.(CTFIT.EQ.'3'))THEN         !ajustamos polinomio
            NFIT=0                !definimos los puntos a utilizar en el ajuste
            IF(LREGIONS)THEN
              DO J=1,NCHAN
                IF(IFCHAN(J))THEN
                  NFIT=NFIT+1
                  XFIT(NFIT)=X(J)
                  YFIT(NFIT)=D(J)
                  IFCHAN_MFIT(J)=.TRUE.
                ELSE
                  IFCHAN_MFIT(J)=.FALSE.
                END IF
              END DO
            ELSE
              DO J=1,NCHAN
                IF((X(J).GE.XMIN).AND.(X(J).LE.XMAX))THEN
                  NFIT=NFIT+1
                  XFIT(NFIT)=X(J)
                  YFIT(NFIT)=D(J)
                  IFCHAN_MFIT(J)=.TRUE.
                ELSE
                  IFCHAN_MFIT(J)=.FALSE.
                END IF
              END DO
            END IF
            CREP='y'
            NN=0
            NTERMS=3                              !polinomio inicial de grado 2
            DO WHILE(CREP.EQ.'y')
              IF(CTFIT.EQ.'1')THEN
                WRITE(*,100)'Polynomial degree '
                WRITE(CDUMMY,*) NTERMS-1
                NTERMS=READILIM(CDUMMY,0,19)
                NTERMS=NTERMS+1
              END IF
              CALL POLFIT(XFIT,YFIT,YFIT,NFIT,NTERMS,0,COEFF,CHISQR)
              DO J=1,NCHAN
                P(J)=COEFF(NTERMS)
                DO K=NTERMS-1,1,-1
                  P(J)=P(J)*X(J)+COEFF(K)
                END DO
              END DO
              NN=NN+1
              IF(NN.GT.10) NN=1
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                IF(LCOLOR(ITERM)) CALL PGSCI(4+NN)
                CALL PGLINE(NCHAN,X,P)
                IF(LCOLOR(ITERM)) CALL PGSCI(1)
              END DO
              IF(CTFIT.EQ.'1')THEN
                WRITE(*,100)'Repeat fit (y/n) '
                CREP(1:1)=READC('n','yn')
              ELSE
                NDIM=4
                WRITE(CDUMMY,*) YRMSTOL
                WRITE(*,100)'YRMSTOL for DOWNHILL '
                YRMSTOL=READF(CDUMMY)
                XX0(1)=1.0
                XX0(2)=COEFF(1)
                XX0(3)=COEFF(2)
                XX0(4)=COEFF(3)
                DXX0(1)=0.1
                DO I=1,3
                  IF(DXX0(I+1).NE.0.0)THEN
                    DXX0(I+1)=COEFF(1)*0.05
                  ELSE
                    DXX0(I+1)=0.001 !que remedio
                  END IF
                END DO
                CALL DOWNHILL(NDIM,XX0,DXX0,FUNKPOLT,1.0,0.5,2.0,
     +           YRMSTOL,XXF,DXXF,NEVAL)
                K1=XXF(1)
                COEFF(1)=XXF(2)
                COEFF(2)=XXF(3)
                COEFF(3)=XXF(4)
                DO J=1,NCHAN
                  Q(J)=K1*T(J)+COEFF(1)+COEFF(2)*X(J)+COEFF(3)*X(J)*X(J)
                END DO
                NN=NN+1
                IF(NN.GT.10) NN=1
                DO ITERM=NTERM,1,-1
                  CALL PGSLCT(IDN(ITERM))
                  IF(LCOLOR(ITERM)) CALL PGSCI(4+NN)
                  CALL PGBIN(NCHAN,X,Q,.TRUE.)
                END DO
                DO J=1,NCHAN
                  Q(J)=K1*T(J)
                END DO
                DO ITERM=NTERM,1,-1
                  CALL PGSLCT(IDN(ITERM))
                  CALL PGBIN(NCHAN,X,Q,.TRUE.)
                END DO
                DO J=1,NCHAN
                  P(J)=COEFF(1)+COEFF(2)*X(J)+COEFF(3)*X(J)*X(J)
                END DO
                DO ITERM=NTERM,1,-1
                  CALL PGSLCT(IDN(ITERM))
                  CALL PGBIN(NCHAN,X,P,.TRUE.)
                  IF(LCOLOR(ITERM)) CALL PGSCI(1)
                END DO
                WRITE(*,100)'>>> k1...:'
                WRITE(*,*) K1
                WRITE(*,100)'>>> a0...:'
                WRITE(*,*) COEFF(1)
                WRITE(*,100)'>>> a1...:'
                WRITE(*,*) COEFF(2)
                WRITE(*,100)'>>> a2...:'
                WRITE(*,*) COEFF(3)
                CREP='n'
                WRITE(*,100)'Press <CR> to continue...'
                READ(*,*)
              END IF
            END DO
C-------->
          ELSEIF(CTFIT.EQ.'2')THEN   !ajustamos power law+template con DOWNHILL
            NDIM=3
            MP=NDIM+1
            JP=NDIM
            CREP='y'
            NN=0
            K1=1.0
            K2=0.1*(STWV+REAL(NCHAN)/2.*DISP)
            ALPHA=1.
            DO WHILE(CREP.EQ.'y')
              WRITE(CDUMMY,*) YRMSTOL
              WRITE(*,100)'YRMSTOL for DOWNHILL '
              YRMSTOL=READF(CDUMMY)
              WRITE(*,101)'* Enter initial parameters to start '//
     +         'DOWNHILL:'
              WRITE(CDUMMY,*) K1
              WRITE(*,100)'k1 '
              XX0(1)=READF(CDUMMY)
              WRITE(CDUMMY,*) 0.05*K1
              WRITE(*,100)'length-scale for K1 (DOWNHILL) '
              DXX0(1)=READF(CDUMMY)
              WRITE(CDUMMY,*) K2
              WRITE(*,100)'k2 '
              XX0(2)=READF(CDUMMY)
              WRITE(CDUMMY,*) 0.05*K2
              WRITE(*,100)'length-scale for K2 (DOWNHILL) '
              DXX0(2)=READF(CDUMMY)
              WRITE(CDUMMY,*) ALPHA
              WRITE(*,100)'alpha '
              XX0(3)=READF(CDUMMY)
              WRITE(CDUMMY,*) 0.05*ALPHA
              WRITE(*,100)'length-scale for alpha (DOWNHILL) '
              DXX0(3)=READF(CDUMMY)
              CALL DOWNHILL(NDIM,XX0,DXX0,FUNKPL,1.0,0.5,2.0,YRMSTOL,
     +         XXF,DXXF,NEVAL)
              K1=XXF(1)
              K2=XXF(2)
              ALPHA=XXF(3)
              DO J=1,NCHAN
                Q(J)=K1*T(J)+K2*(STWV+REAL(J-1)*DISP)**(ALPHA-2.0)
              END DO
              NN=NN+1
              IF(NN.GT.10) NN=1
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                IF(LCOLOR(ITERM)) CALL PGSCI(4+NN)
                CALL PGBIN(NCHAN,X,Q,.TRUE.)          !k1*template+k2*power_law
              END DO
              DO J=1,NCHAN
                Q(J)=K1*T(J)
              END DO
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                CALL PGBIN(NCHAN,X,Q,.TRUE.)                       !k1*template
              END DO
              DO J=1,NCHAN
                Q(J)=K2*(STWV+REAL(J-1)*DISP)**(ALPHA-2.0)
              END DO
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                CALL PGBIN(NCHAN,X,Q,.TRUE.)                      !k2*power_law
                IF(LCOLOR(ITERM)) CALL PGSCI(1)
              END DO
              WRITE(*,100)'>>> k1...:'
              WRITE(*,*) K1
              WRITE(*,100)'>>> k2...:'
              WRITE(*,*) K2
              WRITE(*,100)'>>> alpha: '
              WRITE(*,*) ALPHA
C
              CALL SUBPHOT(K1,K2,ALPHA,'U')
              CALL SUBPHOT(K1,K2,ALPHA,'B')
              CALL SUBPHOT(K1,K2,ALPHA,'V')
C
              WRITE(*,100)'Repeat fit (y/n) '
              CREP(1:1)=READC('n','yn')
            END DO
C-------->
          END IF
          CALL RPGERASW(0.,1.,0.,0.80)
          GOTO 20
C..............................................................................
        ELSEIF(NB.EQ.9)THEN
          CALL BUTTON(9,'interpol[.]',5)
          WRITE(*,101)'(1) Replace object data by fit to polynomial'
          WRITE(*,101)'(2) Replace object data by current template data'
          WRITE(*,101)'(3) Replace object data by multifit'
          WRITE(*,100)'Option (1/2/3) '
          CINTER0(1:1)=READC(CINTER0,'123')
          IF(CSAVETXT.EQ.'y')THEN
            WRITE(40,150)
            WRITE(40,100)'>>> Current Scan is #'
            WRITE(40,*) NS0
            WRITE(40,101)'Interpolation:'
            IF(CINTER0.EQ.'1')THEN
              WRITE(40,101)'(1) Replace object data by fit to '//
     +         'polynomial'
            ELSEIF(CINTER0.EQ.'2')THEN
              WRITE(40,101)'(2) Replace object data by current '//
     +         'template data'
            ELSEIF(CINTER0.EQ.'3')THEN
              WRITE(40,101)'(3) Replace object data by multifit'
            END IF
          END IF
          IF(CINTER0.EQ.CINTER)THEN
            WRITE(*,100)'Are you using the same regions (y/n) '
            CSAME(1:1)=READC('y','yn')
            IF(CSAVETXT.EQ.'y')THEN
              WRITE(40,101)'Are you using the same regions (y/n)? '//
     +         CSAME
            END IF
          ELSE
            CSAME='n'
          END IF
          IF(CSAVETXT.EQ.'y')THEN
            IF(CSAME.EQ.'n')THEN
              WRITE(*,101)'* Enter comment (optional) for log file '//
     +         '(max. 80 char.):'
              READ(*,'(A)')COMENTARIO
              WRITE(40,101) COMENTARIO
            ELSE
              WRITE(40,101) COMENTARIO
            END IF
          END IF
          CINTER=CINTER0
C-------->
          IF(CINTER.EQ.'1')THEN
            IF(CSAME.EQ.'n')THEN
              WRITE(*,101)'* Enter channels to fit polynomial:'
              IF(CSAVETXT.EQ.'y')WRITE(40,101)'* Enter channels '//
     +         'to fit polynomial:'
              CALL SUBREGIONS(IFCHAN_POLYN,NCHAN,.FALSE.)
            END IF
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              IF(LCOLOR(ITERM)) CALL PGSCI(7)
              DO J=1,NCHAN
                IF(IFCHAN_POLYN(J))THEN
                  CALL PGMOVE(REAL(J)-0.5,S(J))
                  CALL PGDRAW(REAL(J)+0.5,S(J))
                END IF
              END DO
              IF(LCOLOR(ITERM)) CALL PGSCI(1)
            END DO
            IF(CSAME.EQ.'n')THEN
              WRITE(*,101)'* Enter channels to be replaced by '//
     +         'polynomial:'
              IF(CSAVETXT.EQ.'y') WRITE(40,101)'* Enter channels '//
     +         'to be replaced by polynomial:'
              CALL SUBREGIONS(IFCHAN_INTER,NCHAN,.FALSE.)
            END IF
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              IF(LCOLOR(ITERM)) CALL PGSCI(11)
              DO J=1,NCHAN
                IF(IFCHAN_INTER(J))THEN
                  CALL PGMOVE(REAL(J)-0.5,S(J))
                  CALL PGDRAW(REAL(J)+0.5,S(J))
                END IF
              END DO
              IF(LCOLOR(ITERM)) CALL PGSCI(1)
            END DO
            WRITE(*,100)'Pol. degree '
            NTERMS=READILIM('0',0,19)
            IF(CSAVETXT.EQ.'y')THEN
              WRITE(40,100)'Pol. degree:'
              WRITE(40,*)NTERMS
            END IF
            NTERMS=NTERMS+1
            NFIT=0
            DO J=1,NCHAN
              IF(IFCHAN_POLYN(J))THEN
                NFIT=NFIT+1
                XFIT(NFIT)=X(J)
                YFIT(NFIT)=S(J)
              END IF
            END DO
            CALL POLFIT(XFIT,YFIT,YFIT,NFIT,NTERMS,0,COEFF,CHISQR)
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              DO J=1,NCHAN
                IF(IFCHAN_INTER(J))THEN
                  POL=COEFF(NTERMS)
                  DO K=NTERMS-1,1,-1
                    POL=POL*X(J)+COEFF(K)
                  END DO
                  CALL PGMOVE(REAL(J)-0.5,POL)
                  CALL PGDRAW(REAL(J)+0.5,POL)
                END IF
              END DO
            END DO
            WRITE(*,100)'Is this interpolation OK (y/n) '
            COK(1:1)=READC('y','yn')
            IF(CSAVETXT.EQ.'y') WRITE(40,101)'Is this '//
     +       'interpolation OK (y/n)? '//COK
            IF(COK.EQ.'y')THEN
              DO J=1,NCHAN
                IF(IFCHAN_INTER(J))THEN
                  POL=COEFF(NTERMS)
                  DO K=NTERMS-1,1,-1
                    POL=POL*X(J)+COEFF(K)
                  END DO
                  A(J,NS0)=POL*SMEAN*FLUX(J)
                END IF
              END DO
            END IF
C-------->
          ELSEIF(CINTER.EQ.'2')THEN
            IF(CSAME.EQ.'n')THEN
              WRITE(*,101)'* Enter channels to be replaced by template:'
              IF(CSAVETXT.EQ.'y') WRITE(40,101)'* Enter '//
     +         'channels to be replaced by template:'
              CALL SUBREGIONS(IFCHAN_INTER,NCHAN,.FALSE.)
            END IF
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              DO J=1,NCHAN
                IF(IFCHAN_INTER(J))THEN
                  IF(LCOLOR(ITERM)) CALL PGSCI(11)
                  CALL PGMOVE(REAL(J)-0.5,S(J))
                  CALL PGDRAW(REAL(J)+0.5,S(J))
                  IF(LCOLOR(ITERM)) CALL PGSCI(1)
                  CALL PGMOVE(REAL(J)-0.5,T(J))
                  CALL PGDRAW(REAL(J)+0.5,T(J))
                END IF
              END DO
            END DO
            WRITE(*,100)'Is this interpolation OK (y/n) '
            COK(1:1)=READC('y','yn')
            IF(CSAVETXT.EQ.'y') WRITE(40,101)'Is this '//
     +       'interpolation OK (y/n)? '//COK
            IF(COK.EQ.'y')THEN
              DO J=1,NCHAN
                IF(IFCHAN_INTER(J))THEN
                  A(J,NS0)=T(J)*SMEAN*FLUX(J)
                END IF
              END DO
            END IF
C-------->
          ELSEIF(CINTER.EQ.'3')THEN
            IF(LFIT)THEN
              WRITE(*,101)'ERROR: you must deactivate [f]it '//
     +         'button first.'
              IF(CSAVETXT.EQ.'y') WRITE(40,101)'ERROR: you must '//
     +         'deactivate [f]it button first.'
            END IF
            IF(LREGIONS)THEN
              WRITE(*,101)'ERROR: you must deactivate [r]egions '//
     +         'button first.'
              IF(CSAVETXT.EQ.'y') WRITE(40,101)'ERROR: you must '//
     +         'deactivate [r]egions button first.'
            END IF
            IF(LFIT.OR.LREGIONS)THEN
              WRITE(*,100)'Press <CR> to continue...'
              READ(*,*)
            ELSE
              NDIM=7
              IF(CSAME.EQ.'n')THEN
                WRITE(*,101)'Enter channels to fit data: '
                IF(CSAVETXT.EQ.'y') WRITE(40,101)'Enter '//
     +           'channels to fit data: '
                CALL SUBREGIONS(IFCHAN_MFIT,NCHAN,.FALSE.)
                WRITE(*,100)'Central wavelength of the emission line'
                WV0MFIT=READF('@')
                IF(CSAVETXT.EQ.'y') WRITE(40,100)'Central '//
     +           'wavelength of the emission line: '
                WRITE(40,*) WV0MFIT
                XX0(1)=1.0                                             !k1
                XX0(3)=(WV0MFIT-STWV)/DISP+1.0                         !x0
                XX0(2)=S(NINT(XX0(3)))-T(NINT(XX0(3)))                 !a
                XX0(4)=2.0                                             !sigma
                XX0(5)=0.01                                            !a0
                XX0(6)=0.0001                                          !a1
                XX0(7)=0.000001                                        !a2
                WRITE(CDUMMY,*) XX0(1)
                WRITE(*,100)'k1    '
                XX0(1)=READF(CDUMMY)
                WRITE(CDUMMY,*)0.1
                WRITE(*,100)'length-scale (DOWNHILL) '
                DXX0(1)=READF(CDUMMY)
                WRITE(CDUMMY,*) XX0(2)
                WRITE(*,100)'a     '
                XX0(2)=READF(CDUMMY)
                WRITE(CDUMMY,*)0.05*XX0(2)
                WRITE(*,100)'length-scale (DOWNHILL) '
                DXX0(2)=READF(CDUMMY)
                WRITE(CDUMMY,*) XX0(03)
                WRITE(*,100)'x0    '
                XX0(3)=READF(CDUMMY)
                WRITE(CDUMMY,*)0.05*XX0(3)
                WRITE(*,100)'length-scale (DOWNHILL) '
                DXX0(3)=READF(CDUMMY)
                WRITE(CDUMMY,*) XX0(4)
                WRITE(*,100)'sigma '
                XX0(4)=READF(CDUMMY)
                WRITE(CDUMMY,*)0.05*XX0(4)
                WRITE(*,100)'length-scale (DOWNHILL) '
                DXX0(4)=READF(CDUMMY)
                WRITE(CDUMMY,*) XX0(5)
                WRITE(*,100)'a0    '
                XX0(5)=READF(CDUMMY)
                WRITE(CDUMMY,*)0.05*XX0(5)
                WRITE(*,100)'length-scale (DOWNHILL) '
                DXX0(5)=READF(CDUMMY)
                WRITE(CDUMMY,*) XX0(6)
                WRITE(*,100)'a1    '
                XX0(6)=READF(CDUMMY)
                WRITE(CDUMMY,*)0.05*XX0(6)
                WRITE(*,100)'length-scale (DOWNHILL) '
                DXX0(6)=READF(CDUMMY)
                WRITE(CDUMMY,*) XX0(7)
                WRITE(*,100)'a2    '
                XX0(7)=READF(CDUMMY)
                WRITE(CDUMMY,*)0.05*XX0(7)
                WRITE(*,100)'length-scale (DOWNHILL) '
                DXX0(7)=READF(CDUMMY)
                WRITE(CDUMMY,*) YRMSTOL
                WRITE(*,100)'YRMSTOL for DOWNHILL '
                YRMSTOL=READF(CDUMMY)
              ELSE
                DO J=1,NDIM
                  XX0(J)=XXOLD(J)
                  DXX0(J)=DXXOLD(J)
                END DO
              END IF
              WRITE(*,100)'Please wait (running DOWNHILL)...'
              CALL DOWNHILL(NDIM,XX0,DXX0,FUNKMPFIT,1.0,0.5,2.0,YRMSTOL,
     +         XXF,DXXF,NEVAL)
              WRITE(*,*)
              WRITE(*,100)'k1    '
              WRITE(*,*)XXF(1)
              WRITE(*,100)'a     '
              WRITE(*,*)XXF(2)
              WRITE(*,100)'x0    '
              WRITE(*,*)XXF(3),(XXF(3)-1.)*DISP+STWV
              WRITE(*,100)'sigma '
              WRITE(*,*)XXF(4)
              WRITE(*,100)'a0    '
              WRITE(*,*)XXF(5)
              WRITE(*,100)'a1    '
              WRITE(*,*)XXF(6)
              WRITE(*,100)'a2    '
              WRITE(*,*)XXF(7)
              IF(CSAME.EQ.'n')THEN
                DO J=1,JP
                  XXOLD(J)=XXF(J)
                  DXXOLD(J)=DXXF(J)
                END DO
              END IF
              IF(CSAVETXT.EQ.'y')THEN
                WRITE(40,100)'k1    '
                WRITE(40,*)XXF(1)
                WRITE(40,100)'a     '
                WRITE(40,*)XXF(2)
                WRITE(40,100)'x0    '
                WRITE(40,*)XXF(3),(XXF(3)-1.)*DISP+STWV
                WRITE(40,100)'sigma '
                WRITE(40,*)XXF(4)
                WRITE(40,100)'a0    '
                WRITE(40,*)XXF(5)
                WRITE(40,100)'a1    '
                WRITE(40,*)XXF(6)
                WRITE(40,100)'a2    '
                WRITE(40,*)XXF(7)
              END IF
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                DO J=1,NCHAN
                  IF(IFCHAN_MFIT(J))THEN
                    POL=XXF(1)*T(J)+
     +               XXF(5)+XXF(6)*REAL(J)+XXF(7)*REAL(J)*REAL(J)
                    IF(LCOLOR(ITERM)) CALL PGSCI(5)
                    CALL PGMOVE(REAL(J)-0.5,POL)
                    CALL PGDRAW(REAL(J)+0.5,POL)
                    POL=XXF(1)*T(J)
                    IF(LCOLOR(ITERM)) CALL PGSCI(6)
                    CALL PGMOVE(REAL(J)-0.5,POL)
                    CALL PGDRAW(REAL(J)+0.5,POL)
                    POL=XXF(2)*EXP(-(REAL(J)-XXF(3))*(REAL(J)-XXF(3))/
     +                (2.*XXF(4)*XXF(4)))
                    IF(LCOLOR(ITERM)) CALL PGSCI(7)
                    CALL PGPOINT(1,REAL(J),POL,1)
                    POL=XXF(5)+XXF(6)*REAL(J)+XXF(7)*REAL(J)*REAL(J)
                    IF(LCOLOR(ITERM)) CALL PGSCI(8)
                    CALL PGMOVE(REAL(J)-0.5,POL)
                    CALL PGDRAW(REAL(J)+0.5,POL)
                  END IF
                END DO
                IF(LCOLOR(ITERM)) CALL PGSCI(1)
              END DO
C
              WRITE(*,100)'Write ascii file with fitted data (y/n) '
              CASCII(1:1)=READC('n','yn')
              IF(CASCII.EQ.'y')THEN
                WRITE(*,100)'Output file name'
                FILEASCII=OUTFILEX(70,'@',0,0,0.,0.,3,.FALSE.)
                DO J=1,NCHAN
                  IF(IFCHAN_MFIT(J))THEN
                    POL1=XXF(1)*T(J)+
     +               XXF(5)+XXF(6)*REAL(J)+XXF(7)*REAL(J)*REAL(J)
                    POL2=XXF(1)*T(J)
                    POL3=XXF(2)*EXP(-(REAL(J)-XXF(3))*(REAL(J)-XXF(3))/
     +                (2.*XXF(4)*XXF(4)))
                    POL4=XXF(5)+XXF(6)*REAL(J)+XXF(7)*REAL(J)*REAL(J)
                    WRITE(70,*)STWV+(REAL(J)-1.)*DISP,
     +               (STWV+(REAL(J)-1.)*DISP)/(RCVEL1),S(J),T(J),
     +               POL1,POL2,POL3,POL4
                  ELSE
                    WRITE(70,*)STWV+(REAL(J)-1.)*DISP,
     +               (STWV+(REAL(J)-1.)*DISP)/(RCVEL1),S(J),T(J),
     +               -99.,-99.,-99.,-99.
                  END IF
                END DO
                CLOSE(70)
              END IF
C
              WRITE(*,100)'Is this fit OK (y/n) '
              COK(1:1)=READC('y','yn')
              IF(CSAVETXT.EQ.'y') WRITE(40,101)'Is this fit '//
     +         'OK (y/n)? '//COK
              IF(COK.EQ.'y')THEN
                WRITE(*,101)'* Enter channels to be replaced by fit:'
                IF(CSAVETXT.EQ.'y') WRITE(40,101)'* Enter '//
     +           'channels to be replaced by fit:'
                CALL SUBREGIONS(IFCHAN_INTER,NCHAN,.FALSE.)
                DO ITERM=NTERM,1,-1
                  CALL PGSLCT(IDN(ITERM))
                  IF(LCOLOR(ITERM)) CALL PGSCI(11)
                  DO J=1,NCHAN
                    IF(IFCHAN_INTER(J))THEN
                      POL=XXF(1)*T(J)+
     +                 XXF(5)+XXF(6)*REAL(J)+XXF(7)*REAL(J)*REAL(J)
                      CALL PGMOVE(REAL(J)-0.5,POL)
                      CALL PGDRAW(REAL(J)+0.5,POL)
                    END IF
                  END DO
                  IF(LCOLOR(ITERM)) CALL PGSCI(1)
                END DO
                WRITE(*,100)'Is this interpolation OK (y/n) '
                COK(1:1)=READC('y','yn')
                IF(CSAVETXT.EQ.'y') WRITE(40,101)'Is this '//
     +           'interpolation OK (y/n)? '//COK
                IF(COK.EQ.'y')THEN
                  DO J=1,NCHAN
                    IF(IFCHAN_INTER(J))THEN
                      POL=XXF(1)*T(J)+
     +                 XXF(5)+XXF(6)*REAL(J)+XXF(7)*REAL(J)*REAL(J)
                      A(J,NS0)=POL*SMEAN*FLUX(J)
                    END IF
                  END DO
                END IF
              END IF
            END IF
          END IF
C
          CALL BUTTON(9,'interpol[.]',0)
          CALL RPGERASW(0.,1.,0.,0.80)
          GOTO 10
C..............................................................................
        ELSEIF(NB.EQ.10)THEN
          IF(LYMINMAX)THEN
            CALL BUTTON(10,'[y] limits',0)
            LYMINMAX=.FALSE.
            CALL RPGERASW(0.,1.,0.,0.80)
            GOTO 25
          ELSE
            CALL BUTTON(10,'[y] limits',1)
            LYMINMAX=.TRUE.
            WRITE(*,100)'Ymin '
            WRITE(CDUMMY,*) YMIN
            YMIN0=READF(CDUMMY)
            YMIN=YMIN0
            WRITE(*,100)'Ymax '
            WRITE(CDUMMY,*) YMAX
            YMAX0=READF(CDUMMY)
            YMAX=YMAX0
            CALL RPGERASW(0.,1.,0.,0.80)
            GOTO 30
          END IF
C..............................................................................
        ELSEIF(NB.EQ.11)THEN
          CALL BUTTON(11,'[m]easure',5)
          WRITE(*,101)'(1) object-template'
          WRITE(*,101)'(2) object'
          WRITE(*,101)'(3) template'
          IF(LPRESIDUALS)THEN
            WRITE(*,101)'(4) residuals'
            WRITE(*,100)'Option (1..4) '
            CMEAS(1:1)=READC('1','1234')
          ELSE
            WRITE(*,100)'Option (1..3) '
            CMEAS(1:1)=READC('1','123')
          END IF
          IF(CMEAS.EQ.'1')THEN
            CALL MIDEWIDTH(NCHAN,S,ES,T,ET,CERR,CMEAS,NITER,YRMSTOL)
          ELSEIF(CMEAS.EQ.'2')THEN
            CALL MIDEWIDTH(NCHAN,S,ES,S,ES,CERR,CMEAS,NITER,YRMSTOL)
          ELSEIF(CMEAS.EQ.'3')THEN
            CALL MIDEWIDTH(NCHAN,T,ET,T,ET,CERR,CMEAS,NITER,YRMSTOL)
          ELSEIF(CMEAS.EQ.'4')THEN
            CALL MIDEWIDTH(NCHAN,D,ED,D,ED,CERR,CMEAS,NITER,YRMSTOL)
          END IF
          CALL BUTTON(11,'[m]easure',0)
C..............................................................................
        ELSEIF(NB.EQ.12)THEN
          IF(LPRESIDUALS)THEN
            CALL BUTTON(12,'resid[u]als',0)
            LPRESIDUALS=.FALSE.
          ELSE
            CALL BUTTON(12,'resid[u]als',1)
            LPRESIDUALS=.TRUE.
          END IF
          CALL RPGERASW(0.,1.,0.,0.80)
          GOTO 10
C..............................................................................
        ELSEIF(NB.EQ.13)THEN
          IF(LPINDEX)THEN
            CALL BUTTON(13,'[i]ndices',0)
            LPINDEX=.FALSE.
          ELSE
            CALL BUTTON(13,'[i]ndices',1)
            LPINDEX=.TRUE.
            CALL SUBPINDEX(.TRUE.,NINDEX)
          END IF
C..............................................................................
        ELSEIF(NB.EQ.14)THEN
          IF(LPLINES)THEN
            CALL BUTTON(14,'[l]ines',0)
            LPLINES=.FALSE.
          ELSE
            CALL BUTTON(14,'[l]ines',1)
            LPLINES=.TRUE.
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              CALL SUBPLINES(LCOLOR(ITERM))
            END DO
          END IF
C..............................................................................
        ELSEIF(NB.EQ.15)THEN
          IF(LPERROR)THEN
            CALL BUTTON(15,'[e]rrors',0)
            LPERROR=.FALSE.
          ELSE
            CALL BUTTON(15,'[e]rrors',1)
            LPERROR=.TRUE.
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              IF(LCOLOR(ITERM)) CALL PGSCI(12)
              DO J=1,NCHAN
                IF((X(J).GE.XMIN).AND.(X(J).LE.XMAX))THEN
                  CALL PGMOVE(X(J),S(J)-ES(J))
                  CALL PGDRAW(X(J),S(J)+ES(J))
                  CALL PGMOVE(X(J)-.25,S(J)-ES(J))
                  CALL PGDRAW(X(J)+.25,S(J)-ES(J))
                  CALL PGMOVE(X(J)-.25,S(J)+ES(J))
                  CALL PGDRAW(X(J)+.25,S(J)+ES(J))
                  CALL PGMOVE(X(J),T(J)-ET(J))
                  CALL PGDRAW(X(J),T(J)+ET(J))
                  CALL PGMOVE(X(J)-.25,T(J)-ET(J))
                  CALL PGDRAW(X(J)+.25,T(J)-ET(J))
                  CALL PGMOVE(X(J)-.25,T(J)+ET(J))
                  CALL PGDRAW(X(J)+.25,T(J)+ET(J))
                  IF(LPRESIDUALS)THEN
                    CALL PGMOVE(X(J),D(J)-ED(J))
                    CALL PGDRAW(X(J),D(J)+ED(J))
                    CALL PGMOVE(X(J)-.25,D(J)-ED(J))
                    CALL PGDRAW(X(J)+.25,D(J)-ED(J))
                    CALL PGMOVE(X(J)-.25,D(J)+ED(J))
                    CALL PGDRAW(X(J)+.25,D(J)+ED(J))
                  END IF
                END IF
              END DO
              IF(LCOLOR(ITERM)) CALL PGSCI(1)
            END DO
          END IF
C..............................................................................
        ELSEIF(NB.EQ.16)THEN
          IF(LPBANDS)THEN
            CALL BUTTON(16,'UBV [b]ands',0)
            LPBANDS=.FALSE.
          ELSE
            CALL BUTTON(16,'UBV [b]ands',1)
            LPBANDS=.TRUE.
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              CALL SUBPBANDS(LCOLOR(ITERM))
            END DO
          END IF
C..............................................................................
        ELSEIF(NB.EQ.17)THEN
          IF(LPREGUS)THEN
            CALL BUTTON(17,'[d]ef. reg.',0)
            LPREGUS=.FALSE.
          ELSE
            CALL BUTTON(17,'[d]ef. reg.',1)
            LPREGUS=.TRUE.
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              CALL PLOTREGUSED(NCHAN,X,IFCHAN,LCOLOR(ITERM))
            END DO
          END IF
C..............................................................................
        ELSEIF(NB.EQ.18)THEN
          CALL BUTTON(18,'E[x]it',5)
          WRITE(*,100)'Do you really want to exit (y/n) '
          CSURE(1:1)=READC('n','yn')
          IF(CSURE.EQ.'n')THEN
            CALL BUTTON(18,'E[x]it',0)
          ELSE
            GOTO 70
          END IF
C..............................................................................
        END IF
        GOTO 50
C------------------------------------------------------------------------------
C Salvamos fichero(s) de salida
70      WRITE(*,100)'Save output file name..... (y/n) '
        CSAVE(1:1)=READC('n','yn')
        IF(CSAVE.EQ.'y')THEN
          OBJECT=OBJECT0
          FITSFILE=FITSFILE0
          COMMENT=COMMENT0
          AIRMASS=AIRMASS0
          TIMEXPOS=TIMEXPOS0
          WRITE(*,100)'Output file name'
          OUTFILE=OUTFILEX(30,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
          IF(CSAVETXT.EQ.'y')THEN
            WRITE(40,150)
            WRITE(40,101)'Output file name: '//
     +       OUTFILE(1:TRUELEN(OUTFILE))
          END IF
          DO I=1,NSCAN
            WRITE(30) (A(J,I),J=1,NCHAN)
          END DO
          CLOSE(30)
        END IF
C------------------------------------------------------------------------------
C cerramos salida grafica
        CALL PGEND
C cerramos log file si es necesario
        IF(CSAVETXT.EQ.'y')THEN
          WRITE(40,150)
          CLOSE(40)
        END IF
C------------------------------------------------------------------------------
        STOP
100     FORMAT(A,$)
101     FORMAT(A)
150     FORMAT(79('-'))
        END
C
C******************************************************************************
C
C Dibuja las bandas de los indices seleccionados
C Si LASK=.TRUE., la subrutina pregunta por la(s) banda(s) a dibujar
C NINDEX devuelve el indice dibujado
        SUBROUTINE SUBPINDEX(LASK,NINDEX)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INTEGER READILIM
C
        LOGICAL LASK
        INTEGER NINDEX
C
        INTEGER K,IWL
        INTEGER NI,NINDEXT,ITI,NIND1,NIND2
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        INTEGER NCONTI,NABSOR
        REAL RCVEL1,WLMIN,WLMAX
        REAL YMIN,YMAX
        REAL X1,X2
        REAL WV(NWVMAX),WLMIN0,WLMAX0,FWV(NWVMAX/4)
        REAL DY
        CHARACTER*8 CLABEL
        LOGICAL LANY,LINDOK(NINDMAX)
        LOGICAL LCOLOR(MAX_ID_RED)
C
        COMMON/BLKSHOW/LINDOK
        COMMON/BLKGEN1/STWV,DISP
        COMMON/BLKGEN2/RCVEL1,WLMIN,WLMAX
        COMMON/BLKGEN3/YMIN,YMAX
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
        NCONTI=0 !evita warning de compilacion
        NABSOR=0 !evita warning de compilacion
C
        DY=YMAX-YMIN
        CALL SELINDEX(0,WV,FWV,NINDEXT,CLABEL)
        IF(LASK)THEN
          LANY=.FALSE.
          DO NINDEX=1,NINDEXT
            CALL SELINDEX(NINDEX,WV,FWV,ITI,CLABEL)
            IF((ITI.EQ.1).OR.(ITI.EQ.2))THEN
              DO K=1,6
                WV(K)=WV(K)*RCVEL1
              END DO
            ELSEIF(ITI.EQ.3)THEN
              DO K=1,4
                WV(K)=WV(K)*RCVEL1
              END DO
            ELSEIF(ITI.EQ.4)THEN
              DO K=1,4
                WV(K)=WV(K)*RCVEL1
              END DO
            ELSEIF((ITI.GE.101).OR.(ITI.LE.9999))THEN
              NCONTI=(ITI/100)                  !numero de regiones de continuo
              NABSOR=ITI-NCONTI*100         !numero de regiones con absorciones
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
            ELSEIF((ITI.EQ.3).OR.(ITI.EQ.4))THEN
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
        ELSE
          LANY=.TRUE.
        END IF
C
        IF(LANY)THEN
          IF(LASK)THEN
            CALL SHINDEX(LINDOK,0)
            WRITE(*,100)'Index bands to be plotted'
            NINDEX=READILIM('@',-1,NINDEXT)
          END IF
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
                ELSEIF((ITI.EQ.3).OR.(ITI.EQ.4))THEN
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
        SUBROUTINE SUBPBANDS(LCOLOR)
        IMPLICIT NONE
C
        INCLUDE 'redlib.inc'
        LOGICAL LCOLOR
C
        INTEGER NPBAND,I
        REAL WV(NPBANDMAX),RES(NPBANDMAX)
        REAL RCVEL1,WLMIN,WLMAX
        REAL XX(NPBANDMAX)
        REAL YMIN,YMAX
C
        COMMON/BLKGEN1/STWV,DISP
        COMMON/BLKGEN2/RCVEL1,WLMIN,WLMAX
        COMMON/BLKGEN3/YMIN,YMAX
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
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
        SUBROUTINE SUBPLINES(LCOLOR)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        LOGICAL LCOLOR
C
        INTEGER NLINE,NLINEST
        REAL RCVEL1,WLMIN,WLMAX
        REAL YMIN,YMAX
        REAL X0,X0MIN,X0MAX,Y0MIN,Y0MAX
        REAL WV(NLINMAX)
        REAL DY
        CHARACTER*8 CLABEL(NLINMAX)
C
        COMMON/BLKGEN1/STWV,DISP
        COMMON/BLKGEN2/RCVEL1,WLMIN,WLMAX
        COMMON/BLKGEN3/YMIN,YMAX
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        CALL PGQWIN(X0MIN,X0MAX,Y0MIN,Y0MAX)
        DY=YMAX-YMIN
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
ccc          WV(NLINE)=WV(NLINE)*RCVEL1
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
C Define que canales van a ser utilizados para normalizar los espectros.
C Si LADD=.TRUE. la subrutina no inicializa el array IFCHAN, pudiendose asi
C mantener regiones previamente seleccionadas
C
        SUBROUTINE SUBREGIONS(IFCHAN,NCHAN,LADD)
        IMPLICIT NONE
        INCLUDE 'futils.inc'
C
        INTEGER NCHAN
        LOGICAL IFCHAN(NCHAN),LADD
C
        INTEGER J
        INTEGER NC1,NC2,NC0
        REAL XC,YC
        REAL YMIN,YMAX,DY
        CHARACTER*1 CMODE,CH,CSAVETXT
C
        COMMON/BLKMWIDTH0/CSAVETXT
        COMMON/BLKGEN3/YMIN,YMAX
C------------------------------------------------------------------------------
        DY=YMAX-YMIN
C
        IF(.NOT.LADD)THEN
          DO J=1,NCHAN
            IFCHAN(J)=.FALSE.
          END DO
        END IF
C
        WRITE(*,100)'Select region by [m]ouse or [k]eyboard (m/k) '
        CMODE(1:1)=READC('m','mk')
        IF(CSAVETXT.EQ.'y') WRITE(40,101)'Select region by [m]ouse '//
     +   'or [k]eyboard (m/k)? '//CMODE
10      IF(CMODE.EQ.'m')THEN
          WRITE(*,101)'Press <q>/<x>/<mouse right button> to EXIT'
          WRITE(*,100)'Press mouse (limit #1)...'
          IF(CSAVETXT.EQ.'y')THEN
            WRITE(40,101)'Press <q>/<x>/<mouse right button> to EXIT'
            WRITE(40,100)'Press mouse (limit #1)...'
          END IF
          CALL RPGBAND(6,0,0.,0.,XC,YC,CH)
          CALL CHLOWER(CH)
          IF((CH.EQ.'x').OR.(CH.EQ.'q'))THEN
            WRITE(*,*)
            IF(CSAVETXT.EQ.'y') WRITE(40,*)
            GOTO 20
          END IF
          NC1=NINT(XC)
          IF(NC1.LT.1) NC1=1
          IF(NC1.GT.NCHAN) NC1=NCHAN
          WRITE(*,100)'   cursor at '
          WRITE(*,*)NC1
          IF(CSAVETXT.EQ.'y')THEN
            WRITE(40,100)'   cursor at '
            WRITE(40,*)NC1
          END IF
          CALL PGMOVE(REAL(NC1),YMIN)
          CALL PGDRAW(REAL(NC1),YMAX)
          WRITE(*,100)'Press mouse (limit #2)...'
          IF(CSAVETXT.EQ.'y')THEN
            WRITE(40,100)'Press mouse (limit #2)...'
          END IF
          CALL RPGBAND(6,0,0.,0.,XC,YC,CH)
          CALL CHLOWER(CH)
          IF((CH.EQ.'x').OR.(CH.EQ.'q'))THEN
            WRITE(*,*)
            IF(CSAVETXT.EQ.'y') WRITE(40,*)
            GOTO 20
          END IF
          NC2=NINT(XC)
          IF(NC2.LT.1) NC2=1
          IF(NC2.GT.NCHAN) NC2=NCHAN
          WRITE(*,100)'   cursor at '
          WRITE(*,*)NC2
          IF(CSAVETXT.EQ.'y')THEN
            WRITE(40,100)'   cursor at '
            WRITE(40,*)NC2
          END IF
          CALL PGMOVE(REAL(NC2),YMIN)
          CALL PGDRAW(REAL(NC2),YMAX)
          IF(NC1.GT.NC2)THEN
            NC0=NC1
            NC1=NC2
            NC2=NC0
          END IF
          CALL PGRECT(REAL(NC1),REAL(NC2),YMIN,YMIN+DY/100.)
          CALL PGRECT(REAL(NC1),REAL(NC2),YMAX,YMAX-DY/100.)
        ELSE
12        WRITE(*,100)'Channel region (0,0=EXIT) '
          CALL READ2I('0,0',NC1,NC2)
          IF(CSAVETXT.EQ.'y')THEN
            WRITE(40,100)'Channel region (0,0=EXIT):'
            WRITE(40,*)NC1,NC2
          END IF
          IF((NC1.EQ.0).AND.(NC2.EQ.0)) GOTO 20
          IF((NC1.LT.1).OR.(NC2.GT.NCHAN).OR.(NC1.GT.NC2))THEN
            WRITE(*,101)'ERROR: invalid numbers. Try again.'
            IF(CSAVETXT.EQ.'y') WRITE(40,101)'ERROR: invalid '//
     +       'numbers. Try again.'
            GOTO 12
          END IF
        END IF
        DO J=NC1,NC2
          IFCHAN(J)=.TRUE.
        END DO
        GOTO 10
C
20      NC0=0
        DO J=1,NCHAN
          IF(IFCHAN(J)) NC0=NC0+1
        END DO
        IF(NC0.EQ.0)THEN
          WRITE(*,101)'ERROR: number of channels to be used = 0! '//
     +     'Try again.'
          IF(CSAVETXT.EQ.'y') WRITE(40,101)'ERROR: number of '//
     +     'channels to be used = 0!  Try again.'
          GOTO 10
        END IF
C
100     FORMAT(A,$)
101     FORMAT(A)
        END
C
C******************************************************************************
C Dibuja las regiones que seran utilizadas para los ajustes
        SUBROUTINE PLOTREGUSED(NCHAN,X,IFCHAN,LCOLOR)
        IMPLICIT NONE
C
        INTEGER NCHAN
        REAL X(NCHAN)
        LOGICAL IFCHAN(NCHAN)
        LOGICAL LCOLOR
C
        INTEGER J,J1
        REAL YMIN,YMAX
        REAL DY
        LOGICAL LEXIT
C
        COMMON/BLKGEN3/YMIN,YMAX
C------------------------------------------------------------------------------
        J1=0 !avoid compilation warning
        DY=YMAX-YMIN
C
        IF(LCOLOR) CALL PGSCI(7)
        J=0
        DO WHILE(J.LT.NCHAN)
          LEXIT=.FALSE.
          J=J+1
          DO WHILE(.NOT.LEXIT)
            IF(IFCHAN(J))THEN
              CALL PGMOVE(X(J),YMIN+DY/40.)
              CALL PGDRAW(X(J),YMAX-DY/40.)
              J1=J
              IF(J.EQ.NCHAN)THEN
                CALL PGMOVE(X(J),YMIN+DY/40.)
                CALL PGDRAW(X(J),YMAX-DY/40.)
                CALL PGRECT(X(J1),X(J),YMIN+DY/40.,YMIN+DY/35.)
                CALL PGRECT(X(J1),X(J),YMAX-DY/40.,YMAX-DY/35.)
              END IF
              LEXIT=.TRUE.
            ELSE
              J=J+1
              IF(J.GT.NCHAN)THEN
                LEXIT=.TRUE.
              ELSE
                LEXIT=.FALSE.
              END IF
            END IF
          END DO
          IF(J.LT.NCHAN)THEN
            LEXIT=.FALSE.
            J=J+1
          ELSE
            LEXIT=.TRUE.
          END IF
          DO WHILE(.NOT.LEXIT)
            IF(IFCHAN(J))THEN
              IF(J.GT.NCHAN)THEN
                CALL PGMOVE(X(NCHAN),YMIN+DY/40.)
                CALL PGDRAW(X(NCHAN),YMAX-DY/40.)
                CALL PGRECT(X(J1),X(NCHAN),YMIN+DY/40.,YMIN+DY/35.)
                CALL PGRECT(X(J1),X(NCHAN),YMAX-DY/40.,YMAX-DY/35.)
                LEXIT=.TRUE.
              ELSE
                IF(J.LT.NCHAN)THEN
                  LEXIT=.FALSE.
                  J=J+1
                ELSE
                  CALL PGMOVE(X(NCHAN),YMIN+DY/40.)
                  CALL PGDRAW(X(NCHAN),YMAX-DY/40.)
                  CALL PGRECT(X(J1),X(NCHAN),YMIN+DY/40.,YMIN+DY/35.)
                  CALL PGRECT(X(J1),X(NCHAN),YMAX-DY/40.,YMAX-DY/35.)
                  LEXIT=.TRUE.
                END IF
              END IF
            ELSE
              CALL PGMOVE(X(J-1),YMIN+DY/40.)
              CALL PGDRAW(X(J-1),YMAX-DY/40.)
              CALL PGRECT(X(J1),X(J-1),YMIN+DY/40.,YMIN+DY/35.)
              CALL PGRECT(X(J1),X(J-1),YMAX-DY/40.,YMAX-DY/35.)
              LEXIT=.TRUE.
            END IF
          END DO
        END DO
        IF(LCOLOR) CALL PGSCI(1)
C
        END
C
C******************************************************************************
C Escribe las regiones en el log file
        SUBROUTINE WRITEREG(NCHAN,IFCHAN)
        IMPLICIT NONE
C
        INTEGER NCHAN
        LOGICAL IFCHAN(NCHAN)
C
        INTEGER J,J1
        INTEGER L
        CHARACTER*50 CDUMMY
        LOGICAL LEXIT
C------------------------------------------------------------------------------
        J=0
        DO WHILE(J.LT.NCHAN)
          LEXIT=.FALSE.
          J=J+1
          DO WHILE(.NOT.LEXIT)
            IF(IFCHAN(J))THEN
              J1=J
              IF(J.EQ.NCHAN)THEN
                WRITE(CDUMMY,200)'(',J1,',',J,')'
                CALL RMBLANK(CDUMMY,CDUMMY,L)
                WRITE(40,100)CDUMMY(1:L)//' '
              END IF
              LEXIT=.TRUE.
            ELSE
              J=J+1
              IF(J.GT.NCHAN)THEN
                LEXIT=.TRUE.
              ELSE
                LEXIT=.FALSE.
              END IF
            END IF
          END DO
          IF(J.LT.NCHAN)THEN
            LEXIT=.FALSE.
            J=J+1
          ELSE
            LEXIT=.TRUE.
          END IF
          DO WHILE(.NOT.LEXIT)
            IF(IFCHAN(J))THEN
              IF(J.GT.NCHAN)THEN
                WRITE(CDUMMY,200)'(',J1,',',NCHAN,')'
                CALL RMBLANK(CDUMMY,CDUMMY,L)
                WRITE(40,100)CDUMMY(1:L)//' '
                LEXIT=.TRUE.
              ELSE
                IF(J.LT.NCHAN)THEN
                  LEXIT=.FALSE.
                  J=J+1
                ELSE
                  WRITE(CDUMMY,200)'(',J1,',',NCHAN,')'
                  CALL RMBLANK(CDUMMY,CDUMMY,L)
                  WRITE(40,100)CDUMMY(1:L)//' '
                  LEXIT=.TRUE.
                END IF
              END IF
            ELSE
              WRITE(CDUMMY,200)'(',J1,',',J-1,')'
              CALL RMBLANK(CDUMMY,CDUMMY,L)
              WRITE(40,100)CDUMMY(1:L)//' '
              LEXIT=.TRUE.
            END IF
          END DO
        END DO
        WRITE(40,*)
C
100     FORMAT(A,$)
200     FORMAT(A1,I8,A1,I8,A1)
        END
C
C******************************************************************************
C Realiza una fotometria aproximada y comparativa entre el ajuste a al ley
C de potencias, la template y el espectro promedio.
C CBAND es la banda en la que realizamos la fotometria.
        SUBROUTINE SUBPHOT(K1,K2,ALPHA,CBAND)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        REAL K1,K2,ALPHA           !parametros del ajuste a la ley de potencias
        CHARACTER*1 CBAND
C
        REAL S(NCMAX),T(NCMAX)
        LOGICAL IFCHAN(NCMAX)
C
        INTEGER I,J,J1,J2
        INTEGER IMIN(NCMAX),IMAX(NCMAX)
        INTEGER NPBAND
        REAL YMIN,YMAX
        REAL X1,X2
        REAL WV(NPBANDMAX),RES(NPBANDMAX)
        REAL RCVEL1,WLMIN,WLMAX
        REAL WL0
        REAL RESCHAN(NCMAX)                      !curva respuesta en cada pixel
        REAL FLUXS,FLUXT,FLUXQ
        LOGICAL LOK
C
        COMMON/BLKGEN1/STWV,DISP
        COMMON/BLKGEN2/RCVEL1,WLMIN,WLMAX
        COMMON/BLKGEN3/YMIN,YMAX
        COMMON/BLKFUNK1/IFCHAN
        COMMON/BLKFUNK2/S,T
        COMMON/BLKFUNK3/NCHAN
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        CALL SELBANDS(CBAND,NPBAND,WV,RES)
        DO I=1,NPBAND
          WV(I)=WV(I)*RCVEL1
        END DO
        DO J=1,NCHAN
          IMIN(J)=0
          IMAX(J)=0
        END DO
        DO I=1,NPBAND-1
          X1=(WV(I)-STWV)/DISP+1.
          X2=(WV(I+1)-STWV)/DISP+1.
          J1=NINT(X1)
          J2=NINT(X2)
          LOK=.TRUE.
          IF(J1.LT.1) J1=MIN(1,J2)
          IF(J2.GT.NCHAN) J2=MAX(NCHAN,J1)
          IF(J2.LT.1) LOK=.FALSE.
          IF(J1.GT.NCHAN) LOK=.FALSE.
          IF(LOK)THEN
            DO J=J1,J2
              IMIN(J)=I
              IMAX(J)=I+1
            END DO
          END IF
        END DO
C
        DO J=1,NCHAN
          IF(IFCHAN(J))THEN
            WL0=STWV+REAL(J-1)*DISP
            IF((IMIN(J).NE.0).AND.(IMAX(J).NE.0))THEN
              RESCHAN(J)=RES(IMIN(J))+(RES(IMAX(J))-RES(IMIN(J)))*
     +         (WL0-WV(IMIN(J)))/(WV(IMAX(J))-WV(IMIN(J)))
              CALL PGPOINT(1,REAL(J),RESCHAN(J)*0.90*YMAX,-1)
            ELSE
              RESCHAN(J)=0.
            END IF
          ELSE
            RESCHAN(J)=0.
          END IF
        END DO
C
        FLUXS=0.
        FLUXT=0.
        FLUXQ=0.
        DO J=1,NCHAN
          FLUXS=FLUXS+S(J)*RESCHAN(J)
          FLUXT=FLUXT+K1*T(J)*RESCHAN(J)
          FLUXQ=FLUXQ+K2*RESCHAN(J)*(STWV+REAL(J-1)*DISP)**(ALPHA-2.0)
        END DO
C
        WRITE(*,100)'> '//CBAND//' flux in Object spectrum........: '
        WRITE(*,*)FLUXS
        WRITE(*,100)'> '//CBAND//' flux in k1*Template spectrum...: '
        WRITE(*,*)FLUXT
        WRITE(*,100)'> '//CBAND//' flux in k2*Power law spectrum..: '
        WRITE(*,*)FLUXQ
        IF(FLUXS.GT.0.)THEN
          WRITE(*,100)'> '//CBAND//' flux: (Object-Fit)/Object (%)..: '
          WRITE(*,*) 100.0*(FLUXS-FLUXT-FLUXQ)/FLUXS
          WRITE(*,100)'> '//CBAND//' flux: fraction (%) of template.: '
          WRITE(*,*) 100.0*FLUXT/FLUXS
          WRITE(*,100)'> '//CBAND//' flux: fraction (%) of power law: '
          WRITE(*,*) 100.0*FLUXQ/FLUXS
        END IF
        WRITE(*,*)
C
100     FORMAT(A,$)
C
        END
C
C******************************************************************************
C Ajustamos el espectro del objeto a un factor (k1) por el espectro de la
C template mas otro factor (k2) por una ley de potencias (alpha)
C XX(1)=k1, XX(2)=k2, XX(3)=alpha
        REAL FUNCTION FUNKPL(XX)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        REAL XX(3)
C
        INTEGER J
        REAL S(NCMAX),T(NCMAX)
        REAL FF
        LOGICAL IFCHAN(NCMAX)
C
        COMMON/BLKGEN1/STWV,DISP
        COMMON/BLKFUNK1/IFCHAN
        COMMON/BLKFUNK2/S,T
        COMMON/BLKFUNK3/NCHAN
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        IF(ABS(XX(3)).GT.8.)THEN
          FUNKPL=1.E30
          RETURN
        END IF
        FUNKPL=0.
        DO J=1,NCHAN
          IF(IFCHAN(J))THEN
            FF=XX(1)*T(J)+XX(2)*(STWV+REAL(J-1)*DISP)**(XX(3)-2.0)
            FUNKPL=FUNKPL+(S(J)-FF)*(S(J)-FF)
          END IF
        END DO
C
        END
C
C******************************************************************************
C Ajustamos el espectro del objeto a un factor (k1) por el espectro de la
C template mas una gaussiana mas un polinomio de grado 2
C XX(1)=k1, XX(2)=a, XX(3)=x0, XX(4)=sigma, XX(5)=a0, XX(6)=a1, XX(7)=a2
        REAL FUNCTION FUNKMPFIT(XX)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        REAL XX(7)
C
        INTEGER J
        REAL S(NCMAX),T(NCMAX)
        REAL FF
        LOGICAL IFCHAN_MFIT(NCMAX)
C
        COMMON/BLKFUNKMPFIT/IFCHAN_MFIT
        COMMON/BLKFUNK2/S,T
        COMMON/BLKFUNK3/NCHAN
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        FUNKMPFIT=0.
        DO J=1,NCHAN
          IF(IFCHAN_MFIT(J))THEN
            FF=XX(1)*T(J)+
     +   XX(2)*EXP(-(REAL(J)-XX(3))*(REAL(J)-XX(3))/(2.*XX(4)*XX(4)))+
     +   XX(5)+XX(6)*REAL(J)+XX(7)*REAL(J)*REAL(J)
            FUNKMPFIT=FUNKMPFIT+(S(J)-FF)*(S(J)-FF)
          END IF
        END DO
C
        END
C
C******************************************************************************
C Ajustamos el espectro del objeto a un factor (k1) por el espectro de la
C template mas un polinomio de segundo grado
C XX(1)=k1, XX(2)=a0, XX(3)=a1, XX(4)=a2
        REAL FUNCTION FUNKPOLT(XX)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        REAL XX(4)
C
        INTEGER J
        REAL S(NCMAX),T(NCMAX)
        REAL FF
        LOGICAL IFCHAN_MFIT(NCMAX)
C
        COMMON/BLKFUNKMPFIT/IFCHAN_MFIT
        COMMON/BLKFUNK2/S,T
        COMMON/BLKFUNK3/NCHAN
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        FUNKPOLT=0.
        DO J=1,NCHAN
          IF(IFCHAN_MFIT(J))THEN
            FF=XX(1)*T(J)+XX(2)+XX(3)*REAL(J)+XX(4)*REAL(J*J)
            FUNKPOLT=FUNKPOLT+(S(J)-FF)*(S(J)-FF)
          END IF
        END DO
C
        END
C
C******************************************************************************
C Mide las lineas
        SUBROUTINE MIDEWIDTH(NCHAN,Y,YE,YT,YTE,CERR,CMEAS,NITER,YRMSTOL)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'futils.inc'
        INTEGER TRUELEN
        INTEGER READILIM
        REAL READF
C
        REAL Y(NCHAN),YE(NCHAN)
        REAL YT(NCHAN),YTE(NCHAN)
        CHARACTER*1 CERR,CMEAS
        INTEGER NITER
        REAL YRMSTOL
C
        REAL PI
        PARAMETER(PI=3.141592654)
C
        INTEGER J,JMIN,JMAX,K
        INTEGER NFIT,NPLOT,NTERMS_CONT
        INTEGER NS0
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        REAL X(NCMAX)
        REAL XFIT(NCMAX),YFIT(NCMAX),EYFIT(NCMAX)
        REAL XPLOT(NCMAX),YPLOT(NCMAX)
        REAL A(20),CHISQR
        REAL XX,YY,POL,EPOL
        REAL WIDTH,ERRWIDTH
        REAL AMP,X0,SIGMA
        REAL EAMP,EX0,ESIGMA
        REAL EEAMP,EEX0,EESIGMA
        REAL RCVEL1,WLMIN,WLMAX,WV0
        REAL FLUX1,FLUX2
        REAL ERRFLUX1,ERRFLUX2
        REAL FWHM,FLUXG,EFLUXG
        CHARACTER*1 CBEFORE,CSAVETXT,CCONFIRM,CMEASLAST
        CHARACTER*50 CDUMMY
        CHARACTER*8 CLABEL0
        LOGICAL IFCHAN_CONT(NCMAX)
        LOGICAL IFCHAN_LINE(NCMAX)
        LOGICAL LBUFFER
        LOGICAL LCOLOR(MAX_ID_RED)
C
        COMMON/BLKGEN1/STWV,DISP
        COMMON/BLKGEN2/RCVEL1,WLMIN,WLMAX
        COMMON/BLKMWIDTH0/CSAVETXT
        COMMON/BLKMWIDTH1/LBUFFER
        COMMON/BLKMWIDTH2/IFCHAN_CONT,IFCHAN_LINE
        COMMON/BLKMWIDTH3/NTERMS_CONT
        COMMON/BLKMWIDTH4/NS0
        COMMON/BLKMWIDTH5/CMEASLAST
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
        JMIN=0
        JMAX=0
        EPOL=0.0
C
        IF(LBUFFER)THEN
          IF((CMEASLAST.EQ.CMEAS).OR.(CMEASLAST.NE.'1'))THEN
            WRITE(*,100)'Use previous region selection (y/n) '
            CBEFORE(1:1)=READC('n','yn')
          ELSE
            CBEFORE='n'
            LBUFFER=.TRUE.
          END IF
        ELSE
          CBEFORE='n'
          LBUFFER=.TRUE.
        END IF
        CMEASLAST=CMEAS
C
        IF(CBEFORE.EQ.'n')THEN
          DO J=1,NCHAN
            IFCHAN_CONT(J)=.FALSE.
            IFCHAN_LINE(J)=.FALSE.
          END DO
        END IF
C
        DO J=1,NCHAN
          X(J)=REAL(J)
        END DO
C
        IF(CMEAS.EQ.'1')GOTO 10
C
        IF(CBEFORE.EQ.'n')THEN
          WRITE(*,*)
          WRITE(*,101)'* Select channels to be employed to fit '//
     +     'continuum:'
          CALL SUBREGIONS(IFCHAN_CONT,NCHAN,.FALSE.)
        END IF
        NFIT=0
        IF(CERR.EQ.'y') EPOL=0.
        DO J=1,NCHAN
          IF(IFCHAN_CONT(J))THEN
            NFIT=NFIT+1
            XFIT(NFIT)=X(J)
            YFIT(NFIT)=Y(J)
            IF(CERR.EQ.'y') EPOL=EPOL+YE(J)
          END IF
        END DO
        IF(CERR.EQ.'y') EPOL=EPOL/REAL(NFIT)
        IF(CBEFORE.EQ.'n')THEN
          WRITE(*,100)'Pol. degree '
          IF(NFIT.GT.19)THEN
            NTERMS_CONT=READILIM('0',0,19)
          ELSE
            NTERMS_CONT=READILIM('0',0,NFIT-1)
          END IF
          NTERMS_CONT=NTERMS_CONT+1
        END IF
        CALL POLFIT(XFIT,YFIT,YFIT,NFIT,NTERMS_CONT,0,A,CHISQR)
C determinamos region para dibujar gaussiana
        DO J=1,NCHAN
          IF(IFCHAN_CONT(J))THEN
            JMAX=J
          END IF
        END DO
        DO J=NCHAN,1,-1
          IF(IFCHAN_CONT(J))THEN
            JMIN=J
          END IF
        END DO
C
        NPLOT=0
        DO J=JMIN,JMAX
          NPLOT=NPLOT+1
          XPLOT(NPLOT)=REAL(J)
          YPLOT(NPLOT)=A(NTERMS_CONT)
          DO K=NTERMS_CONT-1,1,-1
            YPLOT(NPLOT)=YPLOT(NPLOT)*XPLOT(NPLOT)+A(K)
          END DO
        END DO
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM)) CALL PGSCI(5)
          CALL PGBIN(NPLOT,XPLOT,YPLOT,.TRUE.)
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
        END DO
C
10      IF(CBEFORE.EQ.'n')THEN
          WRITE(*,*)
          WRITE(*,101)'* Select channels to be employed to measure '//
     +     'the line:'
          CALL SUBREGIONS(IFCHAN_LINE,NCHAN,.FALSE.)
        END IF
C determinamos region para dibujar gaussiana
        IF(CMEAS.EQ.'1')THEN
          DO J=1,NCHAN
            IF(IFCHAN_LINE(J))THEN
              JMAX=J
            END IF
          END DO
          DO J=NCHAN,1,-1
            IF(IFCHAN_LINE(J))THEN
              JMIN=J
            END IF
          END DO
        END IF
C
        NFIT=0
        DO J=1,NCHAN
          IF(IFCHAN_LINE(J))THEN
            NFIT=NFIT+1
            XFIT(NFIT)=X(J)
            YFIT(NFIT)=Y(J)
          END IF
        END DO
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM)) CALL PGSCI(7)
          CALL PGBIN(NFIT,XFIT,YFIT,.TRUE.)
        END DO
        FLUX1=0.
        FLUX2=0.
        WIDTH=0.
        IF(CERR.EQ.'y')THEN
          ERRFLUX1=0.
          ERRFLUX2=0.
          ERRWIDTH=0.
        END IF
        NFIT=0
        DO J=1,NCHAN
          IF(IFCHAN_LINE(J))THEN
            XX=REAL(J)
            IF(CMEAS.EQ.'1')THEN
              POL=YT(J)
              IF(CERR.EQ.'y') EPOL=YTE(J)
            ELSE
              POL=A(NTERMS_CONT)
              DO K=NTERMS_CONT-1,1,-1
                POL=POL*XX+A(K)
              END DO
            END IF
            NFIT=NFIT+1
            YFIT(NFIT)=Y(J)-POL
            FLUX1=FLUX1+(POL-Y(J))*DISP
            FLUX2=FLUX2+Y(J)*DISP
            WIDTH=WIDTH+DISP*(POL-Y(J))/POL
            IF(CERR.EQ.'y')THEN
              EYFIT(NFIT)=SQRT(YE(J)*YE(J)+EPOL*EPOL)
              ERRFLUX1=ERRFLUX1+(YE(J)*YE(J)+EPOL*EPOL)
              ERRFLUX2=ERRFLUX2+YE(J)*YE(J)
              ERRWIDTH=ERRWIDTH+(YE(J)*YE(J)+EPOL*EPOL)/(POL*POL)+
     +         ((POL-Y(J))/(POL*POL))*((POL-Y(J))/(POL*POL))*EPOL*EPOL
            END IF
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              CALL PGMOVE(XX,POL)
              CALL PGDRAW(XX,Y(J))
            END DO
          END IF
        END DO
        IF(CERR.EQ.'y')THEN
          ERRFLUX1=DISP*SQRT(ERRFLUX1)
          ERRFLUX2=DISP*SQRT(ERRFLUX2)
          ERRWIDTH=DISP*SQRT(ERRWIDTH)
        ELSE
          ERRFLUX1=0.
          ERRFLUX2=0.
          ERRWIDTH=0.
        END IF
        WRITE(CDUMMY,*)YRMSTOL
        WRITE(*,100)'YRMSTOL for DOWNHILL '
        YRMSTOL=READF(CDUMMY)
        IF(CERR.EQ.'y')THEN
          WRITE(CDUMMY,*) NITER
          WRITE(*,100)'No of iterations to evaluate errors '
          NITER=READILIM(CDUMMY,2,10000)
          CALL GAUSSFIT(NFIT,XFIT,YFIT,EYFIT,X0,SIGMA,AMP,
     +     EX0,ESIGMA,EAMP,EEX0,EESIGMA,EEAMP,YRMSTOL,NITER)
        ELSE
          CALL GAUSSFIT(NFIT,XFIT,YFIT,EYFIT,X0,SIGMA,AMP,
     +     EX0,ESIGMA,EAMP,EEX0,EESIGMA,EEAMP,YRMSTOL,0)
        END IF
C
        XX=X0    !calculamos el valor del continuo en el centro de la gaussiana
        IF(CMEAS.EQ.'1')THEN
          YY=YT(NINT(XX))
        ELSE
          YY=A(NTERMS_CONT)
          DO K=NTERMS_CONT-1,1,-1
            YY=YY*XX+A(K)
          END DO
        END IF
C
        WRITE(*,100)'>>> EW under gauss.(Angstroms)...'//
     +   '[object frame]: '
        WRITE(*,*) -DISP/YY*SQRT(2.*PI)*AMP*SIGMA,
     +   -DISP/YY*SQRT(2.*PI)*SQRT(SIGMA*SIGMA*EAMP*EAMP+
     +   AMP*AMP*ESIGMA*ESIGMA)
        WRITE(*,100)'>>> EW under spectrum (Angstroms)'//
     +   '[object frame]: '
        WRITE(*,*) WIDTH,ERRWIDTH
        WRITE(*,100)'>>> EW under gauss.(Angstroms).....'//
     +   '[rest frame]: '
        WRITE(*,*) -DISP/YY*SQRT(2.*PI)*AMP*SIGMA/RCVEL1,
     +   -DISP/YY*SQRT(2.*PI)*SQRT(SIGMA*SIGMA*EAMP*EAMP+
     +   AMP*AMP*ESIGMA*ESIGMA)/RCVEL1
        WRITE(*,100)'>>> EW under spectrum (Angstroms)..'//
     +   '[rest frame]: '
        WRITE(*,*) WIDTH/RCVEL1,ERRWIDTH/RCVEL1
        WRITE(*,100)'>>> Integrated Flux (spec.-cont.)...'//
     +   '[any frame]: '
        WRITE(*,*) FLUX1,ERRFLUX1
        WRITE(*,100)'>>> Integrated Flux (spectrum)......'//
     +   '[any frame]: '
        WRITE(*,*) FLUX2,ERRFLUX2
        FLUXG=-AMP*SQRT(2.*PI)*SIGMA*DISP
        EFLUXG=-SQRT(2.*PI)*DISP*SQRT(AMP*AMP*ESIGMA*ESIGMA+
     +   SIGMA*SIGMA*EAMP*EAMP)
        WRITE(*,100)'>>> Integrated Flux (gaussian)......'//
     +   '[any frame]: '
        WRITE(*,*) FLUXG,EFLUXG
        WRITE(*,*)
        FWHM=SIGMA*DISP*SQRT(-8.*ALOG(0.5))
        WRITE(*,100)'FWHM (Angstroms): '
        WRITE(*,*) FWHM,SQRT(-8.*ALOG(0.5))*DISP*ESIGMA
        WRITE(*,101)'Coefficients from fit: '//
     +   'y=a*exp[-((x-x0)**2)/(2 sig**2)]'
        WRITE(*,100)'> a   (+ error, rmsDOWNHILL): '
        WRITE(*,*)AMP,EAMP,EEAMP
        WRITE(*,100)'> x0  (+ error, rmsDOWNHILL): '
        WRITE(*,*)X0,EX0,EEX0
        WRITE(*,100)'> sig (+ error, rmsDOWNHILL): '
        WRITE(*,*)SIGMA,ESIGMA,EESIGMA
        WRITE(*,100)'> Central wavelength [object frame]: '
        WRITE(*,*)STWV+(X0-1.)*DISP
        WRITE(*,100)'> Central wavelength [rest frame]..: '
        WRITE(*,*)(STWV+(X0-1.)*DISP)/RCVEL1
        CALL FINDNEARLINE((STWV+(X0-1.)*DISP)/RCVEL1,CLABEL0,WV0)
        WRITE(*,100)'> Nearest line is: '//CLABEL0(1:TRUELEN(CLABEL0))
        WRITE(*,*)WV0
C
        NPLOT=0
        DO J=JMIN,JMAX
          NPLOT=NPLOT+1
          XPLOT(NPLOT)=REAL(J)
          IF(CMEAS.EQ.'1')THEN
            POL=YT(J)
          ELSE
            POL=A(NTERMS_CONT)
            DO K=NTERMS_CONT-1,1,-1
              POL=POL*XPLOT(NPLOT)+A(K)
            END DO
          END IF
          YPLOT(NPLOT)=POL+AMP*EXP(-(XPLOT(NPLOT)-X0)*
     +     (XPLOT(NPLOT)-X0)/(2.*SIGMA*SIGMA))
        END DO
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM)) CALL PGSCI(6)
          CALL PGLINE(NPLOT,XPLOT,YPLOT)
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
        END DO
C
        IF(CSAVETXT.EQ.'y')THEN
          WRITE(*,100)'Save results into log file (y/n) '
          CCONFIRM(1:1)=READC('y','yn')
          IF(CCONFIRM.EQ.'y')THEN
            WRITE(40,150)
            IF(CMEAS.EQ.'1')THEN
              WRITE(40,100)'> [object-template]'
            ELSEIF(CMEAS.EQ.'2')THEN
              WRITE(40,100)'> [object] '
            ELSEIF(CMEAS.EQ.'3')THEN
              WRITE(40,100)'> [template]'
            ELSEIF(CMEAS.EQ.'4')THEN
              WRITE(40,100)'> [residuals] '
            END IF
            WRITE(40,110)'  Scan #',NS0
            WRITE(40,100)'Continuum (channels): '
            CALL WRITEREG(NCHAN,IFCHAN_CONT)
            WRITE(40,100)'Line..... (channels): '
            CALL WRITEREG(NCHAN,IFCHAN_LINE)
            WRITE(40,100)'>>> EW under gauss.(Angstroms)...'//
     +       '[object frame]: '
            WRITE(40,*) -DISP/YY*SQRT(2.*PI)*AMP*SIGMA,
     +       -DISP/YY*SQRT(2.*PI)*SQRT(SIGMA*SIGMA*EAMP*EAMP+
     +       AMP*AMP*ESIGMA*ESIGMA)
            WRITE(40,100)'>>> EW under spectrum (Angstroms)'//
     +       '[object frame]: '
            WRITE(40,*) WIDTH,ERRWIDTH
            WRITE(40,100)'>>> EW under gauss.(Angstroms).....'//
     +       '[rest frame]: '
            WRITE(40,*) -DISP/YY*SQRT(2.*PI)*AMP*SIGMA/RCVEL1,
     +       -DISP/YY*SQRT(2.*PI)*SQRT(SIGMA*SIGMA*EAMP*EAMP+
     +       AMP*AMP*ESIGMA*ESIGMA)/RCVEL1
            WRITE(40,100)'>>> EW under spectrum (Angstroms)..'//
     +       '[rest frame]: '
            WRITE(40,*) WIDTH/RCVEL1,ERRWIDTH/RCVEL1
            WRITE(40,100)'>>> Integrated Flux (spec.-cont.)...'//
     +       '[any frame]: '
            WRITE(40,*) FLUX1,ERRFLUX1
            WRITE(40,100)'>>> Integrated Flux (spectrum)......'//
     +       '[any frame]: '
            WRITE(40,*) FLUX2,ERRFLUX2
            FLUXG=-AMP*SQRT(2.*PI)*SIGMA*DISP
            EFLUXG=-SQRT(2.*PI)*DISP*SQRT(AMP*AMP*ESIGMA*ESIGMA+
     +       SIGMA*SIGMA*EAMP*EAMP)
            WRITE(40,100)'>>> Integrated Flux (gaussian)......'//
     +       '[any frame]: '
            WRITE(40,*) FLUXG,EFLUXG
            WRITE(40,*)
            FWHM=SIGMA*DISP*SQRT(-8.*ALOG(0.5))
            WRITE(40,100)'FWHM (Angstroms): '
            WRITE(40,*) FWHM,SQRT(-8.*ALOG(0.5))*DISP*ESIGMA
            WRITE(40,101)'Coefficients from fit: '//
     +       'y=a*exp[-((x-x0)**2)/(2 sig**2)]'
            WRITE(40,100)'> a   (+ error, rmsDOWNHILL): '
            WRITE(40,*)AMP,EAMP,EEAMP
            WRITE(40,100)'> x0  (+ error, rmsDOWNHILL): '
            WRITE(40,*)X0,EX0,EEX0
            WRITE(40,100)'> sig (+ error, rmsDOWNHILL): '
            WRITE(40,*)SIGMA,ESIGMA,EESIGMA
            WRITE(40,100)'> Central wavelength [object frame]: '
            WRITE(40,*)STWV+(X0-1.)*DISP
            WRITE(40,100)'> Central wavelength [rest frame]..: '
            WRITE(40,*)(STWV+(X0-1.)*DISP)/RCVEL1
            CALL FINDNEARLINE((STWV+(X0-1.)*DISP)/RCVEL1,CLABEL0,WV0)
            WRITE(40,100)'> Nearest line is: '//
     +       CLABEL0(1:TRUELEN(CLABEL0))
            WRITE(40,*)WV0
          END IF
        END IF
C
100     FORMAT(A,$)
101     FORMAT(A)
110     FORMAT(A,I6)
150     FORMAT(79('-'))
        END
C
C******************************************************************************
C Calcula la linea mas cercana a WV1, devolviendo la etiqueta de esa linea y
C su longitud de onda. Se supone que WV1 esta dado en un sistema de referencia
C en reposo (z=0).
        SUBROUTINE FINDNEARLINE(WV1,CLABEL0,WV2)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        REAL WV1
        CHARACTER*8 CLABEL0
        REAL WV2
C
        INTEGER NLINE,NLINEST
        REAL WV(NLINMAX)
        REAL MINWV
        CHARACTER*8 CLABEL(NLINMAX)
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        CALL SELLINES(0,NLINEST,WV,CLABEL)                     !serie de Balmer
        CLABEL0=CLABEL(1)
        WV2=WV(1)
        MINWV=ABS(WV(1)-WV1)
        DO NLINE=2,NLINEST
          IF(ABS(WV(NLINE)-WV1).LT.MINWV)THEN
            MINWV=ABS(WV(NLINE)-WV1)
            WV2=WV(NLINE)
            CLABEL0=CLABEL(NLINE)
          END IF
        END DO
        CALL SELLINES(1,NLINEST,WV,CLABEL)              !lineas tipicas emision
        DO NLINE=1,NLINEST
          IF(ABS(WV(NLINE)-WV1).LT.MINWV)THEN
            MINWV=ABS(WV(NLINE)-WV1)
            WV2=WV(NLINE)
            CLABEL0=CLABEL(NLINE)
          END IF
        END DO
C
        END
