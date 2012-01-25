C------------------------------------------------------------------------------
C Version 06-December-2007                                        file: movel.f
C @ ncl & fjg
C------------------------------------------------------------------------------
Comment
C
C Program: movel.f
C Classification: arithmetic & manipulations
C Description: Computation of radial velocity and velocity dispersion using 
C              the MOVEL and OPTEMA algorithms
C
Comment
C------------------------------------------------------------------------------
C
        PROGRAM MOVEL
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER TRUELEN
        INTEGER READILIM
        REAL READF
        REAL FPOLY
        EXTERNAL CHISQG
        REAL CHISQG
        REAL RANRED
        REAL FMEDIAN1
C
        REAL PI,C
        PARAMETER (PI=3.141593)
        PARAMETER (C=299792.46)
        INTEGER NMAX_TICKS
        PARAMETER (NMAX_TICKS=100)
C Numero maximo de templates (ojo, esta relacionado con el numero maximo
C de coeficientes que pueden ajustarse; no cambiar)
        INTEGER NMAX_TEMPLATES
        PARAMETER(NMAX_TEMPLATES=20)
C Numero maximo de regiones con lineas de emision a enmascarar
        INTEGER NMAX_EML
        PARAMETER(NMAX_EML=100)
C Numero maximo de iteraciones
        INTEGER NMAX_SIMUL
        PARAMETER (NMAX_SIMUL=1000)
C
        INTEGER I,J,L,K,MJ,II,III,IDUM
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        INTEGER NS0,NC1,NC2,NC1RVEL,NC2RVEL,NC1INI,NC2INI
        INTEGER NSCAN2,NCHAN2,NCH,NP
        INTEGER NCHANEFF
        INTEGER NF1,NF2,I0,J0,NPTS,IVEL0
        INTEGER NL,IMODE,NP0
        INTEGER KMIN,KMAX,KMAXPL
        INTEGER NVAR,NEVAL
        INTEGER NITER,IITER,ITEMPL
        INTEGER NSIMUL,NSEED,ISIMUL,MAXITER
        INTEGER NTEMPLATES,NITEMP,NLINES,IVEL0PLOT
        INTEGER NSCMIN,NSCMAX,NSCENTER,IWORKING_SCAN,NSAVE
        INTEGER LDUM
        INTEGER IEML,NEML
        INTEGER IEML1(NMAX_EML),IEML2(NMAX_EML),LIN1,LIN2
        INTEGER IEML1_SGL(NMAX_EML),IEML2_SGL(NMAX_EML)
        INTEGER IEML1_STL(NMAX_EML),IEML2_STL(NMAX_EML)
        INTEGER IEML1_NCH(NMAX_EML),IEML2_NCH(NMAX_EML)
        INTEGER NTICKS,NTICKS_SMALL,MM,PP
        INTEGER NDEG,NTERMS
        INTEGER IFCOSB
        INTEGER NSKIP_MASK_TEMPLATES(NMAX_TEMPLATES)
        INTEGER LINFO,LFOUND
        REAL XF(NMAXFFT),YF(NMAXFFT)
        REAL STWV2,DISP2,WF,WI,WMIN,WI0
        REAL XMIN,XMAX,YMIN,YMAX
        REAL XMING,XMAXG,YMING,YMAXG,DXG,DYG
        REAL XMI(3),XMA(3),YMI(3),YMA(3)
        REAL GAM,VEL,SIG
        REAL BETA,RCVEL1
        REAL X(NMAXFFT)
        REAL SORIG1(NCMAX),SORIG2(NCMAX),SORIGX(NCMAX),EORIG1(NCMAX)
        REAL SPT(NCMAX),SPG(NCMAX)
        REAL STL(NMAXFFT),SGL(NMAXFFT),SML(NMAXFFT)
        REAL STL2(NMAXFFT),SGL2(NMAXFFT),SML2(NMAXFFT)
        REAL SIGC,VELO,FSCALE,ARG,PCEN,FAC
        REAL XL(NMAXFFT),YL(NMAXFFT),SIGMAY(NMAXFFT)
        REAL CO(20),CHISQR
        REAL PSEUDO_WEIGHT,PSEUDO_POWER,YRMSTOL
        REAL YLT(NMAXFFT),YLM(NMAXFFT),YLG(NMAXFFT),YR
        REAL ST(NMAXFFT),SG(NMAXFFT),SM(NMAXFFT),ST2(NMAXFFT)
        REAL COSBELL(NMAXFFT)
        REAL STR(NMAXFFT),SGR(NMAXFFT),SMR(NMAXFFT),SR(NMAXFFT)
        REAL STI(NMAXFFT),SGI(NMAXFFT),SMI(NMAXFFT),SI(NMAXFFT)
        REAL PTL(NMAXFFT),PGL(NMAXFFT),PML(NMAXFFT),PGL2(NMAXFFT)
        REAL PNOISE0,PNOISE(NMAXFFT)
        REAL XGAUS,SIGGAUS,AMPGAUS,YGAUS,EEX0,EESIGMA,EEAMP,EEY0
        REAL WIEN(NMAXFFT) 
        REAL LCUT,DLAM
        REAL XX0(3),DXX0(3),XX(3),DXX(3),XTH(3)
        REAL QR(NMAXFFT),QI(NMAXFFT),FFTA(NMAXFFT),FFTB(NMAXFFT)
        REAL B0R(NMAXFFT),B0I(NMAXFFT)
        REAL B1R(NMAXFFT),B1I(NMAXFFT)
        REAL B2R(NMAXFFT),B2I(NMAXFFT)
        REAL RESR(NMAXFFT),RESI(NMAXFFT),SGR2(NMAXFFT),SGI2(NMAXFFT)
        REAL GAMF,VELF,SIGF,VEL_PRECISION,VELF0,SIGF0
        REAL DELTA_VEL,DELTA_SIG
        REAL PI2,ERR1(NCMAX),ERR2(NCMAX),R1,R2,SXX
        REAL GAMINI,VELINI,SIGINI,EGAMF,EVELF,ESIGF
        REAL GAMSIM(NMAX_SIMUL),VELSIM(NMAX_SIMUL),SIGSIM(NMAX_SIMUL)
        REAL GAMFINAL,VELFINAL,SIGFINAL
        REAL TEMP(NCMAX,NMAX_TEMPLATES),ETEMP(NCMAX,NMAX_TEMPLATES)
        REAL WEI(NMAX_TEMPLATES),TEMPX(NCMAX,NMAX_TEMPLATES)
        REAL SIGT,VELT,SMPLOT(NMAXFFT),SGPLOT(NMAXFFT),RESID(NMAXFFT)
        REAL WLINES(NLINMAX),ZV,WLIN
        REAL GALAX(NCMAX,NSMAX),EGALAX(NCMAX,NSMAX)
        REAL GAMCEN,VELCEN,SIGCEN
        REAL GALGAM(NSMAX),GALVEL(NSMAX),GALSIG(NSMAX)
        REAL EGALGAM(NSMAX),EGALVEL(NSMAX),EGALSIG(NSMAX)
        REAL FINTGAUSS
        REAL XTICK(NMAX_TICKS),XTICK_SMALL(NMAX_TICKS)
        DOUBLE PRECISION MEAN,MEANT,MEANG,MEANM
        DOUBLE PRECISION DW,DLOGW,DV
        DOUBLE PRECISION W(NCMAX+1)
        DOUBLE PRECISION PT(NMAXFFT),PG(NMAXFFT)
        DOUBLE PRECISION PM(NMAXFFT),PG2(NMAXFFT)
        DOUBLE PRECISION FD,DBR,DBI,RESRMS,RESDUM
        CHARACTER*1 CSHOW,CLOGFILE,CITER,CSHOW2,CCONT
        CHARACTER*1 CPLOT_EML,CPAUSE
        CHARACTER*1 CERR,CREMV,CSAVE,CEML,CPLOTS
        CHARACTER*50 CDUMMY
        CHARACTER*255 INFILE1,INFILE2,ERRFILE2,ERRFILE1
        CHARACTER*255 OUTFILE,FILEHEAD,ERRFILE
        CHARACTER*8 LABLINES(NLINMAX)
        CHARACTER*5 CDUM
        CHARACTER*15 CSCANBINNING(NSMAX)
        CHARACTER*9 XDUM,XDUMF
        CHARACTER*255 CINFO
        CHARACTER*(NMAX_TEMPLATES) CMASK_TEMPLATES
        LOGICAL LNULL_ERR_TEMPLATES,LEXIST
        LOGICAL LCOLOR(MAX_ID_RED),LPAUSE
        LOGICAL FLAGIT,FLAGCA,ALLPLOT
        LOGICAL LSIMULATIONS
        LOGICAL LOPTIMAL_TEMPLATE
        LOGICAL LPROBLEM_IS_A_FRAME
        LOGICAL LREPLACE_CONTINUUM_REGIONS
        LOGICAL LPLOT_ALL,LWRITE_LOGFILE,LOGEXIST
        LOGICAL LOOP,LOUT
        LOGICAL IFCHAN(NMAXFFT)
        LOGICAL LPARIMPAR
        LOGICAL LITERATE_MORE
        LOGICAL LREPITE_TEMPLATE
C
        COMMON/BLK_X_GLOBAL/X
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
        COMMON/BLKFITG1_MOVEL/NP0
        COMMON/BLKFITG2_MOVEL/XF,YF
        COMMON/BLKFIT/NCH,KMIN,KMAX,STR,STI,SR,SI,WIEN
        COMMON/BLKFFTB/FFTA,FFTB
        COMMON/BLKTEMPL1/NTEMPLATES,NCHAN,NP,NF1,NF2,SIG,VEL,GAM
        COMMON/BLKTEMPL1B/DLOGW,DV
        COMMON/BLKTEMPL2/W,SPG,TEMPX,WEI
        COMMON/BLKTEMPL3/NEML,IEML1_SGL,IEML2_SGL
        COMMON/BLKTEMPL4/NTERMS,NDEG,CCONT
        COMMON/BLKTEMPL5/YRMSTOL
        COMMON/BLKTEMPL6/PSEUDO_WEIGHT,PSEUDO_POWER
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
        THISPROGRAM='movel'
        CALL WELCOME('06-December-2007')
C inicializamos
        NSEED=-1                           !semilla generador numero aleatorios
        LNULL_ERR_TEMPLATES=.FALSE.
        LCUT=150.0
        DO II=1,NMAX_TEMPLATES
          CMASK_TEMPLATES(II:II)='0' !0=no, 1=yes
        END DO
C
        DO J=1,NMAXFFT
          X(J)=REAL(J)
        END DO
C------------------------------------------------------------------------------
C evitamos WARNINGS de compilacion
        !esta variable solo se utiliza cuando LPROBLEM_IS_A_FRAME=.TRUE.
        NSCENTER=0
C
        IVEL0PLOT=0
C
        DLAM=0.0
C
        WI0=0.0
C
        NC1INI=0
        NC2INI=0
C
        VELF0=0.0
        SIGF0=0.0
C
        VELT=0.0
        SIGT=0.0
C
        GAMINI=0.0
        VELINI=0.0
        SIGINI=0.0
        GAMCEN=0.0
        VELCEN=0.0
        SIGCEN=0.0
C------------------------------------------------------------------------------
C abrimos salida grafica
        CALL PIDEGTER(NTERM,IDN,LCOLOR)
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGSCH(1.2)
          CALL PGSUBP(1,1)
        END DO
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C decidimos si dibujamos todos los plots
        WRITE(*,100)'Show all plots (y/n) '
        CSHOW(1:1)=READC('y','yn')
C------------------------------------------------------------------------------
C decidimos si todos los plots en cada iteracion
        IF(CSHOW.EQ.'y')THEN
          WRITE(*,100)'PAUSE between plots (y/n) '
          CPAUSE(1:1)=READC('n','yn')
          LPAUSE=(CPAUSE.EQ.'y')
          LPLOT_ALL=.TRUE.
          WRITE(*,100)'Show all plots in every iteration (y/n) '
          CSHOW2(1:1)=READC('n','yn')
        ELSE
          LPAUSE=.FALSE.
          LPLOT_ALL=.FALSE.
          CSHOW2='n'
        END IF
C------------------------------------------------------------------------------
C decidimos si dibujamos lineas de emision en el plot final
        IF(CSHOW.EQ.'y')THEN
          WRITE(*,101)'* Plot of emission lines in final fit:'
          WRITE(*,101)'(0): none'
          WRITE(*,101)'(1): intrinsic (object) emission lines'
          WRITE(*,101)'(2): sky lines'
          WRITE(*,101)'(3): intrinsic + sky lines'
          WRITE(*,100)'Option (0...3) '
          CPLOT_EML(1:1)=READC('1','0123')
        ELSE
          CPLOT_EML='0'
        END IF
C------------------------------------------------------------------------------
C podemos salvar resultados intermedios en un fichero
        WRITE(*,*)
        WRITE(*,100)'Write log file with intermediate results (y/n) '
        CLOGFILE(1:1)=READC('y','yn')
        IF(CLOGFILE.EQ.'y')THEN
          LWRITE_LOGFILE=.TRUE.
          LOGEXIST=.TRUE.
          I=0
          DO WHILE(LOGEXIST)
            I=I+1
            IF(I.EQ.100000)THEN
              STOP 'FATAL ERROR: xdum99999 reached. Delete files.'
            END IF
            WRITE(CDUM,'(I5.5)') I
            XDUM='xdum'//CDUM
            CALL RMBLANK(XDUM,XDUMF,LDUM)
            INQUIRE(FILE=XDUMF(1:LDUM),EXIST=LOGEXIST)
          END DO
          WRITE(*,101) 'Log file will be '//XDUMF
          OPEN(27,FILE=XDUMF,STATUS='NEW',FORM='FORMATTED')
        ELSE
          LWRITE_LOGFILE=.FALSE.
        END IF
C------------------------------------------------------------------------------
C podemos dedicir si vamos a estimar errores
        PI2=2.*PI
        WRITE(*,100)'Are you working with error files (y/n) '
        CERR(1:1)=READC('n','yn') 
        IF(CERR.EQ.'y')THEN
          LSIMULATIONS=.TRUE.
          WRITE(*,100)'Number of simulations for error estimation '
          NSIMUL=READILIM('40',2,NMAX_SIMUL)
        ELSE
          LSIMULATIONS=.FALSE.
          NSIMUL=0
        END IF
C------------------------------------------------------------------------------
C elegimos método para eliminar el continuo
        WRITE(*,*)
        WRITE(*,101)'* Low frequencies can be removed by either:'
        WRITE(*,101)'(1) dividing each spectrum by its fitted continuum'
        WRITE(*,100)'(2) idem + introducing galaxy continuum + dividing'
        WRITE(*,101)' by straight-line fit'
        WRITE(*,100)'Option (1/2) '
        CREMV(1:1)=READC('1','12') 
        IF(CREMV.EQ.'1')THEN
          LREPLACE_CONTINUUM_REGIONS=.TRUE.
        ELSE
          LREPLACE_CONTINUUM_REGIONS=.FALSE.
        END IF
        WRITE(*,101)'* Type of continuum fit:'
        WRITE(*,101)'(1) normal polynomial'
        WRITE(*,101)'(2) pseudo-continuum'
        WRITE(*,100)'Option (1/2) '
        CCONT(1:1)=READC('2','12')
        WRITE(*,100)'Polynomial degree for fit...'
        NDEG=READILIM('5',0,19)
        NTERMS=NDEG+1
        IF(CCONT.EQ.'2')THEN
          WRITE(*,100) 'Weight for pseudocontinuum fit......'
          PSEUDO_WEIGHT=READF('100.0')
          WRITE(*,100) 'Power  for pseudocontinuum fit........'
          PSEUDO_POWER=READF('2.0')
        END IF
        WRITE(*,100) 'YRMSTOL for DOWNHILL...............'
        YRMSTOL=READF('1.0E-5')
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C leemos el/los espectro/s de referencia (introducimos una opcion para trabajar
C con errores nulos en la template ---pero con errores de la imagen PROBLEM---)
        WRITE(*,*)
        WRITE(*,100)'REFERENCE spectrum file name'
        INFILE1=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        IF(LSIMULATIONS)THEN
          LOOP=.TRUE.
          DO WHILE(LOOP)
            WRITE(*,100)'REFERENCE error file name (NONE: errors=0)'
            CALL GUESSEF(INFILE1,ERRFILE1)
            WRITE(*,100) ' ['//ERRFILE1(1:TRUELEN(ERRFILE1))//'] ? '
            READ(*,101) ERRFILE1
            IF(ERRFILE1(1:4).EQ.'NONE')THEN
              LOOP=.FALSE.
            ELSE
              IF(TRUELEN(ERRFILE1).EQ.0)THEN
                CALL GUESSEF(INFILE1,ERRFILE1)
              END IF
              INQUIRE(FILE=ERRFILE1,EXIST=LEXIST)
              IF(LEXIST)THEN
                LOOP=.FALSE.
              ELSE
                WRITE(*,101) 'ERROR: this file does not exist.'//
     +           ' Try again.'
              END IF
            END IF
          END DO
          IF(ERRFILE1(1:4).EQ.'NONE')THEN
            LNULL_ERR_TEMPLATES=.TRUE.
          ELSE
            LNULL_ERR_TEMPLATES=.FALSE.
            ERRFILE1=INFILEX(21,ERRFILE1,NSCAN2,NCHAN2,STWV2,DISP2,11,
     +       .TRUE.)
            IF(NSCAN2.NE.NSCAN)THEN
              WRITE(*,100) 'NSCAN2, NSCAN: '
              WRITE(*,*) NSCAN2,NSCAN
              WRITE(*,101) 'FATAL ERROR: both numbers must be identical'
              CLOSE(20)
              CLOSE(21)
              STOP
            END IF
            IF(NCHAN2.NE.NCHAN)THEN
              WRITE(*,100) 'NCHAN2, NCHAN: '
              WRITE(*,*) NCHAN2,NCHAN
              WRITE(*,101) 'FATAL ERROR: both numbers must be identical'
              CLOSE(20)
              CLOSE(21)
              STOP
            END IF
            IF(STWV2.NE.STWV)THEN
              WRITE(*,100) 'STWV2, STWV: '
              WRITE(*,*) STWV2,STWV
              WRITE(*,101) 'FATAL ERROR: both numbers must be identical'
              CLOSE(20)
              CLOSE(21)
              STOP
            END IF
            IF(DISP2.NE.DISP)THEN
              WRITE(*,100) 'DISP2, DISP: '
              WRITE(*,*) DISP2,DISP
              WRITE(*,101) 'FATAL ERROR: both numbers must be identical'
              CLOSE(20)
              CLOSE(21)
              STOP
            END IF
          END IF
        END IF
        IF(NSCAN.GT.1)THEN
          WRITE(*,101)'TEMPLATE file contains more than 1 spectrum.'
          WRITE(*,100)'Scan number (0=compute optimal template)'
          NS0=READILIM('@',0,NSCAN)
        ELSE
          NS0=1
        END IF
        IF(NS0.EQ.0)THEN
          IF(NSCAN.GT.NMAX_TEMPLATES)THEN
            CLOSE(20)
            IF((LSIMULATIONS).AND.(.NOT.LNULL_ERR_TEMPLATES)) CLOSE(21)
            IF(LWRITE_LOGFILE) CLOSE(27)
            WRITE(*,101) 'FATAL ERROR: NSCAN.GT.NMAX_TEMPLATES'
            WRITE(*,100) 'NSCAN, NMAX_TEMPLATES: '
            WRITE(*,*) NSCAN,NMAX_TEMPLATES
            STOP
          END IF
          NTEMPLATES=0
          DO WHILE(NTEMPLATES.LT.2)
            DO I=1,NSCAN
              CMASK_TEMPLATES(I:I)='1' !0=no, 1=yes
            END DO
            WRITE(*,100)'Mask for templates (0=no, 1=yes) '
            CMASK_TEMPLATES(1:NSCAN)=READC(CMASK_TEMPLATES(1:NSCAN),
     +       '01')
            NTEMPLATES=0
            DO I=1,NSCAN
              IF(CMASK_TEMPLATES(I:I).EQ.'1') NTEMPLATES=NTEMPLATES+1
            END DO
            IF(NTEMPLATES.LT.2)THEN
              WRITE(*,101) 'ERROR: number of templates must be'//
     +         ' at least 2'
            END IF
          END DO
          WRITE(*,100) 'Number of templates to be used: '
          WRITE(*,*) NTEMPLATES
          LOPTIMAL_TEMPLATE=.TRUE.
          !calculamos cuantos scans tendremos que saltarnos para leer la
          !template correspondiente
          I=0
          II=1
          NSKIP_MASK_TEMPLATES(II)=0
          DO WHILE(I.LE.NSCAN)
            I=I+1
            IF(CMASK_TEMPLATES(I:I).EQ.'0')THEN
              NSKIP_MASK_TEMPLATES(II)=NSKIP_MASK_TEMPLATES(II)+1
            ELSE
              II=II+1
              NSKIP_MASK_TEMPLATES(II)=0
            END IF
          END DO
        ELSE
          NTEMPLATES=1
          LOPTIMAL_TEMPLATE=.FALSE.
        END IF
!       NITEMP=0 !Creo que sobra (NCL)
C
        DO II=1,NTEMPLATES
          IF(.NOT.LOPTIMAL_TEMPLATE)THEN
            DO I=1,NS0
              READ(20) (SORIG1(J),J=1,NCHAN)
            END DO
            IF(LSIMULATIONS)THEN
              IF(LNULL_ERR_TEMPLATES)THEN
                DO J=1,NCHAN
                  ERR1(J)=0.0
                END DO
              ELSE
                DO I=1,NS0
                  READ(21) (ERR1(J),J=1,NCHAN)
                END DO
              END IF
            END IF
          ELSE
            IF(NSKIP_MASK_TEMPLATES(II).GT.0)THEN
              DO III=1,NSKIP_MASK_TEMPLATES(II)
                READ(20) (SORIG1(J),J=1,NCHAN)
              END DO
            END IF
            READ(20) (SORIG1(J),J=1,NCHAN)
            IF(LSIMULATIONS)THEN
              IF(LNULL_ERR_TEMPLATES)THEN
                DO J=1,NCHAN
                  ERR1(J)=0.0
                END DO
              ELSE
                IF(NSKIP_MASK_TEMPLATES(II).GT.0)THEN
                  DO III=1,NSKIP_MASK_TEMPLATES(II)
                    READ(21) (ERR1(J),J=1,NCHAN)
                  END DO
                END IF
                READ(21) (ERR1(J),J=1,NCHAN)
              END IF
            END IF
          END IF
          !normalizamos datos y errores
          MEAN=0.D0
          DO J=1,NCHAN
            MEAN=MEAN+DBLE(SORIG1(J))
          END DO
          MEAN=MEAN/DBLE(NCHAN)
          DO J=1,NCHAN
            SORIG1(J)=SORIG1(J)/REAL(MEAN)
          END DO
          IF(LSIMULATIONS)THEN
            DO J=1,NCHAN
              ERR1(J)=ERR1(J)/REAL(MEAN)
            END DO
          END IF
          IF(LOPTIMAL_TEMPLATE)THEN
            DO J=1,NCHAN
              TEMP(J,II)=SORIG1(J)
              TEMPX(J,II)=SORIG1(J)
            END DO
            IF(LSIMULATIONS)THEN
              DO J=1,NCHAN
                ETEMP(J,II)=ERR1(J)
              END DO
            END IF
          END IF
        END DO
        CLOSE(20)
        IF((LSIMULATIONS).AND.(.NOT.LNULL_ERR_TEMPLATES)) CLOSE(21)
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C leemos el/los espectros problema
        WRITE(*,100)'PROBLEM spectrum file name..'
        INFILE2=INFILEX(20,'@',NSCAN2,NCHAN2,STWV2,DISP2,1,.FALSE.)
        IF((NCHAN.NE.NCHAN2).OR.(STWV.NE.STWV2).OR.(DISP.NE.DISP2))THEN
          WRITE(*,101)'NCHAN (template, problem): '
          WRITE(*,*) NCHAN,NCHAN2
          WRITE(*,100)'SWTV  (template, problem): '
          WRITE(*,*) STWV,STWV2
          WRITE(*,100)'DISP  (template, problem): '
          WRITE(*,*) DISP,DISP2
          WRITE(*,101)'FATAL ERROR: header information in PROBLEM '//
     +     'spectrum is different.'
          CLOSE(20)
          STOP
        END IF
        IF(LSIMULATIONS)THEN
          WRITE(*,100)'PROBLEM error file name.....'
          CALL GUESSEF(INFILE2,ERRFILE2)
          ERRFILE2=INFILEX(21,ERRFILE2,NSCAN2,NCHAN2,STWV2,DISP2,
     c     21,.TRUE.)
        END IF
        IF(NSCAN2.GT.1)THEN
          WRITE(*,101)'This file contains more than 1 spectrum.'
          WRITE(*,100)'Indicate the scan number to be employed '
          WRITE(*,101)'(0 to measure all/some spectra)'
          WRITE(*,100)'Scan'
          NS0=READILIM('@',0,NSCAN2)
          IF(NS0.EQ.0)THEN
            WRITE(*,101)'NOTE: the central spectrum will be'//
     +       ' analyzed first'
          END IF
        ELSE
          NS0=1
        END IF
        IF(NS0.EQ.0)THEN
          LPROBLEM_IS_A_FRAME=.TRUE.
          WRITE(*,100) 'Central spectrum'
          WRITE(CDUMMY,*) NSCAN2/2+1
          NSCENTER=READILIM(CDUMMY,1,NSCAN2)
          WRITE(CDUMMY,'(A,I10)') '1,',NSCAN2
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          LOOP=.TRUE.
          DO WHILE(LOOP)
            WRITE(*,100) 'Minimum and maximum spectra to use '
            CALL READ2I(CDUMMY(1:L),NSCMIN,NSCMAX)
            IF((NSCMIN.LT.1).OR.(NSCMAX.GT.NSCAN2).OR.
     +       (NSCMIN.GT.NSCMAX))THEN
              WRITE(*,101)'ERROR: scan range out of range. Try again.'
            ELSE
              LOOP=.FALSE.
            END IF
          END DO
          IF(LPLOT_ALL)THEN
            WRITE(*,100)'Plots for all spectra (y/n) '
            CPLOTS(1:1)=READC('n','yn') 
            IF(CPLOTS.EQ.'n') LPLOT_ALL=.FALSE.
          END IF
          DO I=1,NSCAN2
            READ(20) (GALAX(J,I),J=1,NCHAN)
          END DO
          IF(LSIMULATIONS)THEN
            DO I=1,NSCAN2
              READ(21) (EGALAX(J,I),J=1,NCHAN)
            END DO
          END IF
        ELSE
          LPROBLEM_IS_A_FRAME=.FALSE.
          DO I=1,NS0
            READ(20) (SORIG2(J),J=1,NCHAN)
          END DO
          IF(LSIMULATIONS)THEN
            DO I=1,NS0
              READ(21) (ERR2(J),J=1,NCHAN)
            END DO
          END IF
        END IF
        CLOSE(20)
        IF(LSIMULATIONS) CLOSE(21)
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C Podemos enmmascar regiones problemáticas del espectro a medir sustituyendolas
C por el espectro de referencia en esas mismas regiones
        WRITE(*,*)
        WRITE(*,101)'Spectral regions of the PROBLEM '//
     +   'spectrum can be replaced by the corresponding'
        WRITE(*,101)'MODEL regions. This is useful to'//
     +   ' avoid badly removed sky lines and intrinsic'
        WRITE(*,101)'emission lines of the PROBLEM'//
     +   ' spectrum. These regions will also be ignored'
        WRITE(*,101)'in the continuum fitting of'//
     +   ' the PROBLEM spectrum for the removal of low spatial'
        WRITE(*,101)'frequencies.'
        WRITE(*,100)'Are you using this replacement option (y/n) '
        CEML(1:1)=READC('n','yn')
        IF(CEML.EQ.'y')THEN
          NEML=1
          LOOP=.TRUE.
          DO WHILE(LOOP)
            WRITE(*,100) 'Channel regions to be interpolated '//
     +       '(0,0 to exit) '
            CALL READ2I('0,0',IEML1(NEML),IEML2(NEML))
            IF((IEML1(NEML).EQ.0).AND.(IEML2(NEML).EQ.0))THEN
              LOOP=.FALSE.
            ELSE
              IF((IEML1(NEML).LT.1).OR.(IEML2(NEML).GT.NCHAN).OR.
     +           (IEML1(NEML).GT.IEML2(NEML)))THEN
                WRITE(*,101) 'ERROR: invalid entry. Try again.'
              ELSE
                NEML=NEML+1
                IF(NEML.EQ.NMAX_EML)THEN
                  IF(LWRITE_LOGFILE) CLOSE(27)
                  STOP 'FATAL ERROR: NEML.EQ.NMAX_EML. Redim NMAX_EML.'
                END IF
              END IF
            END IF
          END DO
          NEML=NEML-1
        ELSE
          NEML=0
        END IF
C------------------------------------------------------------------------------
C pedimos numero de iteraciones
        MAXITER=31 !maximum number of iterations
        WRITE(*,*)
        WRITE(*,101)'* Iterations (to refine parameters and '//
     +   'optimal template):'
        WRITE(*,101)'-1=FREE -> user must decide to iterate'
        WRITE(*,101)' 0=AUTO -> end iteration at fixed velocity '//
     +   'precision'
        WRITE(*,101)' n=FIXED-> n=1,2,... fixed number of iterations'
        WRITE(*,100)'Option'
        NITER=READILIM('0',-1,1000)
        IF(NITER.EQ.0)THEN
          WRITE(*,100)'Precision in both velocity and sigma (km/s) '
          VEL_PRECISION=READF('0.1')
          WRITE(*,100)'Maximum number of iterations in any case'
          WRITE(CDUMMY,*) MAXITER
          MAXITER=READILIM(CDUMMY,1,1000)
        ELSE
          IF(NITER.GT.0) MAXITER=NITER
          VEL_PRECISION=0.
        END IF
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C IWORKING_SCAN=0 para el scan central o para el unico scan de la imagen 
C PROBLEM. En el resto de los casos, si la imagen PROBLEM tiene más de un scan,
C IWORKING_SCAN es el número de scan que vamos a medir.
        IWORKING_SCAN=0 !empezamos siempre por el scan central
C------------------------------------------------------------------------------
C Volvemos a este punto cada vez que cambiemos de scan
1000    CONTINUE
        !si tenemos varias templates para calcular una template óptima,
        !volvemos a indicar que habrá que calcular dicha template óptima
        !en el nuevo scan
        IF(NTEMPLATES.GT.1) LOPTIMAL_TEMPLATE=.TRUE.
        NITEMP=0
        ISIMUL=0
        IF(LPROBLEM_IS_A_FRAME)THEN
          IF(IWORKING_SCAN.EQ.0)THEN
            IDUM=NSCENTER
          ELSE
            IDUM=IWORKING_SCAN
          END IF
          IF(IWORKING_SCAN.EQ.NSCENTER)THEN
            !si es el scan central, ya lo hemos hecho y seguimos al siguiente
            IWORKING_SCAN=IWORKING_SCAN+1
            GOTO 1000
          END IF
          DO J=1,NCHAN
            SORIG2(J)=GALAX(J,IDUM)
          END DO
          IF(LSIMULATIONS)THEN
            DO J=1,NCHAN
              ERR2(J)=EGALAX(J,IDUM)
            END DO
          END IF
        END IF
        !normalizamos datos y errores (evitando regiones a enmascarar)
        DO J=1,NCHAN
          IFCHAN(J)=.TRUE.
        END DO
        IF(NEML.GT.0)THEN
          DO IEML=1,NEML
            DO J=IEML1(IEML),IEML2(IEML)
              IFCHAN(J)=.FALSE.
            END DO
          END DO
        END IF
        NCHANEFF=0
        MEAN=0.D0
        DO J=1,NCHAN
          IF(IFCHAN(J))THEN
            MEAN=MEAN+DBLE(SORIG2(J))
            NCHANEFF=NCHANEFF+1
          END IF
        END DO
        MEAN=MEAN/DBLE(NCHANEFF)
        DO J=1,NCHAN
          SORIG2(J)=SORIG2(J)/REAL(MEAN)
        END DO
        IF(LSIMULATIONS)THEN
          DO J=1,NCHAN
            ERR2(J)=ERR2(J)/REAL(MEAN)
          END DO
        END IF
C------------------------------------------------------------------------------
C Por aquí podemos pasar dos veces: la primera vez mostrando todas los 
C espectros de TEMPLATE disponibles y la segunda vez con la TEMPLATE óptima. 
C Si sólo tenemos un espectro TEMPLATE, sólo pasamos por aquí la primera vez.
750     CONTINUE
!       NC1=1     !creo que no hace falta (NCL)
!       NC2=NCHAN !creo que no hace falta (NCL)
C------------------------------------------------------------------------------
C plot TEMPLATE and PROBLEM spectra
        IF((CSHOW.EQ.'y').AND.(NITEMP.LE.1).AND.(ISIMUL.EQ.0)
     +   .AND.((IWORKING_SCAN.EQ.0).OR.(.NOT.LOPTIMAL_TEMPLATE)))THEN
          !para el PROBLEM spectrum, calculamos los límites ignorando
          !regiones a interpolar
          DO J=1,NCHAN
            IFCHAN(J)=.TRUE.
          END DO
          IF(NEML.GT.0)THEN
            DO IEML=1,NEML
              DO J=IEML1(IEML),IEML2(IEML)
                IFCHAN(J)=.FALSE.
              END DO
            END DO
          END IF
          CALL FINDMML(NCHAN,1,NCHAN,X,XMING,XMAXG)
          CALL FINDMMLMASK(NCHAN,1,NCHAN,SORIG2,IFCHAN,YMING,YMAXG)
          DXG=0.05*(XMAXG-XMING)
          XMING=XMING-DXG
          XMAXG=XMAXG+DXG
          DYG=0.05*(YMAXG-YMING)
          YMING=YMING-DYG
          YMAXG=YMAXG+DYG
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGPAGE
            CALL AUTOPLOT(NCHAN,X,SORIG2,1,NCHAN,             !PROBLEM spectrum
     +       'channel (RAW DATA)','normalized flux',' ',
     +       .FALSE.,XMING,XMAXG,YMING,YMAXG,0.05,
     +       0,.TRUE.,'BCNTS','CMTS',
     +       101,2,
     +       0.0,1.0,0.05,1.0)
            IF(LOPTIMAL_TEMPLATE)THEN
              DO II=1,NTEMPLATES
                DO J=1,NCHAN
                  SORIG1(J)=TEMP(J,II)
                END DO
                CALL AUTOPLOT(NCHAN,X,SORIG1,1,NCHAN,         !TEMPLATE spectra
     +           ' ',' ',' ',
     +           .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +            0,.FALSE.,'BCNTS','B',
     +            101,3,
     +            0.0,1.0,0.05,1.0)
              END DO
            ELSE
              CALL AUTOPLOT(NCHAN,X,SORIG1,1,NCHAN,          !TEMPLATE spectrum
     +         ' ',' ',' ',
     +         .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +          0,.FALSE.,'BCNTS','BNTS',
     +          101,3,
     +          0.0,1.0,0.05,1.0)
            END IF
            CALL PGSCI(2)
            CALL PGMTEXT('T',2.,0.,0.,'PROBLEM: '//INFILE2)
            CALL PGSCI(3)
            CALL PGMTEXT('T',1.,0.,0.,'TEMPLATE: '//INFILE1)
            CALL PGSCI(1)
            IF(NEML.GT.0)THEN
              CALL PGSLS(4)
              DO IEML=1,NEML
                CALL PGMOVE(REAL(IEML1(IEML)),YMIN)
                CALL PGDRAW(REAL(IEML1(IEML)),YMAX)
                CALL PGMOVE(REAL(IEML2(IEML)),YMIN)
                CALL PGDRAW(REAL(IEML2(IEML)),YMAX)
              END DO
              CALL PGSLS(1)
            END IF
            CALL PGIDEN_RED
          END DO
          IF(LPAUSE)THEN
            WRITE(*,100) 'Press <CR> to continue...'
            READ(*,*)
          END IF
        END IF
C------------------------------------------------------------------------------
C Simulaciones de Monte-Carlo; se realizan sólo en el caso de que así se
C haya solicitado, cuando el número de simulación no es cero (ISIMUL.GT.0) y
C cuando, además, NITEMP=0.
333     CONTINUE
        IF((LSIMULATIONS).AND.(ISIMUL.GT.0).AND.(NITEMP.EQ.0))THEN
          !simulaciones para TEMPLATE
          IF(NTEMPLATES.EQ.1)THEN
            DO J=1,NCHAN
              R1=RANRED(NSEED)
              R2=RANRED(NSEED)
              SXX=1.414213562*ERR1(J)*SQRT(-1.*LOG(1.-R1))*COS(PI2*R2)
              SPT(J)=SORIG1(J)+SXX
            END DO
          ELSE
            DO II=1,NTEMPLATES
              DO J=1,NCHAN
                R1=RANRED(NSEED)
                R2=RANRED(NSEED)
                SXX=1.414213562*ETEMP(J,II)*SQRT(-1.*LOG(1.-R1))*
     +           COS(PI2*R2)
                TEMPX(J,II)=TEMP(J,II)+SXX
              END DO
            END DO
          END IF
          !simulaciones para PROBLEM
          DO J=1,NCHAN
            R1=RANRED(NSEED)
            R2=RANRED(NSEED)
            SXX=1.414213562*ERR2(J)*SQRT(-1.*LOG(1.-R1))*COS(PI2*R2)
            SPG(J)=SORIG2(J)+SXX
          END DO
        ELSE !.....................no toca (o todavía no toca) hacer simulación
          IF(.NOT.LOPTIMAL_TEMPLATE)THEN
            DO J=1,NCHAN
              SPT(J)=SORIG1(J)
            END DO
          END IF
          IF(NITEMP.EQ.0)THEN
            DO J=1,NCHAN
              SPG(J)=SORIG2(J)
            END DO
          END IF
        END IF
C------------------------------------------------------------------------------
C volvemos a este punto cada vez que tengamos que recalcular la template óptima
C (en ese caso, se reinicializan las iteraciones IITER)
444     IITER=0
        FLAGIT=.TRUE.
        ALLPLOT=.FALSE.
        IF((CSHOW.EQ.'y').AND.(ISIMUL.EQ.0)) ALLPLOT=.TRUE.
C------------------------------------------------------------------------------
C volvemos a este punto cada vez que volvemos a iterar
222     CONTINUE
        IF(CSHOW2.NE.'y')THEN
          IF((ISIMUL.GT.0).OR.(IITER.GT.0)) ALLPLOT=.FALSE.
          IF((NTEMPLATES.GT.1).AND.(NITEMP.GT.1)) ALLPLOT=.FALSE.
        END IF
C------------------------------------------------------------------------------
C introduccion de los primeros parametros
        IF((FLAGIT).AND.(NITEMP.LE.1))THEN
          IF((ISIMUL.EQ.0).AND.(NITEMP.EQ.0))THEN
            IF(IWORKING_SCAN.EQ.0)THEN
              !pedimos los "guess parameters"
              WRITE(*,*)
              WRITE(*,101)'* Introduce first guess parameters:'
              WRITE(*,100)'--> gamma (relative line-strength).........'
              GAM=READF('1.0')
              GAMINI=GAM
              WRITE(*,100)'--> radial velocity (in km/s)..............'
              VEL=READF('0.0')
              VELINI=VEL
              WRITE(*,100)'--> sigma (velocity dispersion, in km/s) '
              SIG=READF('300.0')
              SIGINI=SIG
              WRITE(*,*)
            ELSE
              !tomamos como valores iniciales los ajustados en el scan central
              GAM=GAMCEN
              GAMINI=GAM
              VEL=VELCEN
              VELINI=VEL
              SIG=SIGCEN
              SIGINI=SIG
            END IF
          ELSE
            !usamos los "guess parameters"
            GAM=GAMINI
            VEL=VELINI
            SIG=SIGINI
          END IF
        ELSE
          !usamos los valores de la última iteración
          GAM=GAMF
          VEL=VELF
          SIG=SIGF
        END IF
C------------------------------------------------------------------------------
C evitamos que se elija una zona para la que no haya solape entre el espectro 
C de la template y el espectro del objeto problema
        BETA=VEL/C
        RCVEL1=(1.0+BETA)/SQRT(1.0-BETA*BETA)            !correcion relativista
        IF(VEL.GE.0.0)THEN
          WF=STWV+DISP*REAL(NCHAN-1)
          WI=WF/RCVEL1
          NC2RVEL=NINT((WI-STWV)/DISP)+1
          IF(NC2RVEL.LT.1)THEN
            WRITE(*,101) 'FATAL ERROR: with this radial velocity, '//
     +       'the TEMPLATE and PROBLEM spectra'
            WRITE(*,101) '             do not overlap'
            IF(LWRITE_LOGFILE) CLOSE(27)
            STOP
          END IF
          NC1RVEL=1
        ELSE
          WF=STWV
          WI=WF/RCVEL1
          NC1RVEL=NINT((WI-STWV)/DISP)+1
          IF(NC1RVEL.GT.NCHAN)THEN
            WRITE(*,101) 'FATAL ERROR: with this radial velocity, '//
     +       'the TEMPLATE and PROBLEM spectra'
            WRITE(*,101) '             do not overlap'
            IF(LWRITE_LOGFILE) CLOSE(27)
            STOP
          END IF
          NC2RVEL=NCHAN
        END IF
        FLAGCA=.FALSE.
C definimos region a utilizar de los espectros
        IF(FLAGIT)THEN
          IF((ISIMUL.EQ.0).AND.(NITEMP.EQ.0).AND.
     +     (IWORKING_SCAN.EQ.0))THEN
            LOOP=.TRUE.
            DO WHILE(LOOP)
              WRITE(CDUMMY,'(I10,A1,I10)')NC1RVEL,',',NC2RVEL
              CALL RMBLANK(CDUMMY,CDUMMY,L)
              WRITE(*,101)'* Introduce 1st & last channel to be'//
     +         ' employed'
              WRITE(*,100)'NC1,NC2 (in rest-frame = template-frame) '
              CALL READ2I(CDUMMY(1:L),NC1,NC2)
              IF(NC1.LT.NC1RVEL)THEN
                WRITE(*,100) 'NC1,NC1RVEL: '
                WRITE(*,*) NC1,NC1RVEL
                WRITE(*,101)'ERROR: NC1.LT.NC1RVEL. Try again.'
              ELSEIF(NC2.GT.NC2RVEL)THEN
                WRITE(*,100) 'NC2,NC2RVEL: '
                WRITE(*,*) NC2,NC2RVEL
                WRITE(*,101)'ERROR: NC2.GT.NC2RVEL. Try again.'
              ELSEIF(NC1.GT.NC2)THEN
                WRITE(*,100) 'NC1,NC2: '
                WRITE(*,*) NC1,NC2
                WRITE(*,101)'ERROR: NC1.GT.NC2. Try again.'
              ELSE
                NC1INI=NC1
                NC2INI=NC2
                LOOP=.FALSE.
              END IF
            END DO
            WRITE(*,*)
          ELSE
            !re-establecemos los valores introducidos la primera vez
            NC1=NC1INI
            NC2=NC2INI
          END IF
          FLAGCA=.TRUE.
        END IF
C Comprobamos si, con la nueva velocidad radial, estamos trabajando con un 
C intervalo correcto de canales
        IF(NC1.LT.NC1RVEL)THEN
          WRITE(*,100)'NC1,NC1RVEL: '
          WRITE(*,*) NC1,NC1RVEL
          WRITE(*,101)'WARNING: NC1.LT.NC1RVEL. Forcing NC1=NC1RVEL.'
          FLAGCA=.TRUE.
          NC1=NC1RVEL
        ELSEIF(NC2.GT.NC2RVEL)THEN
          WRITE(*,100)'NC2,NC2RVEL: '
          WRITE(*,*) NC2,NC2RVEL
          WRITE(*,101)'WARNING: NC2.GT.NC2RVEL. Forcing NC2=NC2RVEL.'
          FLAGCA=.TRUE.
          NC2=NC2RVEL
        END IF
C------------------------------------------------------------------------------
C Si ha habido alguna modificación de los límites NC1,NC2, recalculamos los 
C parámetros afectados
        IF(FLAGCA)THEN
C..............................................................................
C Pasamos el espectro a escala logaritmica
          NCHAN2=NC2-NC1+1         !numero de canales utiles en escala original
          K=0
          LOOP=.TRUE.
          DO WHILE(LOOP)
            K=K+1
            NCH=2**K     !numero de canales utiles con los que vamos a trabajar
            IF(NCHAN2.LT.NCH) LOOP=.FALSE.
          END DO
          DLAM=DISP*REAL(NC2-NC1+1)  !incremento en ldo (se usara en el filtro)
          WMIN=STWV-DISP/2.            !l.d.o. borde izquierdo del primer pixel
          DW=DBLE(DISP)
          DO MJ=1,NCHAN+1        !modificado en 1 respecto a la rutina original
            W(MJ)=DBLE(WMIN)+DBLE(MJ-1)*DW !l.d.o. borde izquierdo del pixel MJ
          END DO
C..............................................................................
C Wavelength bin (a partir de la zona util)
          DLOGW=(DLOG(W(NC2+1)/W(NC1)))/DBLE(NCH)
          DV=DBLE(C)*DLOGW                               !velocity bin (km/sec)
          IF((ISIMUL.EQ.0).AND.(NITEMP.EQ.0).AND.
     +     (IWORKING_SCAN.EQ.0))THEN
            WRITE(*,100)'>>> NCHAN (original scale)................: '
            WRITE(*,*) NCHAN
            WRITE(*,100)'>>> Useful channel region (original scale): '
            WRITE(*,*) NCHAN2
            WRITE(*,100)'>>> NCH (nearest power of 2)..............: '
            WRITE(*,*) NCH
            WRITE(*,100)'>>> Wavelength bin (log. scale)...........: '
            WRITE(*,*)DLOGW
            WRITE(*,100)'>>> Velocity bin (km/sec, log. scale).....: '
            WRITE(*,*)DV
          END IF
          DO MJ=1,NCHAN+1
            W(MJ)=DLOG(W(MJ))   !log. neperianos l.d.o. borde izquierdo pixel J
          END DO
C..............................................................................
C Como el wavelength bin se ha calculado para repartir el intervalo espectral 
C con el que vamos a trabajar [NC1,NC2] en todos los canales NCH, y como 
C NC1 y NC2 pueden ser diferentes de 1 y NCHAN, respectivamente, con la escala 
C así calculada se obtiene que al volver de nuevo a la escala original, el 
C espectro que se extendía entre 1 y NCHAN tiene entonces un tamaño que es NP, 
C que puede ser mayor que NCH.
          NP=INT((W(NCHAN+1)-W(1))/DLOGW)
          IF(NP.GT.NMAXFFT)THEN
            WRITE(*,100) 'NP, NMAXFFT: '
            WRITE(*,*) NP,NMAXFFT
            WRITE(*,101) 'FATAL ERROR: NP.GT.NMAXFFT. Redim NMAXFFT.'
            IF(LWRITE_LOGFILE) CLOSE(27)
            STOP
          END IF
C..............................................................................
C compute the NCH pixels in which we are going to have the useful region
          WI0=STWV+DISP*REAL(NC1-1)              !l.d.o. del canal NC1 original
          NF1=INT((LOG(WI0)-REAL(W(1)))/REAL(DLOGW))+1
          NF2=NF1+NCH-1
          IF(ISIMUL.EQ.0.AND.NITEMP.EQ.0.AND.IWORKING_SCAN.EQ.0)THEN
            WRITE(*,100)'>>> NP (full spectrum in log. scale)......: '
            WRITE(*,*) NP
            WRITE(*,100)'>>> NF1,NF2 (=NC1,NC2 in last scale)......: '
            WRITE(*,*) NF1,NF2
          END IF
C..............................................................................
C regiones a enmascarar en el sistema de la galaxia (se calcula aqui porque
C puede hacer falta para la determinacion de la template optima)
          IF(NEML.GT.0)THEN
            DO IEML=1,NEML
              WLIN=STWV+REAL(IEML1(IEML)-1)*DISP
              IEML1_SGL(IEML)=INT((LOG(WLIN)-REAL(W(1)))/REAL(DLOGW))+1
              WLIN=STWV+REAL(IEML2(IEML)-1)*DISP
              IEML2_SGL(IEML)=INT((LOG(WLIN)-REAL(W(1)))/REAL(DLOGW))+1
            END DO
          END IF
C..............................................................................
C Dibujamos la region seleccionada en el plot con los espectros iniciales
          IF((CSHOW.EQ.'y').AND.(NITEMP.LE.1).AND.(ISIMUL.EQ.0)
     +     .AND.(.NOT.LOPTIMAL_TEMPLATE))THEN
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              CALL PGSLS(2)
              CALL PGSCI(3)
              CALL PGMOVE(REAL(NC1),YMIN)
              CALL PGDRAW(REAL(NC1),YMAX)
              CALL PGMOVE(REAL(NC2),YMIN)
              CALL PGDRAW(REAL(NC2),YMAX)
              CALL PGSCI(1)
              CALL PGSLS(1)
            END DO
          END IF
C..............................................................................
C Here we compute the ideal template by finding the best linear combination
C of the input templates. 
C Notar que sólo cuando se recalcula la template óptima (CALL FITTEMPL) cuando 
C trabajamos con el scan central o único (IWORKING_SCAN.EQ.0), y cuando no 
C estamos haciendo simulaciones (ISIMUL.EQ.0). En los demás scans (y en las 
C simulaciones correspondientes) se utilizan los mismos pesos WEI() calculados
C la primera vez.
          IF(LOPTIMAL_TEMPLATE)THEN
            IF((ISIMUL.EQ.0).AND.(IWORKING_SCAN.EQ.0))THEN
              CALL FITTEMPL
            END IF
            DO J=1,NCHAN
              SORIG1(J)=0.
              DO I=1,NTEMPLATES
                SORIG1(J)=SORIG1(J)+TEMPX(J,I)*WEI(I)
              END DO
            END DO
            LOPTIMAL_TEMPLATE=.FALSE.
            NITEMP=NITEMP+1
            ITEMPL=(NITEMP+1)/2
            GOTO 750
          ELSE
            !conversion to log scale
            CALL LOGSCALE(NCHAN,NP,DLOGW,W,SPT,STL)
            CALL LOGSCALE(NCHAN,NP,DLOGW,W,SPG,SGL)
          END IF
C..............................................................................
        END IF
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C create the model spectra
        SIGC=SIG/REAL(DV)
!       VELO=VEL/REAL(DV)                                                  !mal
        BETA=VEL/C
        RCVEL1=(1.0+BETA)/SQRT(1.0-BETA*BETA)            !correcion relativista
        VELO=ALOG(RCVEL1)/DLOGW
        FSCALE=1./SQRT(2*PI)/SIGC
        ARG=-1./2./SIGC/SIGC
        I0=INT(10.*SIGC+0.5) !integramos en +/- 10 veces sigma
        !antes de hacer la integral, miramos si I0=0; en ese caso no es
        !necesario y evitamos errores de redondeo al llamar a FINTGAUSS
        IF(I0.EQ.0)THEN
          DO I=1,NP
            PCEN=REAL(I)-VELO
            J0=INT(PCEN+0.5)
            IF((J0.GE.1).AND.(J0.LE.NP))THEN
              SML(I)=STL(J0)
            ELSE
              SML(I)=0.0
            END IF
          END DO
        ELSE !integramos numéricamente
          DO I=1,NP
            SML(I)=0.
            PCEN=REAL(I)-VELO
            DO J=-I0,I0
              J0=INT(PCEN+0.5)+J
              IF((J0.GE.1).AND.(J0.LE.NP))THEN
                FAC=FSCALE*
     +           FINTGAUSS(REAL(J0)-0.5,REAL(J0)+0.5,20,PCEN,ARG)
                SML(I)=SML(I)+STL(J0)*FAC
              END IF
            END DO
          END DO
        END IF
C------------------------------------------------------------------------------
C calculamos las regiones a interpolar (IEML?_SGL: sistema de la
C galaxia, IEML?_STL: sistema de la template)
        IVEL0=INT(VELO+0.5)
        IF(ISIMUL.EQ.0) IVEL0PLOT=IVEL0
        IF(NEML.GT.0)THEN
          DO IEML=1,NEML
            WLIN=STWV+REAL(IEML1(IEML)-1)*DISP
            IEML1_SGL(IEML)=INT((LOG(WLIN)-REAL(W(1)))/REAL(DLOGW))+1
            IEML1_STL(IEML)=IEML1_SGL(IEML)-IVEL0
            IEML1_NCH(IEML)=IEML1_STL(IEML)-NF1+1
            WLIN=STWV+REAL(IEML2(IEML)-1)*DISP
            IEML2_SGL(IEML)=INT((LOG(WLIN)-REAL(W(1)))/REAL(DLOGW))+1
            IEML2_STL(IEML)=IEML2_SGL(IEML)-IVEL0
            IEML2_NCH(IEML)=IEML2_STL(IEML)-NF1+1
            IF((FLAGIT).AND.(NITEMP.LE.1).AND.(ISIMUL.EQ.0))THEN
              WRITE(*,*)
              WRITE(*,'(A,I2.2,A)') '> PROBLEM region #',IEML,
     +         ' to be replaced by TEMPLATE:'
              WRITE(*,'(A,2(1X,I6))') 'raw data....................: ',
     +         IEML1(IEML),IEML2(IEML)
              WRITE(*,'(A,2(1X,I6))') 'logscale [1:np].............: ',
     +         IEML1_SGL(IEML),IEML2_SGL(IEML)
              WRITE(*,'(A,2(1X,I6))') 'logscale [1:np], unreshifted: ',
     +         IEML1_STL(IEML),IEML2_STL(IEML)
              WRITE(*,'(A,2(1X,I6))') 'logscale [1:nch]............: ',
     +         IEML1_NCH(IEML),IEML2_NCH(IEML)
            END IF
          END DO
        END IF
C------------------------------------------------------------------------------
C plot of the spectra in logarithmic scale
        IF(ALLPLOT)THEN
          !Como hemos hecho un rebinning de los datos, los espectros, 
          !incialmente normalizados a uno, dejan de estar normalizados
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGSUBP(2,2)
            CALL PGPAGE
            !para el PROBLEM spectrum, calculamos los límites ignorando
            !regiones a interpolar
            DO J=1,NP
              IFCHAN(J)=.TRUE.
            END DO
            IF(NEML.GT.0)THEN
              DO IEML=1,NEML
                DO J=IEML1_SGL(IEML),IEML2_SGL(IEML)
                  IFCHAN(J)=.FALSE.
                END DO
              END DO
            END IF
            CALL FINDMML(NP,1,NP,X,XMING,XMAXG)
            CALL FINDMMLMASK(NP,1,NP,SGL,IFCHAN,YMING,YMAXG)
            DXG=0.05*(XMAXG-XMING)
            XMING=XMING-DXG
            XMAXG=XMAXG+DXG
            DYG=0.05*(YMAXG-YMING)
            YMING=YMING-DYG
            YMAXG=YMAXG+DYG
            CALL AUTOPLOT(NP,X,SGL,1,NP, !..............................PROBLEM
     +       'channel (LOG SCALE [1:NP])','rebinned flux',' ',
     +       .FALSE.,XMING,XMAXG,YMING,YMAXG,0.05,
     +       0,.TRUE.,'BCNTS','CMTS',
     +       101,2,
     +       0.,1.,0.05,1.0)
            CALL AUTOPLOT(NP,X,STL,1,NP, !.............................TEMPLATE
     +       ' ',' ' ,' ',
     +       .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +       0,.FALSE.,'BCNTS','BNTS',
     +       101,3,
     +       0.,1.,0.05,1.0)
            CALL AUTOPLOT(NP,X,SML,1,NP, !................................MODEL
     +       ' ',' ',' ',
     +       .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +       0,.FALSE.,'BCNTS','BC',
     +       101,4,
     +       0.,1.,0.05,1.0)
            CALL PGSCI(4)
            CALL PGMTEXT('T',3.,0.,0.,'MODEL')
            CALL PGSCI(2)
            CALL PGMTEXT('T',2.,0.,0.,'PROBLEM: '//INFILE2)
            CALL PGSCI(3)
            CALL PGMTEXT('T',1.,0.,0.,'TEMPLATE: '//INFILE1)
            CALL PGSCI(1)
            CALL PGIDEN_RED
            CALL PGSLS(2)
            CALL PGSCI(3)
            CALL PGMOVE(REAL(NF1),YMIN)
            CALL PGDRAW(REAL(NF1),YMAX)
            CALL PGMOVE(REAL(NF2),YMIN)
            CALL PGDRAW(REAL(NF2),YMAX)
            CALL PGSCI(1)
            IF(NEML.GT.0)THEN !...........................regiones a enmascarar
              CALL PGSLS(4)
              CALL PGSCI(2)
              DO IEML=1,NEML
                CALL PGMOVE(REAL(IEML1_SGL(IEML)),YMIN)
                CALL PGDRAW(REAL(IEML1_SGL(IEML)),YMAX)
                CALL PGMOVE(REAL(IEML2_SGL(IEML)),YMIN)
                CALL PGDRAW(REAL(IEML2_SGL(IEML)),YMAX)
              END DO
              CALL PGSCI(4)
              DO IEML=1,NEML
                CALL PGMOVE(REAL(IEML1_STL(IEML)),YMIN)
                CALL PGDRAW(REAL(IEML1_STL(IEML)),YMAX)
                CALL PGMOVE(REAL(IEML2_STL(IEML)),YMIN)
                CALL PGDRAW(REAL(IEML2_STL(IEML)),YMAX)
              END DO
              CALL PGSCI(1)
            END IF
            CALL PGSLS(1)
          END DO
          IF(LPAUSE)THEN
            WRITE(*,100) 'Press <CR> to continue...'
            READ(*,*)
          END IF
        END IF
C------------------------------------------------------------------------------
c se corrige de redshift galaxia y modelo (en unidades enteras de pixels); 
C válido para velocidades positivas o negativas
        DO I=1,NP
          SGL2(I)=0.0
          SML2(I)=0.0
        END DO
C
        IVEL0=INT(VELO+0.5)
C
        DO I=1,NP
          IF(((I+IVEL0).GE.1).AND.((I+IVEL0).LE.NP))THEN
            SGL2(I)=SGL(I+IVEL0)
            SML2(I)=SML(I+IVEL0)
          END IF
        END DO
C------------------------------------------------------------------------------
C se ajusta el continuo (o pseudo-continuo)
C TEMPLATE: ajustamos todo el intervalo entre NF1 y NF2
        IF(FLAGCA)THEN
          J=0
          DO I=NF1,NF2
            J=J+1
            XL(J)=X(I)
            YL(J)=STL(I)
            SIGMAY(J)=1.
          END DO
          NPTS=J
          IF(CCONT.EQ.'1')THEN
            CALL POLFIT(XL,YL,SIGMAY,NPTS,NTERMS,0,CO,CHISQR)
          ELSE
            CALL PSEUDOFIT(XL,YL,NPTS,NTERMS,YRMSTOL,
     +       PSEUDO_WEIGHT,PSEUDO_POWER,.TRUE.,CO)
          END IF
          DO I=1,NP
            YLT(I)=FPOLY(NDEG,CO,X(I))
          END DO
        END IF
C PROBLEM: en este caso vamos a ignorar las regiones que hay que enmascarar
        DO I=1,NP
          IFCHAN(I)=.TRUE.
        END DO
        IF(NEML.GT.0)THEN
          DO IEML=1,NEML
            DO I=IEML1_STL(IEML),IEML2_STL(IEML) !sist. de ref. de la template
              IF((I.GE.1).AND.(I.LE.NP))THEN
                IFCHAN(I)=.FALSE.
              END IF
            END DO
          END DO
        END IF
        J=0
        DO I=NF1,NF2
          IF(IFCHAN(I))THEN
            J=J+1
            XL(J)=X(I)
            YL(J)=SGL2(I)
            SIGMAY(J)=1.
          END IF
        END DO
        NPTS=J
        IF(CCONT.EQ.'1')THEN
          CALL POLFIT(XL,YL,SIGMAY,NPTS,NTERMS,0,CO,CHISQR)
        ELSE
          CALL PSEUDOFIT(XL,YL,NPTS,NTERMS,YRMSTOL,
     +     PSEUDO_WEIGHT,PSEUDO_POWER,.TRUE.,CO)
        END IF
        DO I=1,NP
          YLG(I)=FPOLY(NDEG,CO,X(I))
        END DO
C MODEL: ajustamos todo el intervalo entre NF1 y NF2
        J=0
        DO I=NF1,NF2
          J=J+1
          XL(J)=X(I)
          YL(J)=SML2(I)
          SIGMAY(J)=1.
        END DO
        NPTS=J
        IF(CCONT.EQ.'1')THEN
          CALL POLFIT(XL,YL,SIGMAY,NPTS,NTERMS,0,CO,CHISQR)
        ELSE
          CALL PSEUDOFIT(XL,YL,NPTS,NTERMS,YRMSTOL,
     +     PSEUDO_WEIGHT,PSEUDO_POWER,.TRUE.,CO)
        END IF
        DO I=1,NP
          YLM(I)=FPOLY(NDEG,CO,X(I))
        END DO
C------------------------------------------------------------------------------
c Plots the redshifted spectra:
        IF(ALLPLOT)THEN
          !para el PROBLEM spectrum, calculamos los límites ignorando
          !regiones a interpolar
          DO J=1,NP
            IFCHAN(J)=.TRUE.
          END DO
          IF(NEML.GT.0)THEN
            DO IEML=1,NEML
              DO J=IEML1_STL(IEML),IEML2_STL(IEML)
                IFCHAN(J)=.FALSE.
              END DO
            END DO
          END IF
          CALL FINDMML(NP,1,NP,X,XMING,XMAXG)
          CALL FINDMMLMASK(NP,1,NP,SGL2,IFCHAN,YMING,YMAXG)
          DXG=0.05*(XMAXG-XMING)
          XMING=XMING-DXG
          XMAXG=XMAXG+DXG
          DYG=0.05*(YMAXG-YMING)
          YMING=YMING-DYG
          YMAXG=YMAXG+DYG
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGPAGE
            CALL AUTOPLOT(NP,X,SGL2,1,NP, !.............................PROBLEM
     +       'channel (LOG SCALE [1:NP], UNREDSHIFTED)','rebinned flux',
     >        ' ',
     +       .FALSE.,XMING,XMAXG,YMING,YMAXG,0.05,
     +       0,.TRUE.,'BCNTS','BNTS',
     +       101,2,
     +       0.,1.,0.05,1.0)
            XMI(1)=XMING
            XMA(1)=XMAXG
            YMI(1)=YMING
            YMA(1)=YMAXG
            CALL AUTOPLOT(NP,X,STL,1,NP, !.............................TEMPLATE
     +       ' ',' ',' ',
     +       .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +       0,.FALSE.,'BCNTS','CMTS',
     +       101,3,
     +       0.,1.,0.05,1.0)
            XMI(2)=XMIN
            XMA(2)=XMAX
            YMI(2)=YMIN
            YMA(2)=YMAX
            CALL AUTOPLOT(NP,X,SML2,1,NP, !...............................MODEL
     +       ' ',' ',' ',
     +       .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +       0,.FALSE.,'BCNTS','BC',
     +       101,4,
     +       0.,1.,0.05,1.0)
            XMI(3)=XMIN
            XMA(3)=XMAX
            YMI(3)=YMIN
            YMA(3)=YMAX
            CALL PGSCI(4)
            CALL PGMTEXT('T',3.,0.,0.,'MODEL')
            CALL PGSCI(2)
            CALL PGMTEXT('T',2.,0.,0.,'PROBLEM: '//INFILE2)
            CALL PGSCI(3)
            CALL PGMTEXT('T',1.,0.,0.,'TEMPLATE: '//INFILE1)
            CALL PGSCI(1)
            CALL AUTOPLOT(NP,X,YLG,1,NP, !....................PROBLEM continuum
     +       ' ',' ',' ',
     +       .FALSE.,XMI(1),XMA(1),YMI(1),YMA(1),0.05,
     +       0,.FALSE.,' ',' ',
     +       101,2,
     +       0.,1.,0.05,1.0)
            CALL AUTOPLOT(NP,X,YLT,1,NP, !...................TEMPLATE continuum
     +       ' ',' ',' ',
     +       .FALSE.,XMI(2),XMA(2),YMI(2),YMA(2),0.05,
     +       0,.FALSE.,' ',' ',
     +       101,3,
     +       0.,1.,0.05,1.0)
            CALL AUTOPLOT(NP,X,YLM,1,NP, !......................MODEL continuum
     +       ' ',' ',' ',
     +       .FALSE.,XMI(3),XMA(3),YMI(3),YMA(3),0.05,
     +       0,.FALSE.,' ',' ',
     +       101,4,
     +       0.,1.,0.05,1.0)
            CALL PGSLS(2)
            CALL PGSCI(3)
            CALL PGMOVE(REAL(NF1),YMIN)
            CALL PGDRAW(REAL(NF1),YMAX)
            CALL PGMOVE(REAL(NF2),YMIN)
            CALL PGDRAW(REAL(NF2),YMAX)
            CALL PGSCI(1)
            IF(NEML.GT.0)THEN !...........................regiones a enmascarar
              CALL PGSLS(4)
              DO IEML=1,NEML
                CALL PGMOVE(REAL(IEML1_STL(IEML)),YMIN)
                CALL PGDRAW(REAL(IEML1_STL(IEML)),YMAX)
                CALL PGMOVE(REAL(IEML2_STL(IEML)),YMIN)
                CALL PGDRAW(REAL(IEML2_STL(IEML)),YMAX)
              END DO
            END IF
            CALL PGSLS(1)
            CALL PGIDEN_RED
          END DO
          IF(LPAUSE)THEN
            WRITE(*,100) 'Press <CR> to continue...'
            READ(*,*)
          END IF
        END IF
C------------------------------------------------------------------------------
C Eliminacion de bajas frecuencias. Tenemos 2 posibilidades:
C En la primera (LREPLACE_CONTINUUM_REGIONS=.TRUE.) cada espectro lo dividimos 
C por su continuo. En la segunda (LREPLACE_CONTINUUM_REGIONS=.FALSE.) ponemos a
C todos los espectros el continuo de la galaxia y luego ajustamos una recta y 
C se la quitamos.
        IF(LREPLACE_CONTINUUM_REGIONS)THEN
          DO I=1,NP
            STL2(I)=STL(I)/YLT(I)
            SGL2(I)=SGL2(I)/YLG(I)
            SML2(I)=SML2(I)/YLM(I)
          END DO
        ELSE
c Hacemos que el continuo de modelo y estrella coincidan con el de la galaxia
          DO I=1,NP
            STL2(I)=STL(I)/YLT(I)*YLG(I)
            SML2(I)=SML2(I)/YLM(I)*YLG(I)
          END DO
C y ajustamos una recta al espectro de la galaxia y dividimos los tres 
c espectros por el
          J=0
          DO I=NF1,NF2
            IF(IFCHAN(I))THEN
              J=J+1
              XL(J)=X(I)
              YL(J)=SGL2(I)
              SIGMAY(J)=1.
            END IF
          END DO
          NPTS=J
          CALL POLFIT(XL,YL,SIGMAY,NPTS,2,0,CO,CHISQR)
          DO I=1,NP
            YR=CO(1)+X(I)*CO(2)
            STL2(I)=STL2(I)/YR
            SGL2(I)=SGL2(I)/YR
            SML2(I)=SML2(I)/YR
          END DO
        END IF
C------------------------------------------------------------------------------
C Extraemos los 2**x=NCH pixeles utiles y restamos la media
        MEANT=0.D0
        MEANG=0.D0
        MEANM=0.D0
        DO I=1,NCH
          ST(I)=STL2(I+NF1-1)
          SG(I)=SGL2(I+NF1-1)
          SM(I)=SML2(I+NF1-1)
          MEANT=MEANT+DBLE(ST(I))
          MEANG=MEANG+DBLE(SG(I))
          MEANM=MEANM+DBLE(SM(I))
        END DO
        MEANT=MEANT/DBLE(NCH)
        MEANG=MEANG/DBLE(NCH)
        MEANM=MEANM/DBLE(NCH)
        DO I=1,NCH
          ST(I)=ST(I)-REAL(MEANT)
          SG(I)=SG(I)-REAL(MEANG)
          SM(I)=SM(I)-REAL(MEANM)
        END DO
C------------------------------------------------------------------------------
C Se introduce gamma (si se hace antes se destruye en la normalizacion)
        DO I=1,NCH
          SM(I)=SM(I)*GAM
        END DO
C------------------------------------------------------------------------------
C Para dibujar al final del programa
        IF(ISIMUL.EQ.0)THEN
          DO I=1,NCH
            SGPLOT(I)=SG(I)
            SMPLOT(I)=SM(I)
          END DO
        END IF
C------------------------------------------------------------------------------
C Para eliminar regiones problemáticas se substituye el modelo por la galaxia
c en las zonas elegidas
        IF(NEML.GT.0)THEN
          DO IEML=1,NEML
            DO J=IEML1_NCH(IEML),IEML2_NCH(IEML)
              SG(J)=SM(J)
            END DO
          END DO
        END IF
C------------------------------------------------------------------------------
        IF(ALLPLOT)THEN
          CALL FINDMML(NCH,1,NCH,SG,YMIN,YMAX)
          CALL FINDMML(NCH,1,NCH,ST,YMING,YMAXG)
          YMIN=AMIN1(YMIN,YMING)
          YMAX=AMAX1(YMAX,YMAXG)
          CALL FINDMML(NCH,1,NCH,SM,YMING,YMAXG)
          YMIN=AMIN1(YMIN,YMING)
          YMAX=AMAX1(YMAX,YMAXG)
          DYG=(YMAX-YMIN)*0.05
          YMIN=YMIN-DYG
          YMAX=YMAX+DYG
          XMIN=-1.*REAL(NCH)/100.
          XMAX=REAL(NCH)*1.01
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGPAGE
            CALL AUTOPLOT(NCH,X,SG,1,NCH, !.............................PROBLEM
     +       'channel (USEFUL LOG SCALE [1:NCH])',
     +        'flux ZERO-MEAN DETRENDED',' ',
     +       .FALSE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +       0,.TRUE.,'BCNTS','BCNTS',
     +       101,2,
     +       0.,1.,0.05,1.0)
            CALL AUTOPLOT(NCH,X,ST,1,NCH, !............................TEMPLATE
     +       ' ',' ',' ',
     +       .FALSE.,XMIN,XMAX,YMIN,YMAX,0.00,
     +       0,.FALSE.,' ',' ',
     +       101,3,
     +       0.,1.,0.05,1.0)
            CALL AUTOPLOT(NCH,X,SM,1,NCH, !...............................MODEL
     +       ' ',' ',' ',
     +       .FALSE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +       0,.FALSE.,' ',' ',
     +       101,4,
     +       0.,1.,0.05,1.0)
            CALL PGSLS(2)
            CALL PGMOVE(XMIN,0.)
            CALL PGDRAW(XMAX,0.)
            CALL PGSLS(1)
            CALL PGSCI(4)
            CALL PGMTEXT('T',3.,0.,0.,'MODEL')
            CALL PGSCI(2)
            CALL PGMTEXT('T',2.,0.,0.,'PROBLEM: '//INFILE2)
            CALL PGSCI(3)
            CALL PGMTEXT('T',1.,0.,0.,'TEMPLATE: '//INFILE1)
            CALL PGSCI(1)
            IF(NEML.GT.0)THEN !...........................regiones a enmascarar
              CALL PGSLS(4)
              DO IEML=1,NEML
                CALL PGMOVE(REAL(IEML1_NCH(IEML)),YMIN)
                CALL PGDRAW(REAL(IEML1_NCH(IEML)),YMAX)
                CALL PGMOVE(REAL(IEML2_NCH(IEML)),YMIN)
                CALL PGDRAW(REAL(IEML2_NCH(IEML)),YMAX)
              END DO
              CALL PGSLS(1)
            END IF
            CALL PGIDEN_RED
          END DO
          IF(LPAUSE)THEN
            WRITE(*,100) 'Press <CR> to continue...'
            READ(*,*)
          END IF
        END IF
C------------------------------------------------------------------------------
C Cosbell the spectra (fraction of spectra should be 1/16, 1/8, 1/4,...)
        IF((FLAGCA).AND.(ISIMUL.EQ.0).AND.(NITEMP.LE.1).AND.
     +   (IWORKING_SCAN.EQ.0))THEN
          WRITE(*,*)
          WRITE(*,101)'* Fraction f of spectra to be cosbelled in '//
     +     'the form 1/f (ex. 1/16, 1/8, 1/4):'
          WRITE(*,100)'fraction '
          IFCOSB=READILIM('8',1,NCH)
          NL=NINT(REAL(NCH)/REAL(IFCOSB))
          WRITE(CDUMMY,*) NL
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          WRITE(*,101)'>>> Cosine bell will be applied over the first'//
     +     ' and last '//CDUMMY(1:L)//' channels'
          DO J=1,NL
            COSBELL(J)=0.5*(1.-COS(PI*REAL(J)/REAL(NL)))
          END DO
          DO J=NL+1,NCH-NL
            COSBELL(J)=1.0
          END DO
          DO J=NCH-NL+1,NCH
            COSBELL(J)=0.5*(1.-COS(PI*REAL(NCH-J)/REAL(NL)))
          END DO
        END IF
C..............................................................................
C aplicamos la campana de coseno
        DO I=1,NCH
          SG(I)=SG(I)*COSBELL(I)
          ST2(I)=ST(I)*COSBELL(I)
          SM(I)=SM(I)*COSBELL(I)
        END DO
C..............................................................................
C dibujamos
        IF(ALLPLOT)THEN
          CALL FINDMML(NCH,1,NCH,SG,YMIN,YMAX)
          CALL FINDMML(NCH,1,NCH,ST2,YMING,YMAXG)
          YMIN=AMIN1(YMIN,YMING)
          YMAX=AMAX1(YMAX,YMAXG)
          CALL FINDMML(NCH,1,NCH,SM,YMING,YMAXG)
          YMIN=AMIN1(YMIN,YMING)
          YMAX=AMAX1(YMAX,YMAXG)
          DYG=(YMAX-YMIN)*0.05
          YMIN=YMIN-DYG
          YMAX=YMAX+DYG
          XMIN=-1.*REAL(NCH)/100.
          XMAX=REAL(NCH)*1.01
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGPAGE
            CALL AUTOPLOT(NCH,X,SG,1,NCH, !.............................PROBLEM
     +       'channel (USEFUL LOG SCALE [1:NCH]',
     +        'flux COSBELLED',' ',
     +       .FALSE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +       0,.TRUE.,'BCNTS','BCNTS',
     +       101,2,
     +       0.,1.,0.05,1.0)
            CALL AUTOPLOT(NCH,X,ST2,1,NCH, !...........................TEMPLATE
     +       ' ',' ',' ',
     +       .FALSE.,XMIN,XMAX,YMIN,YMAX,0.00,
     +       0,.FALSE.,' ',' ',
     +       101,3,
     +       0.,1.,0.05,1.0)
            CALL AUTOPLOT(NCH,X,SM,1,NCH, !...............................MODEL
     +       ' ',' ',' ',
     +       .FALSE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +       0,.FALSE.,' ',' ',
     +       101,4,
     +       0.,1.,0.05,1.0)
            CALL PGSLS(2)
            CALL PGMOVE(XMIN,0.)
            CALL PGDRAW(XMAX,0.)
            CALL PGSLS(1)
            CALL AUTOPLOT(NCH,X,COSBELL,1,NCH, !........................COSBELL
     +       ' ',' ',' ',
     +       .FALSE.,XMIN,XMAX,-0.1,1.1,0.05,
     +       0,.FALSE.,' ',' ',
     +       101,5,
     +       0.,1.,0.05,1.0)
            CALL PGSCI(2)
            CALL PGMTEXT('T',2.,0.,0.,'PROBLEM: '//INFILE2)
            CALL PGSCI(3)
            CALL PGMTEXT('T',1.,0.,0.,'TEMPLATE: '//INFILE1)
            CALL PGSCI(4)
            CALL PGMTEXT('T',3.,0.,0.,'MODEL')
            CALL PGSCI(5)
            CALL PGMTEXT('T',3.,1.,1.,'COSBELL')
            CALL PGSCI(1)
            IF(NEML.GT.0)THEN !...........................regiones a enmascarar
              CALL PGSLS(4)
              DO IEML=1,NEML
                CALL PGMOVE(REAL(IEML1_NCH(IEML)),YMIN)
                CALL PGDRAW(REAL(IEML1_NCH(IEML)),YMAX)
                CALL PGMOVE(REAL(IEML2_NCH(IEML)),YMIN)
                CALL PGDRAW(REAL(IEML2_NCH(IEML)),YMAX)
              END DO
              CALL PGSLS(1)
            END IF
            CALL PGIDEN_RED
          END DO
          IF(LPAUSE)THEN
            WRITE(*,100) 'Press <CR> to continue...'
            READ(*,*)
          END IF
        END IF
C------------------------------------------------------------------------------
C calculamos transformadas de Fourier
        IMODE=1 !directa
        DO I=1,NCH
           SGR(I)=SG(I)
           SGI(I)=0.
           STR(I)=ST2(I)
           STI(I)=0.
           SMR(I)=SM(I)
           SMI(I)=0.
        END DO
        CALL CFFT(NCH,SGR,SGI,IMODE)
        CALL CFFT(NCH,STR,STI,IMODE)
        CALL CFFT(NCH,SMR,SMI,IMODE)
C------------------------------------------------------------------------------
C calculamos los espectros de potencias
        DO I=1,NCH
          PG(I)=DBLE(SGR(I))*DBLE(SGR(I))+DBLE(SGI(I))*DBLE(SGI(I))
          IF(PG(I).GT.0.D0)THEN
            PGL(I)=REAL(DLOG10(PG(I)))
          ELSE
            PGL(I)=PGL(I-1)
          END IF
          PT(I)=DBLE(STR(I))*DBLE(STR(I))+DBLE(STI(I))*DBLE(STI(I))
          IF(PT(I).GT.0.D0)THEN
            PTL(I)=REAL(DLOG10(PT(I)))
          ELSE
            PTL(I)=PTL(I-1)
          END IF
          PM(I)=DBLE(SMR(I))*DBLE(SMR(I))+DBLE(SMI(I))*DBLE(SMI(I))
          IF(PM(I).GT.0.D0)THEN
            PML(I)=REAL(DLOG10(PM(I)))
          ELSE
            PML(I)=PML(I-1)
          END IF
        END DO
C------------------------------------------------------------------------------
C Calculo del filtro de Wiener:
C..............................................................................
C Noise: (tomamos la mediana de la segunda mitad del espectro de frecuencias)
        DO I=NCH/4+1,NCH/2
          PNOISE(I-NCH/4)=PG(I)-PM(I)
        END DO
        PNOISE0=FMEDIAN1(NCH/4,PNOISE)
        IF(PNOISE0.LT.0.)THEN
          DO I=NCH/4+1,NCH/2
            PNOISE(I-NCH/4)=PG(I)
          END DO
          PNOISE0=FMEDIAN1(NCH/4,PNOISE)
          WRITE(*,101) 'WARNING: error in the determination of noise '//
     +     '(noise forced to 1.E-10)'
          PNOISE0=1.E-10
        END IF
C..............................................................................
C Signal: (el ajuste es en escala logaritmica)
        XF(NCH/4)=0.
        YF(NCH/4)=PML(1)
        DO I=1,NCH/4-1
          XF(NCH/4+I)=REAL(I)
          XF(NCH/4-I)=-1.*REAL(I)
          YF(NCH/4+I)=PML(I+1)
          YF(NCH/4-I)=PML(I+1)
        END DO
        NP0=NCH/2-1
        CALL GAUSCFIT_MOVEL(XGAUS,SIGGAUS,AMPGAUS,YGAUS,
     +   EEX0,EESIGMA,EEAMP,EEY0,YRMSTOL) 
        DO I=1,NCH/2
          J=I-1
          YL(I)=YGAUS+AMPGAUS*EXP(-1.*((REAL(J)-XGAUS)*(REAL(J)-XGAUS)/
     +      (2.*SIGGAUS*SIGGAUS)))
        END DO
C..............................................................................
C Filter:
        DO I=1,NCH/2
          FD=10.D0**(2.D0*DBLE(YL(I)))
          WIEN(I)=REAL(FD/(FD+DBLE(PNOISE0)*DBLE(PNOISE0)))
        END DO
C..............................................................................
C El filtro es cero para frecuencias correspondientes a menos de LCUT Angstroms
        IF((FLAGCA).AND.(ISIMUL.EQ.0).AND.(NITEMP.LE.1).AND.
     +   (IWORKING_SCAN.EQ.0))THEN
          WRITE(*,*)
          WRITE(*,101)'* Low frequency cut-off:'
          DO I=1,10 !display first 10 frequencies
            WRITE(CDUMMY,*) DLAM/REAL(I)
            CALL RMBLANK(CDUMMY,CDUMMY,L)
            WRITE(*,'(A,I2.2,A)') 'frequency #',I,' --> '//
     +       CDUMMY(1:L)//' Angstroms'
          END DO
          WRITE(*,100)'Introduce maximum wavelength (in A) '
          WRITE(CDUMMY,*) LCUT
          LCUT=READF(CDUMMY)
        END IF
        LOOP=.TRUE.
        I=0
        DO WHILE(LOOP)
          I=I+1
          IF(DLAM/REAL(I).GT.LCUT)THEN
            WIEN(I)=0.
          ELSE
            LOOP=.FALSE.
          END IF
          IF(I.EQ.NCH/2) LOOP=.FALSE.
        END DO
        IF((FLAGCA).AND.(ISIMUL.EQ.0).AND.(NITEMP.EQ.1).AND.
     +   (IWORKING_SCAN.EQ.0))THEN
          WRITE(*,'(A,I6)') 'Low frequency cut-off: ',I
        END IF
C..............................................................................
C Calculamos la frecuencia minima y maxima del filtro
        KMIN=I
        LOOP=.TRUE.
        DO WHILE(LOOP)
          I=I+1
          IF(I.EQ.NCH/2)THEN
            LOOP=.FALSE.
            IF(WIEN(I).GT.0) I=I+1
          END IF
          IF(WIEN(I).LE.0.) LOOP=.FALSE.
        END DO
        KMAX=I-1
        IF((FLAGCA).AND.(ISIMUL.EQ.0).AND.(NITEMP.LE.1).AND.
     +   (IWORKING_SCAN.EQ.0))THEN
          WRITE(*,'(A,I6)') 'Maximum frequency....: ',KMAX
        END IF
C------------------------------------------------------------------------------
C Plots of power spectra and filter
        IF(ALLPLOT)THEN
          CALL FINDMML(NCH/2,1,NCH/2,PML,YMIN,YMAX)
          CALL FINDMML(NCH/2,1,NCH/2,PGL,YMING,YMAXG)
          YMIN=AMIN1(YMIN,YMING)
          YMAX=AMAX1(YMAX,YMAXG)
          DYG=(YMAX-YMIN)/20.
          YMAX=YMAX+DYG
          YMIN=YMIN-DYG
          XMIN=-1.*REAL(NCH)/100.
          XMAX=REAL(NCH/2)*1.01
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGPAGE
C
            CALL AUTOPLOT(NCH/2,X,PGL,1,NCH/2, !..................PROBLEM power
     +       'frequency','log P',' ',
     +       .FALSE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +       0,.TRUE.,'BCNTS','BCNTS',
     +       101,2,
     +       0.,1.,0.05,1.0)
            CALL AUTOPLOT(NCH/2,X,PML,1,NCH/2, !....................MODEL power
     +       ' ',' ',' ',
     +       .FALSE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +       0,.FALSE.,' ',' ',
     +       101,4,
     +       0.,1.,0.05,1.0)
            CALL PGSCI(2)
            CALL PGMTEXT('T',2.,0.,0.,'PROBLEM: '//INFILE2)
            CALL PGSCI(4)
            CALL PGMTEXT('T',1.,0.,0.,'MODEL')
            CALL PGSCI(1)
            CALL PGMTEXT('T',2.,1.,1.,'WIENER FILTER')
C
            CALL PGSCI(6)
            CALL PGMOVE(0.,ALOG10(PNOISE0))
            CALL PGDRAW(REAL(NCH/2),ALOG10(PNOISE0)) !....................Noise
            CALL PGSCI(1)
C
            CALL PGSCI(5)
            CALL PGLINE(NCH/2,X,YL) !..............................ajuste señal
            CALL PGSCI(1)
C
            YMING=0.
            YMAXG=1.1
            CALL AUTOPLOT(NCH/2,X,WIEN,1,NCH/2, !.................filtro WIENER
     +       ' ',' ',' ',
     +       .FALSE.,XMIN,XMAX,YMING,YMAXG,0.05,
     +       0,.FALSE.,' ',' ',
     +       101,8,
     +       0.,1.,0.05,1.0)
C
            CALL PGSCI(1)
            CALL PGSLS(4)
            CALL PGMOVE(0.,1.)
            CALL PGDRAW(REAL(NCH/2),1.)
            CALL PGSLS(1)
            CALL PGIDEN_RED
          END DO
          IF(LPAUSE)THEN
            WRITE(*,100) 'Press <CR> to continue...'
            READ(*,*)
          END IF
        END IF
C------------------------------------------------------------------------------
C Ajuste para el modelo
        NVAR=3
        DXX0(1)=0.1
        DXX0(2)=1.
        DXX0(3)=1.
        DO I=1,NCH/2
          SR(I)=SMR(I)
          SI(I)=SMI(I)
        END DO
        XX0(1)=GAM
        XX0(2)=REAL(IVEL0)-VELO
        XX0(3)=SIGC
        CALL DOWNHILL(NVAR,XX0,DXX0,CHISQG,1.0,0.5,2.0,YRMSTOL,
     +   XX,DXX,NEVAL)
        IF(FLAGCA.AND.ISIMUL.EQ.0.AND.NITEMP.LE.1)THEN
          WRITE(*,*)
          WRITE(*,101) 'RESULTS FOR THE MODEL:'
          WRITE(*,100) 'GAMMA: '
          WRITE(*,*) XX(1)
          WRITE(*,100) 'VEL..: '
!         WRITE(*,*) (REAL(IVEL0)-XX(2))*REAL(DV)                          !mal
          RCVEL1=EXP((REAL(IVEL0)-XX(2))*DLOGW)         !correccion relativista
          WRITE(*,*) (RCVEL1*RCVEL1-1.0)/(RCVEL1*RCVEL1+1.0)*C
          WRITE(*,100) 'SIGMA: '
          WRITE(*,*) XX(3)*REAL(DV)
        END IF
C------------------------------------------------------------------------------
c Calcula la transformada de Fourier de la funcion de ensanchamiento con los
c parametros teoricos
        XTH(1)=GAM
!       XTH(2)=REAL(IVEL0)-VEL/REAL(DV)                                    !mal
        BETA=VEL/C
        RCVEL1=(1.0+BETA)/SQRT(1.0-BETA*BETA)            !correcion relativista
        XTH(2)=REAL(IVEL0)-ALOG(RCVEL1)/DLOGW
        XTH(3)=SIG/REAL(DV)
        CALL FFTBROAD(NCH,XTH)
        DO I=1,NCH/2
           B0R(I)=FFTA(I)
           B0I(I)=FFTB(I)
        END DO
c Dibuja los ajustes. En particular dibuja las partes reales e imaginarias
c del cociente F(M)/F(T) frente a las componentes de la transformada de la
c funcion de ensanchamiento:
        IF(ALLPLOT)THEN
C Calcula la transformada de Fourier de la funcion de ensanchamiento con los
c parametros ajustados
           CALL FFTBROAD(NCH,XX)
           DO I=1,NCH/2
              B1R(I)=FFTA(I)
              B1I(I)=FFTB(I)
           END DO
c Solo dibujamos hasta la frecuencia en la que el filtro es > 0.001
          LOOP=.TRUE.
          I=KMIN-1
          DO WHILE(LOOP)
            I=I+1
            IF(WIEN(I).LT.0.001)THEN
              LOOP=.FALSE.
            END IF
            IF(I.EQ.KMAX)THEN
              IF(LOOP) I=I+1
              LOOP=.FALSE.
            END IF
          END DO
          KMAXPL=I-1
c Cociente:
          DO I=1,KMAXPL
             DBR=DBLE(STR(I))*DBLE(SR(I))+DBLE(STI(I))*DBLE(SI(I))
             DBI=DBLE(STR(I))*DBLE(SI(I))-DBLE(STI(I))*DBLE(SR(I))
             QR(I)=REAL(DBR/PT(I))
             QI(I)=REAL(DBI/PT(I))
          END DO
C
          DO ITERM=NTERM,1,-1
C..............................................................................
            CALL PGSLCT(IDN(ITERM))
            CALL PGPAGE
C..............................................................................
            CALL AUTOPLOT(KMAXPL,X,QI,1,KMAXPL, !...................I(Quotient)
     +       'frequency','I(Quotient)',' ',
     +       .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +       0,.TRUE.,'BCNTS','BCNTS',
     +       101,1,
     +       0.,1.,0.05,0.65)
            CALL PGSCI(4)
            CALL PGMOVE(REAL(KMIN),YMIN)
            CALL PGDRAW(REAL(KMIN),YMAX)
            CALL PGSCI(2)
            CALL PGLINE(NCH/2,X,B1I) !...................................fitted
            CALL PGSCI(3)
            CALL PGLINE(NCH/2,X,B0I) !..............................theoretical
            CALL PGSCI(1)
C..............................................................................
            CALL AUTOPLOT(KMAXPL,X,QR,1,KMAXPL, !...................R(Quotient)
     +       ' ','R(Quotient)',' ',
     +       .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +       0,.FALSE.,'BCTS','BCNTS',
     +       101,1,
     +       0.,1.,0.40,1.00)
            CALL PGSCI(4)
            CALL PGMOVE(REAL(KMIN),YMIN)
            CALL PGDRAW(REAL(KMIN),YMAX)
            CALL PGSCI(2)
            CALL PGLINE(NCH/2,X,B1R) !...................................fitted
            CALL PGSCI(3)
            CALL PGLINE(NCH/2,X,B0R) !..............................theoretical
            CALL PGSCI(1)
C..............................................................................
            CALL PGSCI(2)
            CALL PGMTEXT('T',1.,0.,0.,'FITTED')
            CALL PGSCI(3)
            CALL PGMTEXT('T',1.,1.,1.,'THEORETICAL')
            CALL PGSCI(1)
            CALL PGMTEXT('T',1.,0.5,0.5,'FIT TO THE MODEL')
            CALL PGIDEN_RED
C..............................................................................
          END DO
          IF(LPAUSE)THEN
            WRITE(*,100) 'Press <CR> to continue...'
            READ(*,*)
          END IF
        END IF
C------------------------------------------------------------------------------
C Calcula los residuos como la transformada del modelo menos 
C la esperada (es decir el producto de las trans. de estrella y funcion de
C ensanchamiento teorica)-> Res = F(M) - F(Bteor)*F(Templ)
        DO I=1,NCH/2
          DBR=DBLE(B0R(I))*DBLE(STR(I))-DBLE(B0I(I))*DBLE(STI(I))
          DBI=DBLE(B0R(I))*DBLE(STI(I))+DBLE(B0I(I))*DBLE(STR(I))
          RESR(I)=REAL(DBLE(SMR(I))-DBR)
          RESI(I)=REAL(DBLE(SMI(I))-DBI)
        END DO
C y se le aplica a la transformada de Fourier de la galaxia:
        DO I=1,NCH/2
          SGR2(I)=SGR(I)-RESR(I)
          SGI2(I)=SGI(I)-RESI(I)
        END DO
C------------------------------------------------------------------------------
C Dibuja espectros de potencias
        IF(ALLPLOT)THEN
          DO I=1,NCH/2
            PG2(I)=DBLE(SGR2(I))*DBLE(SGR2(I))+
     +       DBLE(SGI2(I))*DBLE(SGI2(I))
            IF(PG2(I).GT.0.D0)THEN
              PGL2(I)=REAL(DLOG10(PG2(I)))
            ELSE
              PGL2(I)=PGL2(I-1)
            END IF
          END DO
C..............................................................................
          CALL FINDMML(NCH/2,1,NCH/2,PTL,YMIN,YMAX)
          CALL FINDMML(NCH/2,1,NCH/2,PGL2,YMING,YMAXG)
          YMIN=AMIN1(YMIN,YMING)
          YMAX=AMAX1(YMAX,YMAXG)
          DYG=(YMAX-YMIN)/20.
          YMAX=YMAX+DYG
          YMIN=YMIN-DYG
          XMIN=-1.*REAL(NCH)/100.
          XMAX=REAL(NCH/2)*1.01
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGPAGE
            CALL AUTOPLOT(NCH/2,X,PGL2,1,NCH/2, !.................PROBLEM power
     +       'frequency','log P',' ',
     +       .FALSE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +       0,.TRUE.,'BCNTS','BCNTS',
     +       101,2,
     +       0.,1.,0.05,1.0)
            CALL AUTOPLOT(NCH/2,X,PTL,1,NCH/2, !.................TEMPLATE power
     +       ' ',' ',' ',
     +       .FALSE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +       0,.FALSE.,' ',' ',
     +       101,3,
     +       0.,1.,0.05,1.0)
            YMING=0.
            YMAXG=1.1
            CALL AUTOPLOT(NCH/2,X,WIEN,1,NCH/2, !.................filtro WIENER
     +       ' ',' ',' ',
     +       .FALSE.,XMIN,XMAX,YMING,YMAXG,0.05,
     +       0,.FALSE.,' ','CMTS',
     +       101,8,
     +       0.,1.,0.05,1.0)
            CALL PGSCI(1)
            CALL PGSLS(4)
            CALL PGMOVE(0.,1.)
            CALL PGDRAW(REAL(NCH/2),1.)
            CALL PGSLS(1)
            CALL PGSCI(2)
            CALL PGMTEXT('T',2.,0.,0.,'PROBLEM: '//INFILE2)
            CALL PGSCI(3)
            CALL PGMTEXT('T',1.,0.,0.,'TEMPLATE: '//INFILE1)
            CALL PGSCI(1)
            CALL PGMTEXT('T',2.,1.,1.,'POWER SPECTRA')
            CALL PGIDEN_RED
          END DO
          IF(LPAUSE)THEN
            WRITE(*,100) 'Press <CR> to continue...'
            READ(*,*)
          END IF
        END IF
C------------------------------------------------------------------------------
C Ajuste final
        NVAR=3
        DXX0(1)=0.1
        DXX0(2)=1.
        DXX0(3)=1.
        DO I=1,NCH/2
          SR(I)=SGR2(I)
          SI(I)=SGI2(I)
        END DO
        XX0(1)=GAM
        XX0(2)=REAL(IVEL0)-VELO
        XX0(3)=SIGC
        CALL DOWNHILL(NVAR,XX0,DXX0,CHISQG,1.0,0.5,2.0,YRMSTOL,
     +   XX,DXX,NEVAL)
        GAMF=XX(1)
!       VELF=(REAL(IVEL0)-XX(2))*REAL(DV)                                  !mal
        RCVEL1=EXP((REAL(IVEL0)-XX(2))*DLOGW)           !correccion relativista
        VELF=(RCVEL1*RCVEL1-1.0)/(RCVEL1*RCVEL1+1.0)*C
        SIGF=XX(3)*REAL(DV)
        IF(NEVAL.EQ.-1) WRITE(*,100) 'ERROR IN DOWNHILL. NO SOLUTION'
C------------------------------------------------------------------------------
c Dibuja el ajuste
        IF(ALLPLOT)THEN
          !Calcula la transformada de Fourier de la funcion de ensanchamiento 
          !con los parametros ajustados
          CALL FFTBROAD(NCH,XX)
          DO I=1,NCH/2
            B2R(I)=FFTA(I)
            B2I(I)=FFTB(I)
          END DO
          !Cociente:
          DO I=1,KMAXPL
            DBR=DBLE(STR(I))*DBLE(SR(I))+DBLE(STI(I))*DBLE(SI(I))
            DBI=DBLE(STR(I))*DBLE(SI(I))-DBLE(STI(I))*DBLE(SR(I))
            QR(I)=REAL(DBR/PT(I))
            QI(I)=REAL(DBI/PT(I))
          END DO
C..............................................................................
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGPAGE
            CALL AUTOPLOT(KMAXPL,X,QI,1,KMAXPL, !...................I(Quotient)
     +       'frequency','I(Quotient)',' ',
     +       .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +       0,.TRUE.,'BCNTS','BCNTS',
     +       101,1,
     +       0.,1.,0.05,0.65)
C
            CALL PGSCI(4)
            CALL PGMOVE(REAL(KMIN),YMIN)
            CALL PGDRAW(REAL(KMIN),YMAX)
            CALL PGSCI(2)
            CALL PGLINE(NCH/2,X,B2I)
            CALL PGSCI(1)
C
            CALL AUTOPLOT(KMAXPL,X,QR,1,KMAXPL, !...................R(Quotient)
     +       ' ','R(Quotient)',' ',
     +       .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +       0,.FALSE.,'BCTS','BCNTS',
     +       101,1,
     +       0.,1.,0.40,1.00)
C
            CALL PGSCI(4)
            CALL PGMOVE(REAL(KMIN),YMIN)
            CALL PGDRAW(REAL(KMIN),YMAX)
            CALL PGSCI(2)
            CALL PGLINE(NCH/2,X,B2R)
            CALL PGSCI(1)
            CALL PGMTEXT('T',1.0,1.0,1.0,'FINAL FIT')
            WRITE(CINFO,'(3(A,I8))')
     +       'ISCAN=',IWORKING_SCAN,',_IITER=',IITER,
     +       ',_ITEMPL=',ITEMPL
            CALL RMBLANK(CINFO,CINFO,LINFO)
            LOOP=.TRUE. !eliminamos pseudo-espacios en blanco
            DO WHILE(LOOP)
              LFOUND=INDEX(CINFO(1:LINFO),'_')
              IF(LFOUND.EQ.0)THEN
                LOOP=.FALSE.
              ELSE
                CINFO(LFOUND:LFOUND)=' '
              END IF
            END DO
            CALL PGMTEXT('T',1.0,0.0,0.0,CINFO(1:LINFO))
            CALL PGIDEN_RED
          END DO
          IF(LPAUSE)THEN
            WRITE(*,100) 'Press <CR> to continue...'
            READ(*,*)
          END IF
        END IF
C------------------------------------------------------------------------------
C Mostramos resultados de la iteración actual
        IF(IITER.EQ.0.OR.NITER.EQ.-1)THEN
          WRITE(*,*)
          WRITE(*,100)
     +     '#ISCAN GAMMA     VEL       SIG     IITER  ITEMPL ISIMUL'
          IF(NTEMPLATES.GT.1)THEN
            DO II=1,NTEMPLATES
              WRITE(*,'(2X,A,I2.2,$)') 'TEMPW',II
            END DO
          END IF
          WRITE(*,*)
          IF(LWRITE_LOGFILE)THEN
            WRITE(27,100)
     +       '#ISCAN GAMMA     VEL       SIG     IITER  ITEMPL ISIMUL'
            IF(NTEMPLATES.GT.1)THEN
              DO II=1,NTEMPLATES
                WRITE(27,'(2X,A,I2.2,$)') 'TEMPW',II
              END DO
            END IF
            WRITE(27,*)
          END IF
        END IF
        WRITE(*,'(1X,I5,1X,F5.3,3X,F8.2,3X,F6.2,3(3X,I4),$)') 
     +   IWORKING_SCAN,GAMF,VELF,SIGF,IITER,ITEMPL,ISIMUL
        IF(NTEMPLATES.GT.1)THEN
          WRITE(*,100)' '
          DO II=1,NTEMPLATES
            WRITE(*,'(4X,F5.3,$)') WEI(II)
          END DO
        END IF
        WRITE(*,*)
        IF(LWRITE_LOGFILE)THEN
          WRITE(27,'(1X,I5,1X,F5.3,3X,F8.2,3X,F6.2,3(3X,I4),$)') 
     +     IWORKING_SCAN,GAMF,VELF,SIGF,IITER,ITEMPL,ISIMUL
          IF(NTEMPLATES.GT.1)THEN
            WRITE(27,100)' '
            DO II=1,NTEMPLATES
              WRITE(27,'(4X,F5.3,$)') WEI(II)
            END DO
          END IF
          WRITE(27,*)
        END IF
C------------------------------------------------------------------------------
C Decidimos si hay que iterar más o no
        LITERATE_MORE=.FALSE.
        IF(NITER.EQ.-1)THEN !..............................................FREE
          !el usuario decide
          WRITE(*,100)'Iterate (y/n) '
          CITER(1:1)=READC('y','yn')
          IF(CITER.EQ.'y') LITERATE_MORE=.TRUE.
        ELSEIF(NITER.EQ.0)THEN !...........................................AUTO
          IF(IITER.EQ.0)THEN
            !hay que iterar al menos una vez; de lo contrario no podemos 
            !calcular ningún DELTA_VEL ni DELTA_SIG
            LITERATE_MORE=.TRUE.
          ELSE
            !comparamos con los valores de la iteración anterior
            DELTA_VEL=ABS(VELF-VELF0)
            DELTA_SIG=ABS(SIGF-SIGF0)
            !si todavía no hemos alcanzado la precisión requerida, iteramos
            IF((DELTA_VEL.GT.VEL_PRECISION).OR.
     +       (DELTA_SIG.GT.VEL_PRECISION)) LITERATE_MORE=.TRUE.
          END IF
        ELSE !............................................................FIXED
          !tenemos un número fijo de iteraciones NITER
          IF(IITER.LT.NITER) LITERATE_MORE=.TRUE.
        END IF
C Proteccion para evitar iterar indefinidamente
        IF(LITERATE_MORE)THEN
          IF(IITER.GE.MAXITER) LITERATE_MORE=.FALSE.
        END IF
C------------------------------------------------------------------------------
C Seguimos iterando
        IF(LITERATE_MORE)THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!(NCL) el trozo de código que sigue está comentado porque no me gusta influir 
!de esta forma en el proceso de cálculo.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         IF(ISIMUL.GT.0.AND.IITER.GT.2)THEN
!           !esto es para evitar valores extremos; en ese caso se repite la
!           !simulacion
!           !(NCL) No entiendo por qué VELF.LT.0. en la siguiente línea
!           !      Me preocupa que esto no ayude a calcular velocidades
!           !      negativas
!           IF(SIGF/DV.LT.0.02.OR.SIGF.GT.500..OR.VELF.LT.0.)THEN
!             print*,'WARNING(NCL): pasando por VELF.LT.0'
!             NITEMP=0
!             IF(NTEMPLATES.GT.1) LOPTIMAL_TEMPLATE=.TRUE.
!             GOTO 333
!           END IF
!         END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          FLAGIT=.FALSE.
          IITER=IITER+1
          VELF0=VELF
          SIGF0=SIGF
          GOTO 222
        END IF
C------------------------------------------------------------------------------
C Ya hemos calculado GAMF, VELF y SIGF dentro de la precision requerida.
C Estos valores se han obtenido con la template óptima utilizada. Ahora
C estudiamos si hay que volver a calcular la template óptima otra vez.
        LREPITE_TEMPLATE=.TRUE.
        IF(NTEMPLATES.EQ.1)THEN !........................no hay template óptima
          LREPITE_TEMPLATE=.FALSE.
        ELSE !......................................puede haber template óptima
          IF((ISIMUL.GT.0).OR.(IWORKING_SCAN.GT.0))THEN
            !solo se recalcula la template óptima cuando estamos en el scan 
            !central y cuando, además, no estamos simulando errores
            LREPITE_TEMPLATE=.FALSE.
          ELSE
            IF((NITER.GT.0).AND.((NITEMP+1)/2.GE.NITER))THEN
              !hemos alcanzando el numero máximo de iteraciones
              LREPITE_TEMPLATE=.FALSE.
            ELSE
              !proteccion para evitar iterar indefinidamente
              IF((NITEMP+1)/2.EQ.MAXITER) LREPITE_TEMPLATE=.FALSE.
            END IF
          END IF
        END IF
C
        IF(LREPITE_TEMPLATE)THEN
          IF(NITEMP.LE.1)THEN
            !hay que iterar la template al menos una vez mas; de lo contrario 
            !no podemos calcular DELTA_VEL ni DELTA_SIG
          ELSE
            DELTA_VEL=ABS(VELF-VELT)
            DELTA_SIG=ABS(SIGF-SIGT)
            IF((DELTA_VEL.GT.VEL_PRECISION).OR.
     +       (DELTA_SIG.GT.VEL_PRECISION))THEN
            ELSE
              !ya hemos alcanzado la precision requerida
              LREPITE_TEMPLATE=.FALSE.
            END IF
          END IF
        END IF
C------------------------------------------------------------------------------
        IF(LREPITE_TEMPLATE)THEN
          LOPTIMAL_TEMPLATE=.TRUE.
          NITEMP=NITEMP+1
          ITEMPL=(NITEMP+1)/2
          VELT=VELF
          SIGT=SIGF
          GOTO 444
        END IF
C------------------------------------------------------------------------------
C Si no hay que hacer simulaciones, mostramos resultados. En caso contrario, 
C vemos si hay que proseguir con las simulaciones o si, de haber terminado 
C éstas, mostramos resultados.
        IF(.NOT.LSIMULATIONS)THEN !.........................no hay simulaciones
          WRITE(*,*)
          WRITE(*,101) '#*******FINAL VALUES********'
          IF(LWRITE_LOGFILE)THEN
            WRITE(27,101) '#*******FINAL VALUES********'
          END IF
          IF(LPROBLEM_IS_A_FRAME)THEN
            IF(IWORKING_SCAN.EQ.0)THEN
              WRITE(*,101) '#.....CENTRAL SPECTRUM......'
              IF(LWRITE_LOGFILE)THEN
                WRITE(27,101) '#.....CENTRAL SPECTRUM......'
              END IF
            ELSE
              WRITE(*,'(A,I4,A)') '#.....SPECTRUM ',IWORKING_SCAN,
     +         ' ........'
              IF(LWRITE_LOGFILE)THEN
                WRITE(27,'(A,I4,A)') '#.....SPECTRUM ',IWORKING_SCAN,
     +           ' ........'
              END IF
            END IF
          END IF
          WRITE(*,*)
          WRITE(*,101) '#GAMMA     VEL       SIG'
          IF(LWRITE_LOGFILE)THEN
            WRITE(27,101) '#GAMMA     VEL       SIG'
          END IF
          GAMFINAL=GAMF
          VELFINAL=VELF
          SIGFINAL=SIGF
          WRITE(*,'(A1,F5.3,3X,F8.2,3X,F6.2,A)')
     +     '#',GAMFINAL,VELFINAL,SIGFINAL,'   (FINAL)'
          IF(LWRITE_LOGFILE)THEN
            WRITE(27,'(A1,F5.3,3X,F8.2,3X,F6.2,A)') 
     +       '#',GAMFINAL,VELFINAL,SIGFINAL,'   (FINAL)'
          END IF
        ELSE !.................................................hay simulaciones
          IF(ISIMUL.EQ.0)THEN
            !datos originales sin simular
            GAMFINAL=GAMF
            VELFINAL=VELF
            SIGFINAL=SIGF
          ELSE
            !resultado de las simulaciones
            GAMSIM(ISIMUL)=GAMF
            VELSIM(ISIMUL)=VELF
            SIGSIM(ISIMUL)=SIGF
          END IF
          IF(ISIMUL.LT.NSIMUL)THEN      !hay que continuar con las simulaciones
            ISIMUL=ISIMUL+1
            NITEMP=0
            IF(NTEMPLATES.GT.1) LOPTIMAL_TEMPLATE=.TRUE.
            GOTO 333
          END IF
          !ya hemos terminado con las simulaciones; hacemos estadística
          GAMF=0.
          VELF=0.
          SIGF=0.
          DO I=1,NSIMUL
            GAMF=GAMF+GAMSIM(I)
            VELF=VELF+VELSIM(I)
            SIGF=SIGF+SIGSIM(I)
          END DO
          GAMF=GAMF/REAL(NSIMUL)
          VELF=VELF/REAL(NSIMUL)
          SIGF=SIGF/REAL(NSIMUL)
          EGAMF=0.
          EVELF=0.
          ESIGF=0.
          DO I=1,NSIMUL
            EGAMF=EGAMF+(GAMSIM(I)-GAMFINAL)*(GAMSIM(I)-GAMFINAL)
            EVELF=EVELF+(VELSIM(I)-VELFINAL)*(VELSIM(I)-VELFINAL)
            ESIGF=ESIGF+(SIGSIM(I)-SIGFINAL)*(SIGSIM(I)-SIGFINAL)
          END DO
          EGAMF=SQRT(EGAMF/REAL(NSIMUL-1))
          EVELF=SQRT(EVELF/REAL(NSIMUL-1))
          ESIGF=SQRT(ESIGF/REAL(NSIMUL-1))
          WRITE(*,*)
          WRITE(*,101) '#*******FINAL VALUES********'
          IF(LWRITE_LOGFILE)THEN
            WRITE(27,101) '#*******FINAL VALUES********'
          END IF
          IF(LPROBLEM_IS_A_FRAME)THEN
            IF(IWORKING_SCAN.EQ.0)THEN
              WRITE(*,101) '#.....CENTRAL SPECTRUM......'
              IF(LWRITE_LOGFILE)THEN
                WRITE(27,101) '#.....CENTRAL SPECTRUM......'
              END IF
            ELSE
              WRITE(*,'(A,I4,A)') '#.....SPECTRUM ',IWORKING_SCAN,
     +         ' ........'
              IF(LWRITE_LOGFILE)THEN
                WRITE(27,'(A,I4,A)') '#.....SPECTRUM ',IWORKING_SCAN,
     +           ' ........'
              END IF
            END IF
          END IF
          WRITE(*,*)
          WRITE(*,101) '#GAMMA     VEL       SIG'
          IF(LWRITE_LOGFILE) WRITE(27,101) '#GAMMA     VEL       SIG'
          WRITE(*,'(A1,F5.3,3X,F8.2,3X,F6.2,A)')
     +     '#',GAMFINAL,VELFINAL,SIGFINAL,'   (FINAL)'
          IF(LWRITE_LOGFILE)THEN
            WRITE(27,'(A1,F5.3,3X,F8.2,3X,F6.2,A)')
     +       '#',GAMFINAL,VELFINAL,SIGFINAL,'   (FINAL)'
          END IF
          WRITE(*,'(A1,F5.3,3X,F8.2,3X,F6.2,A)') 
     +     '#',GAMF,VELF,SIGF,'   (MEANS)'
          IF(LWRITE_LOGFILE)THEN
            WRITE(27,'(A1,F5.3,3X,F8.2,3X,F6.2,A)') 
     +       '#',GAMF,VELF,SIGF,'   (MEANS)'
          END IF
          WRITE(*,'(A1,F5.3,3X,F8.2,3X,F6.2,A)') 
     +     '#',EGAMF,EVELF,ESIGF,'   (ERROR)'
          IF(LWRITE_LOGFILE)THEN
            WRITE(27,'(A1,F5.3,3X,F8.2,3X,F6.2,A)') 
     +       '#',EGAMF,EVELF,ESIGF,'   (ERROR)'
          END IF
        END IF
C------------------------------------------------------------------------------
C calculamos el r.m.s. en los residuos
        RESRMS=0.D0
        DO I=1,NCH
          RESDUM=DBLE(SGPLOT(I)-SMPLOT(I))
          RESRMS=RESRMS+RESDUM*RESDUM
        END DO
        RESRMS=SQRT(RESRMS/DBLE(NCH-1))
        WRITE(*,*)
        WRITE(*,'(A,E10.5)') '#RMS RES = ',REAL(RESRMS)
        IF(LWRITE_LOGFILE)THEN
          WRITE(27,'(A,E10.5)') '#RMS RES = ',REAL(RESRMS)
        END IF
C------------------------------------------------------------------------------
C Plots con el ajuste final para el scan considerado
        IF(CSHOW.EQ.'y')THEN
          !restamos un continuo solo para dibujar mejor (quitamos del ajuste 
          !las regiones enmascaradas)
          DO J=1,NCH
            IFCHAN(J)=.TRUE.
          END DO
          IF(NEML.GT.0)THEN
            DO IEML=1,NEML
              DO J=IEML1_NCH(IEML),IEML2_NCH(IEML)
                IFCHAN(J)=.FALSE.
              END DO
            END DO
          END IF
          J=0
          DO I=1,NCH
            IF(IFCHAN(I))THEN
              J=J+1
              XL(J)=X(I)
              YL(J)=SGPLOT(I)
              SIGMAY(J)=1.
            END IF
          END DO
          NPTS=J
          IF(NPTS.GE.NTERMS)THEN
            CALL POLFIT(XL,YL,SIGMAY,NPTS,NTERMS,0,CO,CHISQR)
            DO I=1,NCH
              YL(I)=FPOLY(NDEG,CO,X(I))
              SGPLOT(I)=SGPLOT(I)-YL(I)
              SMPLOT(I)=SMPLOT(I)-YL(I)
            END DO
          ELSE
            WRITE(*,100) 'NPTS, NTERMS: '
            WRITE(*,*) NPTS,NTERMS
            WRITE(*,101) 'FATAL ERROR: NPTS.LT.NTERMS'
            IF(LWRITE_LOGFILE) CLOSE(27)
            STOP
          END IF
          !calculamos residuos
          DO I=1,NCH
            RESID(I)=SGPLOT(I)-SMPLOT(I)
          END DO
C..............................................................................
          CALL FINDMML(NCH,1,NCH,SMPLOT,YMIN,YMAX)
          CALL FINDMMLMASK(NCH,1,NCH,SGPLOT,IFCHAN,YMING,YMAXG)
          YMIN=AMIN1(YMIN,YMING)
          YMAX=AMAX1(YMAX,YMAXG)
          DYG=(YMAX-YMIN)*0.05
          YMIN=YMIN-DYG
          YMAX=YMAX+DYG
          !los limites en X los calculamos sin expansion para que los
          !ticks en longitud de onda sean correctos
          XMIN=1.0       !sin expansion
          XMAX=REAL(NCH) !sin expansion
ccc       XMIN=-1.*REAL(NCH)/100.
ccc       XMAX=REAL(NCH)*1.01
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGSUBP(1,1)
            CALL PGPAGE
C dibujamos espectros PROBLEM y MODEL
            CALL AUTOPLOT(NCH,X,SGPLOT,1,NCH, !.........................PROBLEM
     +       ' ',' ',' ',
     +       .FALSE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +       0,.TRUE.,'BC','BCNTS',
     +       101,2,
     +       0.,1.,0.0,0.60)
            CALL AUTOPLOT(NCH,X,SMPLOT,1,NCH, !...........................MODEL
     +       ' ',' ',' ',
     +       .FALSE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +       0,.FALSE.,'BC','BCNTS',
     +       101,3,
     +       0.,1.,0.0,0.60)
            CALL PGMTXT('B',2.5,0.5,0.5,'wavelength (template frame)')
            CALL PGMTXT('L',2.5,0.5,0.5,'rectified flux')
C ticks en longitud de onda
            CALL RETORNA_TICKS(IDN(ITERM),STWV+REAL(NC1-1)*DISP,
     +       STWV+REAL(NC2-1)*DISP,NMAX_TICKS,NTICKS,XTICK,
     +       NTICKS_SMALL,XTICK_SMALL)
            !ticks grandes
            DO I=1,NTICKS
              WLIN=XTICK(I)
              WLIN=(LOG(WLIN)-LOG(WI0))/REAL(DLOGW)+0.5
!             WLIN=WLIN-REAL(IVEL0PLOT)
              CALL PGMOVE(WLIN,YMIN)
              CALL PGDRAW(WLIN,YMIN+(YMAX-YMIN)/25.)
              CALL PGSCH(0.85)
              IF(WLIN.GT.0..AND.WLIN.LT.REAL(NCH))THEN
                MM=INT(XTICK(I)*1000) !le pongo 3 decimales
                PP=-3
                CALL PGNUMB(MM,PP,0,CDUMMY,L)
                CALL PGPTEXT(WLIN,YMIN-(YMAX-YMIN)/13.,0.,0.5,
     +           CDUMMY(1:L))
              END IF
              CALL PGSCH(1.2)
            END DO
            !ticks peque~nos
            DO I=1,NTICKS_SMALL
              WLIN=XTICK_SMALL(I)
              WLIN=(LOG(WLIN)-LOG(WI0))/REAL(DLOGW)+0.5
!             WLIN=WLIN-REAL(IVEL0PLOT)
              CALL PGMOVE(WLIN,YMIN)
              CALL PGDRAW(WLIN,YMIN+(YMAX-YMIN)/50.)
            END DO
C Residuos (rectificados)
            CALL AUTOPLOT(NCH,X,RESID,1,NCH, !.........................RESIDUOS
     +       ' ',' ',' ',
     +       .FALSE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +       0,.FALSE.,'BC','BCNTS',
     +       101,1,
     +       0.,1.,0.40,1.00)
            CALL PGSCH(0.85)
            CALL PGBOX('BCNTS',0.0,0,' ',0.0,0)
            CALL PGSCH(1.2)
            CALL PGMTXT('L',2.5,0.5,0.5,'rectified residuals')
            DO I=1,NEML !.................................regiones a enmascarar
              WLIN=STWV+REAL(IEML1(I)-1)*DISP
              WLIN=(LOG(WLIN)-LOG(WI0))/REAL(DLOGW)+0.5
              LIN1=INT(WLIN+0.5)-IVEL0PLOT
              WLIN=STWV+REAL(IEML2(I)-1)*DISP
              WLIN=(LOG(WLIN)-LOG(WI0))/REAL(DLOGW)+0.5
              LIN2=INT(WLIN+0.5)-IVEL0PLOT
              CALL AUTOPLOT(LIN2-LIN1+1,X,RESID,LIN1,LIN2,
     +         ' ',' ',' ',
     +         .FALSE.,XMIN,XMAX,YMIN,YMAX,0.05,
     +         0,.FALSE.,' ',' ',
     +         101,5,
     +         0.,1.,0.40,1.00)
            END DO
            CALL PGSCI(3)
            CALL PGMTEXT('T',2.0,0.,0.,'MODEL')
            CALL PGSCI(2)
            CALL PGMTEXT('T',3.0,0.,0.,'PROBLEM: '//INFILE2)
            CALL PGSCI(1)
            CALL PGMTEXT('T',3.0,1.0,1.0,'FINAL FIT')
            WRITE(CINFO,'(3(A,I8))')
     +       'ISCAN=',IWORKING_SCAN,',_IITER=',IITER,
     +       ',_ITEMPL=',ITEMPL
            CALL RMBLANK(CINFO,CINFO,LINFO)
            LOOP=.TRUE. !eliminamos pseudo-espacios en blanco
            DO WHILE(LOOP)
              LFOUND=INDEX(CINFO(1:LINFO),'_')
              IF(LFOUND.EQ.0)THEN
                LOOP=.FALSE.
              ELSE
                CINFO(LFOUND:LFOUND)=' '
              END IF
            END DO
            CALL PGSCH(0.85)
            CALL PGMTEXT('T',2.82,1.0,1.0,CINFO(1:LINFO))
            CALL PGSCH(1.20)
C Pinta posicion de las lineas de emision
            ZV=1.+VELFINAL/C
            CALL PGSLS(4)
            CALL PGSCI(4)
            DYG=(YMAX-YMIN)/20.0
            IF((CPLOT_EML.EQ.'1').OR.(CPLOT_EML.EQ.'3'))THEN
              !Lineas de Balmer
              CALL SELLINES(0,NLINES,WLINES,LABLINES)
              LPARIMPAR=.TRUE.
              DO J=1,NLINES
                WLIN=WLINES(J)*ZV
                WLIN=(LOG(WLIN)-LOG(WI0))/REAL(DLOGW)+0.5
                WLIN=WLIN-REAL(IVEL0PLOT)
                CALL PGMOVE(WLIN,YMIN)
                CALL PGDRAW(WLIN,YMAX)
                CALL PGSCH(0.5)
                IF(WLIN.GT.0..AND.WLIN.LT.REAL(NCH))THEN
                  IF(LPARIMPAR)THEN
                    LPARIMPAR=.FALSE.
                    CALL PGPTEXT(WLIN,YMAX+0.5*DYG,0.,0.5,LABLINES(J))
                  ELSE
                    LPARIMPAR=.TRUE.
                    CALL PGPTEXT(WLIN,YMAX+1.5*DYG,0.,0.5,LABLINES(J))
                  END IF
                END IF
                CALL PGSCH(1.2)
              END DO
              CALL PGSCI(8)
              !Lineas de emision tipicas
              CALL SELLINES(1,NLINES,WLINES,LABLINES)
              LPARIMPAR=.TRUE.
              DO J=1,NLINES
                WLIN=WLINES(J)*ZV
                WLIN=(LOG(WLIN)-LOG(WI0))/REAL(DLOGW)+0.5
                WLIN=WLIN-REAL(IVEL0PLOT)
                CALL PGMOVE(WLIN,YMIN)
                CALL PGDRAW(WLIN,YMAX)
                CALL PGSCH(0.5)
                IF(WLIN.GT.0..AND.WLIN.LT.REAL(NCH))THEN
                  IF(LPARIMPAR)THEN
                    LPARIMPAR=.FALSE.
                    CALL PGPTEXT(WLIN,YMAX+0.5*DYG,0.,0.5,LABLINES(J))
                  ELSE
                    LPARIMPAR=.TRUE.
                    CALL PGPTEXT(WLIN,YMAX+1.5*DYG,0.,0.5,LABLINES(J))
                  END IF
                END IF
                CALL PGSCH(1.2)
              END DO
            END IF
            IF((CPLOT_EML.EQ.'2').OR.(CPLOT_EML.EQ.'3'))THEN
              !Lineas de cielo
              CALL PGSCI(7)
              CALL SELLINES(2,NLINES,WLINES,LABLINES)
              DO J=1,NLINES
                WLIN=WLINES(J)*ZV
                WLIN=(LOG(WLIN)-LOG(WI0))/REAL(DLOGW)+0.5
                WLIN=WLIN-REAL(IVEL0PLOT)
                CALL PGMOVE(WLIN,YMIN)
                CALL PGDRAW(WLIN,YMAX)
              END DO
            END IF
            CALL PGSLS(1)
            CALL PGSCI(1)
            CALL PGSCH(1.2)
            CALL PGIDEN_RED
          END DO
          IF(LPAUSE)THEN
            WRITE(*,100) 'Press <CR> to continue...'
            READ(*,*)
          END IF
C..............................................................................
        END IF
C------------------------------------------------------------------------------
C Bloque para controlar la realización de las medidas en una imagen PROBLEM con
C varios scans
        IF(LPROBLEM_IS_A_FRAME)THEN
          IF(IWORKING_SCAN.EQ.0)THEN !.............................scan central
            GALGAM(NSCENTER)=GAMFINAL
            GALVEL(NSCENTER)=VELFINAL
            GALSIG(NSCENTER)=SIGFINAL
            IF(LSIMULATIONS)THEN
              EGALGAM(NSCENTER)=EGAMF
              EGALVEL(NSCENTER)=EVELF
              EGALSIG(NSCENTER)=ESIGF
            END IF
            GAMCEN=GAMFINAL
            VELCEN=VELFINAL
            SIGCEN=SIGFINAL
            !después de medir el scan central, vamos a empezar a medir en el
            !resto de la imagen, recorriendo los scans [NSCMIN,NSCMAX]
            IWORKING_SCAN=NSCMIN-1
          ELSE !..........................................no es el scan central
            GALGAM(IWORKING_SCAN)=GAMFINAL
            GALVEL(IWORKING_SCAN)=VELFINAL
            GALSIG(IWORKING_SCAN)=SIGFINAL
            IF(LSIMULATIONS)THEN
              EGALGAM(IWORKING_SCAN)=EGAMF
              EGALVEL(IWORKING_SCAN)=EVELF
              EGALSIG(IWORKING_SCAN)=ESIGF
            END IF
          END IF
          IWORKING_SCAN=IWORKING_SCAN+1
          IF(.NOT.LPLOT_ALL) CSHOW='n'
          WRITE(*,*)
          IF(IWORKING_SCAN.LE.NSCMAX) GOTO 1000
        END IF
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C Codigo para escribir resultados
        LOOP=.TRUE.
        DO WHILE(LOOP)
          WRITE(*,100)'Output file name (RET if no output)? '
          READ(*,101) OUTFILE
          IF(TRUELEN(OUTFILE).EQ.0)THEN
            LOOP=.FALSE.
          ELSE
            INQUIRE(FILE=OUTFILE,EXIST=LOUT)
            IF(LOUT)THEN
              WRITE(*,101)'ERROR: this file already exist. Try again.'
            ELSE
              LOOP=.FALSE.
            END IF
          END IF
        END DO
        IF(TRUELEN(OUTFILE).GT.0)THEN
          OPEN(22,FILE=OUTFILE,STATUS='NEW',FORM='FORMATTED')
          WRITE(22,101) 'FILENAME: '//INFILE2(1:60)
          WRITE(22,101) 'TEMPLATE: '//INFILE1(1:60)
          IF(NTEMPLATES.GT.1)THEN
            WRITE(22,101) 'OPTIMAL TEMPLATE:'
            DO J=1,NTEMPLATES
              WRITE(22,'(I3,A1,F5.2)') J,':',WEI(J)
            END DO
          END IF
          IF(LSIMULATIONS)THEN
            WRITE(22,'(I4,A)') NSIMUL,' SIMULATIONS'
          ELSE
            WRITE(22,101) '  NO SIMULATIONS'
          END IF
          WRITE(22,101) '----------------------------------------------'
          WRITE(22,'(13X,A5,5X,A3,6X,A3,5X,A3)')
     +     'SCANS','GAM','VEL','SIG'
          IF(LPROBLEM_IS_A_FRAME)THEN
            WRITE(*,100)'Header information file'
            FILEHEAD=INFILEX(23,INFILE2(1:TRUELEN(INFILE2))
     +       //'.log',0,0,.0,.0,3,.FALSE.)
            DO I=1,NSCAN2
              READ(23,'(A15)') CSCANBINNING(I)
            END DO
            CLOSE(23)
            DO I=NSCMIN,NSCMAX
              WRITE(22,'(I4,A15,3X,F5.3,3X,F7.1,3X,F5.1)')
     +         I,CSCANBINNING(I),GALGAM(I),GALVEL(I),GALSIG(I)
              IF(LSIMULATIONS)THEN
                WRITE(22,'(22X,F5.3,3X,F7.1,3X,F5.1)')
     +           EGALGAM(I),EGALVEL(I),EGALSIG(I)
              END IF
            END DO
          ELSE
            WRITE(22,'(12X,I4,6X,F5.3,3X,F7.1,3X,F5.1)') 
     +       NS0,GAMFINAL,VELFINAL,SIGFINAL
            IF(LSIMULATIONS)THEN
              WRITE(22,'(22X,F5.3,3X,F7.1,3X,F5.1)')EGAMF,EVELF,ESIGF
            END IF
          END IF
          WRITE(22,101) '----------------------------------------------'
          CLOSE(22)
        END IF
C------------------------------------------------------------------------------
C Codigo para salvar la template optima. Varias opciones
        WRITE(*,100)'Save optimal template (y/n) '
        CSAVE(1:1)=READC('n','yn')
        IF(CSAVE.EQ.'y')THEN
          LOOP=.TRUE.
          DO WHILE(LOOP)
            WRITE(*,*)
            IF(NTEMPLATES.GT.1) WRITE(*,101)'  1: Composite template'
            WRITE(*,101)'  2: Composite template, corrected for '//
     +       'gamma effect'
            WRITE(*,101)'  3: Final fit (templ+gal) in log scale'
            WRITE(*,101)'  0: Exit'
            WRITE(*,100) 'Choose option '
            NSAVE=READILIM('@',0,3)
C..............................................................................
            IF(NSAVE.EQ.0)THEN
              LOOP=.FALSE.
C..............................................................................
            ELSEIF(NSAVE.EQ.1)THEN
              IF(NTEMPLATES.GT.1)THEN
                OBJECT='OPT.TEMPL.'
                WRITE(*,100)'Output file name'
                OUTFILE=OUTFILEX(16,'@',1,NCHAN,STWV,DISP,1,.FALSE.)
                WRITE(16) (SORIG1(J),J=1,NCHAN)
                CLOSE(16)
                IF(LSIMULATIONS)THEN
                  DO J=1,NCHAN
                    EORIG1(J)=0.
                    DO II=1,NTEMPLATES
                      EORIG1(J)=EORIG1(J)+
     +                 WEI(II)*WEI(II)*ETEMP(J,II)*ETEMP(J,II)
                    END DO
                    EORIG1(I)=SQRT(EORIG1(I))
                  END DO
                  WRITE(*,100)'Output error file name '
                  CALL GUESSEF(OUTFILE,ERRFILE)
                  OUTFILE=OUTFILEX(17,ERRFILE,1,NCHAN,STWV,DISP,
     +             1,.TRUE.)
                  WRITE(17) (EORIG1(J),J=1,NCHAN)
                  CLOSE(17)
                END IF
              END IF
C..............................................................................
c Aplicamos el factor gamma: para ello dividimos por un continuo, restamos la 
C media, aplicamos el factor gamma, sumamos la media y multiplicamos por
C el continuo. El ajuste del continuo se realiza entre NC1 y NC2.
            ELSEIF(NSAVE.EQ.2)THEN
              OBJECT='OPT.TEMPL'
              J=0
              DO I=NC1,NC2
                J=J+1
                XL(J)=X(I)
                YL(J)=SORIG1(I)
                SIGMAY(J)=1.
              END DO
              NPTS=J
              IF(CCONT.EQ.'1')THEN
                CALL POLFIT(XL,YL,SIGMAY,NPTS,NTERMS,0,CO,CHISQR)
              ELSE
                CALL PSEUDOFIT(XL,YL,NPTS,NTERMS,YRMSTOL,
     +           PSEUDO_WEIGHT,PSEUDO_POWER,.TRUE.,CO)
              END IF
              DO I=1,NCHAN
                YL(I)=FPOLY(NDEG,CO,X(I))
                SORIGX(I)=SORIG1(I)/YL(I)
              END DO
              MEAN=0.D0
              DO I=NC1,NC2
                MEAN=MEAN+DBLE(SORIGX(I))
              END DO
              MEAN=MEAN/DBLE(NC2-NC1+1)
              DO I=1,NCHAN
                SORIGX(I)=SORIGX(I)-REAL(MEAN)
              END DO
              IF(LPROBLEM_IS_A_FRAME) GAMFINAL=GAMCEN
              DO I=1,NCHAN
                SORIGX(I)=SORIGX(I)*GAMFINAL
                SORIGX(I)=SORIGX(I)+REAL(MEAN)
                SORIGX(I)=SORIGX(I)*YL(I)
              END DO
              WRITE(*,100)'Output file name'
              OUTFILE=OUTFILEX(16,'@',1,NCHAN,STWV,DISP,1,.FALSE.)
              WRITE(16) (SORIGX(J),J=1,NCHAN)
              CLOSE(16)
C..............................................................................
            ELSEIF(NSAVE.EQ.3)THEN
              IF(CSHOW.NE.'y')THEN
                DO I=1,NCH
                  RESID(I)=SGPLOT(I)-SMPLOT(I)
                END DO
              END IF
              OBJECT='GAL./MODEL/RES.'
              COMMENT='log wavelength scale'
              WRITE(*,100)'Output file name'
              OUTFILE=OUTFILEX(16,'@',3,NCH,0.,0.,1,.FALSE.)
              WRITE(*,100) 'Writing galaxy, model and residuals...'
              WRITE(16) (SGPLOT(J),J=1,NCH)
              WRITE(16) (SMPLOT(J),J=1,NCH)
              WRITE(16) (RESID(J),J=1,NCH)
              CLOSE(16)
              WRITE(*,101) '  OK!'
            END IF
C..............................................................................
          END DO
        END IF
C------------------------------------------------------------------------------
        CALL PGEND
        IF(LWRITE_LOGFILE) CLOSE(27)
        STOP
100     FORMAT(A,$)
101     FORMAT(A)
        END
C
C******************************************************************************
C******************************************************************************
C
        REAL FUNCTION CHISQG(XX)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        REAL XX(3)
C
        REAL PI
        PARAMETER (PI=3.141593)

        INTEGER I,NCH
        INTEGER KMIN,KMAX
        INTEGER NTEMPLATES,NP,NF1,NF2
        REAL STR(NMAXFFT),STI(NMAXFFT)
        REAL SR(NMAXFFT),SI(NMAXFFT),WIEN(NMAXFFT)
        REAL SIG,VEL,GAM,SIGMIN,GAMMIN
        DOUBLE PRECISION DCHISQG,ZQ,ZQ1,ZQ3,X1,X2,X3
        DOUBLE PRECISION ZCOS,ZSIN,ZEXP,A,B,FR,FI
        DOUBLE PRECISION DLOGW,DV

        COMMON/BLKFIT/NCH,KMIN,KMAX,STR,STI,SR,SI,WIEN
        COMMON/BLKTEMPL1/NTEMPLATES,NCHAN,NP,NF1,NF2,SIG,VEL,GAM
        COMMON/BLKTEMPL1B/DLOGW,DV
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C Esto es para que no se obtengan valores negativos o muy pequeños de gamma
        GAMMIN=0.001
        IF(XX(1).LT.GAMMIN)THEN
          CHISQG=1.E25
          RETURN
        END IF
C Esto es para que no se obtengan valores negativos o muy pequeños de sigma
        SIGMIN=0.000001
        IF(XX(3).LT.SIGMIN)THEN
          CHISQG=1.E25
          RETURN
        END IF
C
        X1=DBLE(XX(1))
        X2=DBLE(XX(2))
        X3=DBLE(XX(3))
C
        DCHISQG=0.D0
        ZQ=DBLE(2*PI/REAL(NCH))
        ZQ3=-0.5D0
        DO I=KMIN,KMAX
          ZQ1=ZQ*DBLE(I-1)
          ZCOS=DCOS(ZQ1*X2)
          ZSIN=DSIN(ZQ1*X2)
          ZEXP=X1*DEXP(ZQ3*ZQ1*ZQ1*X3*X3)
          A=ZCOS*ZEXP
          B=ZSIN*ZEXP
          FR=DBLE(SR(I))-DBLE(STR(I))*A+DBLE(STI(I))*B
          FI=DBLE(SI(I))-DBLE(STR(I))*B-DBLE(STI(I))*A
          DCHISQG=DCHISQG+DBLE(WIEN(I))*(FR*FR+FI*FI)
        END DO
        CHISQG=REAL(DCHISQG)
C
        END
C
C******************************************************************************
C
        SUBROUTINE FFTBROAD(NCH,XX)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        REAL XX(3)
C
        REAL PI
        PARAMETER (PI=3.141593)
C
        INTEGER I,NCH
        REAL FFTA(NMAXFFT),FFTB(NMAXFFT)
        DOUBLE PRECISION ZQ,ZQ1,ZQ3,X1,X2,X3
        DOUBLE PRECISION ZCOS,ZSIN,ZEXP
C
        COMMON/BLKFFTB/FFTA,FFTB
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        X1=DBLE(XX(1))
        X2=DBLE(XX(2))
        X3=DBLE(XX(3))
C
        ZQ=DBLE(2*PI/REAL(NCH))
        ZQ3=-0.5D0
        DO I=1,NCH/2
          ZQ1=ZQ*DBLE(I-1)
          ZCOS=DCOS(ZQ1*X2)
          ZSIN=DSIN(ZQ1*X2)
          ZEXP=X1*DEXP(ZQ3*ZQ1*ZQ1*X3*X3)
          FFTA(I)=REAL(ZCOS*ZEXP)
          FFTB(I)=REAL(ZSIN*ZEXP)
        END DO
C
        END
C
C******************************************************************************
C Transforma a escala logaritmica
        SUBROUTINE LOGSCALE(NCHAN,NP,DLOGW,W,SI,SF)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
C
        INTEGER MJ,I,NP
        INTEGER LMIN,LMIN2,LMIN3
        REAL SI(NCMAX),SF(NMAXFFT)
        DOUBLE PRECISION W(NCMAX+1)
        DOUBLE PRECISION S2(NMAXFFT)     !espectro de salida (esc. logaritmica)
        DOUBLE PRECISION S3(NCMAX)           !espectro de entrada (esc. lineal)
        DOUBLE PRECISION DLOGW,WL1,WL2,EXC,RES
        DOUBLE PRECISION DL1,DL2
!       DOUBLE PRECISION XCNT1,XCNT2
!       DOUBLE PRECISION DIFF
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C protecciones
        IF(NCHAN.GT.NCMAX)THEN
          WRITE(*,100) 'NCHAN, NCMAX: '
          WRITE(*,*) NCHAN,NCMAX
          STOP 'FATAL ERROR in LOGSCALE: NCHAN.GT.NCMAX.'
        END IF
        IF(NP.GT.NMAXFFT)THEN
          WRITE(*,100) 'NP, NMAXFFT: '
          WRITE(*,*) NP,NMAXFFT
          STOP 'FATAL ERROR in LOGSCALE: NP.GT.NMAXFFT.'
        END IF
C
        DO MJ=1,NCHAN
          S3(MJ)=DBLE(SI(MJ))
        END DO
C
        DO MJ=1,NP
          S2(MJ)=0.0D+00
        END DO
C
C BIN TO LOG WAVELENGTH
C
C RES es el numero de cuentas del pixel anterior (lineal) no usadas en la
C generacion del ultimo pixel (logaritmico)
        RES=0.0D+00
        WL1=W(1)                           !LOG(l.d.o. inicial del nuevo pixel)
        WL2=WL1+DLOGW                        !LOG(l.d.o. final del nuevo pixel)
        LMIN=1
        LMIN2=2
        IF(WL2.GE.W(LMIN2)) RES=S3(1)
        DO 12 MJ=1,NP
          LMIN3=LMIN2+1
C si el nuevo pixel (logaritmico) contiene mas de un antiguo pixel (lineal)
C se van sumando. Se suma hasta el ultimo pixel (lineal) completo abarcado,
C es decir, siempre se suma por defecto
          DO 13 I=LMIN3,NCHAN
            LMIN2=I-1
            IF(W(I).GT.WL2) GOTO 14
            S2(MJ)=S3(LMIN2)+S2(MJ)
 13       CONTINUE
 14       CONTINUE
          S2(MJ)=S2(MJ)+RES
          IF(W(LMIN2).GE.WL2) GOTO 15
          DL1=W(LMIN2+1)-W(LMIN2)
          DL2=WL2-W(LMIN2)
          EXC=S3(LMIN2)*(DL2/DL1)
          RES=S3(LMIN2)-EXC
          S2(MJ)=S2(MJ)+EXC
          WL1=WL2
          WL2=WL2+DLOGW
          IF(WL2.LE.W(LMIN2+1)) RES=0.0D+00
          LMIN=LMIN2
          LMIN2=LMIN+1
          GOTO 20
 15       DL1=W(LMIN2)-W(LMIN)
          DL2=WL2-WL1
          S2(MJ)=S2(MJ)+S3(LMIN)*(DL2/DL1)
          WL1=WL2
          WL2=WL2+DLOGW
          IF(WL2.GT.W(LMIN2)) GOTO 19
          RES=0.0D+00
          GOTO 20
 19       DL1=W(LMIN2)-W(LMIN)
          DL2=W(LMIN2)-WL1
          RES=S3(LMIN)*(DL2/DL1)
 20       CONTINUE
 12     CONTINUE
C
C ADD UP COUNTS TO CHECK
!
!        XCNT1=0.0D0
!        XCNT2=0.0D0
!        DO I=1,NCHAN
!          XCNT1=XCNT1+S3(I)
!        END DO
!        DO I=1,NP
!          XCNT2=XCNT2+S2(I)
!        END DO
!        WRITE(*,100)'Check number of counts: '
!        DIFF=DABS(XCNT2-XCNT1)
!        WRITE(*,*) XCNT1,XCNT2,DIFF
C
        DO I=1,NP
          SF(I)=REAL(S2(I))
        END DO
C
100     FORMAT(A,$)
        END
C
C******************************************************************************
C
        SUBROUTINE FITTEMPL
C
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
C
        REAL PI,C
        PARAMETER (PI=3.141593)
        PARAMETER (C=299792.46)
        INTEGER NMAX_TEMPLATES
        PARAMETER (NMAX_TEMPLATES=20)
        INTEGER NMAX_EML
        PARAMETER (NMAX_EML=100)
C
        REAL FINTGAUSS
        REAL FPOLY
        EXTERNAL CHISQT
        REAL CHISQT
C
        INTEGER I,II,J,JJ,I0,J0
        INTEGER NTEMPLATES,NP,NF1,NF2,NPTS,NCH,NCHF
        INTEGER NVAR,NEVAL,NTEM
        INTEGER IEML,NEML
        INTEGER IEML1_SGL(NMAX_EML),IEML2_SGL(NMAX_EML)
        INTEGER NTERMS,NDEG
        REAL SPG(NCMAX),TEMPX(NCMAX,NMAX_TEMPLATES)
        REAL STL(NMAXFFT),SML(NMAXFFT),WEI(NMAX_TEMPLATES)
        REAL DWEI(NMAX_TEMPLATES)
        REAL SPT(NCMAX),TEMPL(NMAXFFT,NMAX_TEMPLATES),GAL(NMAXFFT)
        REAL SIG,VEL,GAM,PCEN,FAC,ARG,SIGC,VELO,FSCALE
        REAL X(NMAXFFT),XL(NMAXFFT),YL(NMAXFFT)
        REAL SIGMAY(NMAXFFT),CO(20),CHISQR
        REAL X0(NMAX_TEMPLATES),DX0(NMAX_TEMPLATES)
        REAL XX(NMAX_TEMPLATES),DXX(NMAX_TEMPLATES),FSUM
        REAL YRMSTOL,PSEUDO_WEIGHT,PSEUDO_POWER
        REAL BETA,RCVEL1
        DOUBLE PRECISION DLOGW,DV,W(NCMAX+1),MEAN
        CHARACTER*1 CCONT
        LOGICAL IFCHAN(NMAXFFT),IFCHANCH(NMAXFFT)
C
        COMMON/BLK_X_GLOBAL/X
        COMMON/BLKTEMPL1/NTEMPLATES,NCHAN,NP,NF1,NF2,SIG,VEL,GAM
        COMMON/BLKTEMPL1B/DLOGW,DV
        COMMON/BLKTEMPL2/W,SPG,TEMPX,WEI
        COMMON/BLKTEMPL3/NEML,IEML1_SGL,IEML2_SGL
        COMMON/BLKTEMPL4/NTERMS,NDEG,CCONT
        COMMON/BLKTEMPL5/YRMSTOL
        COMMON/BLKTEMPL6/PSEUDO_WEIGHT,PSEUDO_POWER
        COMMON/BLKTEMPFIT/NCHF,NTEM,GAL,TEMPL
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        WRITE(*,*)
        WRITE(*,100)'-|-|-> Computing optimal template (please wait...)'
C
        NCH=NF2-NF1+1
        SIGC=SIG/REAL(DV)
!       VELO=VEL/REAL(DV)                                                  !mal
        BETA=VEL/C
        RCVEL1=(1.0+BETA)/SQRT(1.0-BETA*BETA)            !correcion relativista
        VELO=ALOG(RCVEL1)/DLOGW
        FSCALE=1./SQRT(2*PI)/SIGC
        ARG=-1./2./SIGC/SIGC
        I0=INT(10.*SIGC+0.5)
C------------------------------------------------------------------------------
C definimos región en la que vamos a buscar la template optima (solo entre 
C NF1 y NF2, eliminando las regiones a enmascarar); IFCHAN se utiliza
C para definir region util en la escala que va de 1 a NP, y IFCHANCH se
C emplea para definir la misma region en la escala que va de 1 a NCH
C
C IFCHAN: escala de 1 a NP
        IF(NF1.GT.1)THEN
          DO J=1,NF1-1
            IFCHAN(J)=.FALSE.
          END DO
        END IF
        DO J=NF1,NF2
          IFCHAN(J)=.TRUE.
        END DO
        IF(NF2.LT.NP)THEN
          DO J=NF2+1,NP
            IFCHAN(J)=.FALSE.
          END DO
        END IF
        !regiones a enmascarar
        IF(NEML.GT.0)THEN
          DO IEML=1,NEML
            DO J=IEML1_SGL(IEML),IEML2_SGL(IEML)
              IFCHAN(J)=.FALSE.
            END DO
          END DO
        END IF
C IFCHANCH: escala de 1 a NCH
        DO J=1,NCH
          IFCHANCH(J)=.TRUE.
        END DO
        !regiones a enmascarar
        IF(NEML.GT.0)THEN
          DO IEML=1,NEML
            DO J=IEML1_SGL(IEML),IEML2_SGL(IEML)
              JJ=J-NF1+1
              IF((JJ.GE.1).AND.(JJ.LE.NCH))THEN
                IFCHANCH(JJ)=.FALSE.
              END IF
            END DO
          END DO
        END IF
C Calculamos el numero de puntos efectivo para el calculo del minimo
        NCHF=0
        DO J=1,NCH
          IF(IFCHANCH(J)) NCHF=NCHF+1
        END DO
        IF(NCHF.LT.NTERMS)THEN
          WRITE(*,100) 'NCHF, NTERMS: '
          WRITE(*,*) NCHF,NTERMS
          WRITE(*,101) 'FATAL ERROR: NCHF.LT.NTERMS in FITTEMPL'
          STOP
        END IF
C------------------------------------------------------------------------------
C Convertimos a escala logaritmica, ensanchamos y movemos en velocidad radial
C las diferentes templates
        DO II=1,NTEMPLATES
          DO I=1,NCHAN
            SPT(I)=TEMPX(I,II)
          END DO
          CALL LOGSCALE(NCHAN,NP,DLOGW,W,SPT,STL)           !escala logaritmica
          IF(I0.EQ.0)THEN !no es necesario hacer la integral numéricamente
            DO I=1,NP
              PCEN=REAL(I)-VELO
              J0=INT(PCEN+0.5)
              IF((J0.GE.1).AND.(J0.LE.NP))THEN
                SML(I)=STL(J0)
              ELSE
                SML(I)=0.0
              END IF
            END DO
          ELSE !integramos numéricamente
            DO I=1,NP
              SML(I)=0.
              PCEN=REAL(I)-VELO                    !movemos en velocidad radial
              DO J=-I0,I0                                          !ensanchamos
                J0=INT(PCEN+0.5)+J
                IF((J0.GE.1).AND.(J0.LE.NP))THEN
                  FAC=FSCALE*
     +             FINTGAUSS(REAL(J0)-0.5,REAL(J0)+0.5,20,PCEN,ARG)
                  SML(I)=SML(I)+STL(J0)*FAC
                END IF
              END DO
            END DO
          END IF
C Se ajusta un polinomio o pseudocontinuo
          J=0
          DO I=1,NP
            IF(IFCHAN(I))THEN
              J=J+1
              XL(J)=X(I)
              YL(J)=SML(I)
              SIGMAY(J)=1.
            END IF
          END DO
          NPTS=J
          IF(CCONT.EQ.'1')THEN
            CALL POLFIT(XL,YL,SIGMAY,NPTS,NTERMS,0,CO,CHISQR)
          ELSE
            CALL PSEUDOFIT(XL,YL,NPTS,NTERMS,YRMSTOL,
     +       PSEUDO_WEIGHT,PSEUDO_POWER,.TRUE.,CO)
          END IF
C Normalizamos por el ajuste anterior
          DO I=1,NP
            YL(I)=FPOLY(NDEG,CO,X(I))
            SML(I)=SML(I)/YL(I)
          END DO
C Se extraen los pixels utiles y se resta la media
          MEAN=0.D0
          J=0
          DO I=1,NCH
            IF(IFCHANCH(I))THEN
              J=J+1
              TEMPL(J,II)=SML(I+NF1-1)
              MEAN=MEAN+DBLE(TEMPL(J,II))
            END IF
          END DO
          MEAN=MEAN/DBLE(NCHF) !ojo, es NCHF y no NCH
          DO J=1,NCHF
            TEMPL(J,II)=TEMPL(J,II)-REAL(MEAN)
          END DO
C y se introduce gamma
          DO J=1,NCHF
            TEMPL(J,II)=TEMPL(J,II)*GAM
          END DO
        END DO
C------------------------------------------------------------------------------
C Se hace lo mismo para la galaxia
        CALL LOGSCALE(NCHAN,NP,DLOGW,W,SPG,STL)             !escala logaritmica
        J=0
        DO I=1,NP
          IF(IFCHAN(I))THEN
            J=J+1
            XL(J)=X(I)
            YL(J)=STL(I)
            SIGMAY(J)=1.
          END IF
        END DO
        NPTS=J
        IF(CCONT.EQ.'1')THEN
          CALL POLFIT(XL,YL,SIGMAY,NPTS,NTERMS,0,CO,CHISQR)
        ELSE
          CALL PSEUDOFIT(XL,YL,NPTS,NTERMS,YRMSTOL,
     +     PSEUDO_WEIGHT,PSEUDO_POWER,.TRUE.,CO)
        END IF
        DO I=1,NP
          YL(I)=FPOLY(NDEG,CO,X(I))
          SML(I)=STL(I)/YL(I)
        END DO
        MEAN=0.D0
        J=0
        DO I=1,NCH
          IF(IFCHANCH(I))THEN
            J=J+1
            GAL(J)=SML(I+NF1-1)
            MEAN=MEAN+DBLE(GAL(J))
          END IF
        END DO
        MEAN=MEAN/DBLE(NCHF) !ojo, es NCHF y no NCH
        DO J=1,NCHF
          GAL(J)=GAL(J)-REAL(MEAN)
        END DO
C------------------------------------------------------------------------------
C Buscamos la template optima
        NTEM=NTEMPLATES
        NVAR=NTEMPLATES
        DO I=1,NVAR
          DX0(I)=1./REAL(NTEMPLATES)
          X0(I)=1./REAL(NTEMPLATES)
        END DO
C
        CALL DOWNHILL(NVAR,X0,DX0,CHISQT,1.0,0.5,2.0,YRMSTOL,
     +   XX,DXX,NEVAL)
C
        FSUM=0.
        DO J=1,NTEMPLATES
          FSUM=FSUM+XX(J)
        END DO
C
        WRITE(*,101) '  OK'
        WRITE(*,*)
!Mejor mostramos los pesos en el listado junto con las medidas
!       WRITE(*,101) 'WEIGHTS FOR THE TEMPLATES'
        DO J=1,NTEMPLATES
          WEI(J)=XX(J)/FSUM
          DWEI(J)=DXX(J)/FSUM
!         WRITE(*,'(I2,A3,F6.3,A3,F5.3,A1)')
!    +     J,' : ',WEI(J),' (',DWEI(J),')'
        END DO
C
100     FORMAT(A,$)
101     FORMAT(A)
        END
C
C******************************************************************************
C
        REAL FUNCTION CHISQT(XX)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INTEGER NMAX_TEMPLATES
        PARAMETER (NMAX_TEMPLATES=20)
        REAL XX(NMAX_TEMPLATES)
C
        INTEGER I,J,NCHF,NTEM
        REAL GAL(NMAXFFT),TEMPL(NMAXFFT,NMAX_TEMPLATES)
        DOUBLE PRECISION DCHISQT,SED,TDIF

        COMMON/BLKTEMPFIT/NCHF,NTEM,GAL,TEMPL
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        DO J=1,NTEM
          IF(XX(J).LT.0.)THEN                     !no admitimos pesos negativos
            CHISQT=1.E25
            RETURN
          END IF
        END DO
C
        DCHISQT=0.D0
        DO I=1,NCHF
          SED=0.D0
          DO J=1,NTEM
            SED=SED+DBLE(XX(J)*TEMPL(I,J))
          END DO
          TDIF=DBLE(GAL(I))-SED
          DCHISQT=DCHISQT+TDIF*TDIF
        END DO
        CHISQT=REAL(DCHISQT)
C
        END
