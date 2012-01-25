C------------------------------------------------------------------------------
C Version 07-September-2007                                   file: midelines.f
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
C Program: midelines
C Classification: measurement
C Description: Measures lines (EW) in a spectrum selecting regions manually.
C
Comment
C
        PROGRAM MIDELINES
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
        INTEGER I,J,L,L1,L2
        INTEGER NS0
        INTEGER NSCAN2,NCHAN2
        INTEGER NB,NBLOCAL
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        INTEGER NTERMS_CONT
        INTEGER NFIT
        INTEGER NITER
        REAL A(NCMAX,NSMAX),ERRA(NCMAX,NSMAX)
        REAL FLUX(NCMAX),FLUXNOR
        REAL X(NCMAX),SS(NCMAX),ESS(NCMAX)
        REAL SKY(NCMAX),SKYEFF(NCMAX),FACTORSKY
        REAL RVGAL
        REAL STWV2,DISP2
        REAL XC,YC
        REAL XMIN,XMAX,YMIN,YMAX,DX,DY
        REAL XMINL,XMAXL
        REAL WLMIN,WLMAX,RCVEL1
        REAL YRMSTOL
        CHARACTER*1 CERR,CFLUX,CRCONT,CH,CSURE,COUT,CCONT,CSKY
        CHARACTER*1 CTYPEFIT
        CHARACTER*50 CDUMMY
        CHARACTER*75 INFILE,ERRINFILE,FLUXFILE,OUTFILE
        CHARACTER*75 SKYFILE
        LOGICAL LBEXIST,LBUFFER
        LOGICAL LLOG,LFIRST
        LOGICAL LPBARS
        LOGICAL LCOLOR(MAX_ID_RED)
        LOGICAL LPLINES
        LOGICAL LYLIMAUTO
        LOGICAL IFCHAN_CONT(NCMAX)
        LOGICAL IFCHAN_LINE(NCMAX)
C
        COMMON/BLKNFIT/NFIT
        COMMON/BLKGEN1/STWV,DISP
        COMMON/BLKGEN2/RCVEL1,WLMIN,WLMAX
        COMMON/BLKGEN3/YMIN,YMAX
        COMMON/BLKGEN4/COUT
        COMMON/BLKMWIDTH1/CRCONT
        COMMON/BLKMWIDTH2/IFCHAN_CONT,IFCHAN_LINE
        COMMON/BLKMWIDTH3/NTERMS_CONT
        COMMON/BLKMWIDTH4/NS0
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
        LLOG=.FALSE.                            !log file does not exist (yet!)
        LBUFFER=.FALSE.
        LFIRST=.TRUE.
        LPBARS=.FALSE.
        LPLINES=.TRUE.
        LYLIMAUTO=.TRUE.
        CRCONT='n'
        YRMSTOL=1.E-4
        THISPROGRAM='midelines'
        CALL WELCOME('07-September-2001')
C------------------------------------------------------------------------------
C abrimos salida grafica
        CALL RPGBEGIN(NTERM,IDN,LCOLOR)
        CALL BUTTSXB(6)
        CALL BUTTSYB(3)
C------------------------------------------------------------------------------
        WRITE(*,100)'Work with error images (y/n) '
        CERR(1:1)=READC('y','yn')
        IF(CERR.EQ.'y')THEN
          NITER=100
        ELSE
          NITER=0
        END IF
        WRITE(*,100) 'Are you flux calibrating the data (y/n) '
        CFLUX(1:1)=READC('n','yn')
        WRITE(*,*)
        WRITE(*,101) 'WARNING: it is important to normalize the data'
        WRITE(*,101) 'in order to avoid errors in the numerical fitting'
        WRITE(*,101) 'of gaussians to the data. This is specially an'
        WRITE(*,101) 'issue when the input spectra are in absolute'
        WRITE(*,101) 'fluxes (~10**-17 erg s**-1 cm**2 A**-1).'
        WRITE(*,100) 'Factor of normalization to be applied '
        FLUXNOR=READF('1.0')
C------------------------------------------------------------------------------
C Leemos fichero(s) de entrada
5       WRITE(*,100)'Input file name'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        DO I=1,NSCAN
          READ(20) (A(J,I),J=1,NCHAN)
          DO J=1,NCHAN
            A(J,I)=A(J,I)/FLUXNOR
          END DO
        END DO
        CLOSE(20)
        IF(CERR.EQ.'y')THEN
          WRITE(*,100)'Input error file name '
          CALL GUESSEF(INFILE,ERRINFILE)
          ERRINFILE=INFILEX(21,ERRINFILE,NSCAN,NCHAN,STWV,DISP,
     +     21,.TRUE.)                                                !match
          DO I=1,NSCAN
            READ(21) (ERRA(J,I),J=1,NCHAN)
            DO J=1,NCHAN
              ERRA(J,I)=ERRA(J,I)/FLUXNOR
            END DO
          END DO
          CLOSE(21)
        END IF
C
        WRITE(*,100) 'Are you overplotting a sky spectrum (y/n) '
        CSKY(1:1)=READC('y','yn')
        IF(CSKY.EQ.'y')THEN
           WRITE(*,100)'Input sky file name (single spectrum)'
          SKYFILE=INFILEX(22,'@',
     +     NSCAN2,NCHAN2,STWV2,DISP2,1,.FALSE.)
          IF(NSCAN2.GT.1)THEN
            WRITE(*,101)'WARNING: this file contains more than 1 '//
     +       'spectrum.'
            WRITE(*,101)'The first spectrum will be used as '//
     +       'sky spectrum.'
          END IF
          IF(NCHAN2.NE.NCHAN)THEN
            CLOSE(22)
            WRITE(*,100)'FATAL ERROR: NCHAN in last image is different!'
            STOP
          END IF
          IF(STWV2.NE.STWV)THEN
            CLOSE(22)
            WRITE(*,100)'FATAL ERROR: STWV in last image is different!'
            STOP
          END IF
          IF(DISP2.NE.DISP)THEN
            CLOSE(22)
            WRITE(*,100)'FATAL ERROR: DISP in last image is different!'
            STOP
          END IF
          READ(22) (SKY(J),J=1,NCHAN)
          DO J=1,NCHAN
            SKY(J)=SKY(J)/FLUXNOR
          END DO
          CLOSE(22)
          WRITE(*,100) 'Factor to be applied to the sky signal '
          FACTORSKY=READF('100.0')
        END IF
C
        WRITE(*,100)'Radial velocity (km/sec)'
        RVGAL=READF('@')
C
        IF(.NOT.LLOG)THEN
          WRITE(*,100)'Create output file (y/n) '
          COUT(1:1)=READC('n','yn')
          IF(COUT.EQ.'y')THEN
            LLOG=.TRUE.
            WRITE(*,100)'Output file name'
            OUTFILE=OUTFILEX(30,'@',0,0,0.,0.,3,.FALSE.)
          END IF
        END IF
C
        IF(LLOG)THEN
          WRITE(30,150)
          WRITE(30,150)
          WRITE(30,'(A,A)')'This  file: ',OUTFILE(1:TRUELEN(OUTFILE))
          WRITE(30,'(A,A)')'Input file: ',INFILE(1:TRUELEN(INFILE))
          IF(CERR.EQ.'y') WRITE(30,'(A,A)')'Error file: ',
     +     ERRINFILE(1:TRUELEN(ERRINFILE))
          IF(CSKY.EQ.'y') WRITE(30,'(A,A)')'Sky file..: ',
     +     SKYFILE(1:TRUELEN(SKYFILE))
          WRITE(30,100)'STWV..: '
          WRITE(30,*) STWV
          WRITE(30,100)'DISP..: '
          WRITE(30,*) DISP
          WRITE(30,101)'Object: '//OBJECT(1:TRUELEN(OBJECT))
          WRITE(30,100)'Radial velocity (km/sec): '
          WRITE(30,*) RVGAL
        END IF
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
C calibracion en flujo
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
          IF(COUT.EQ.'y')THEN
            WRITE(30,101)'Response curve file name: '//FLUXFILE
          END IF
        ELSE
          DO J=1,NCHAN
            FLUX(J)=1.0
          END DO
        END IF
C------------------------------------------------------------------------------
C calculamos el intervalo disponible en longitud de onda
        IF((STWV.EQ.0.0).AND.(DISP.EQ.0.0))THEN
          WRITE(*,101) 'WARNING: there is no wavelength calibration.'
          WRITE(*,100) 'Do you want to continue anyway (y/n) '
          CCONT(1:1)=READC('y','yn')
          IF(CCONT.EQ.'n')THEN
            CALL PGEND
            STOP
          END IF
          STWV=1.
          DISP=1.
        END IF
        WLMIN=STWV-DISP/2.
        WLMAX=STWV+REAL(NCHAN-1)*DISP+DISP/2.0
        WRITE(*,100)'>>> WLMIN..: '
        WRITE(*,*) WLMIN
        WRITE(*,100)'>>> WLMAX..: '
        WRITE(*,*) WLMAX
        RCVEL1=1.0+RVGAL/C
        RCVEL1=RCVEL1/SQRT(1.0-(RVGAL*RVGAL)/(C*C))     !correccion relativista
C------------------------------------------------------------------------------
        IF(LFIRST)THEN
          LFIRST=.FALSE.
        ELSE
          CALL PGERAS
        END IF
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
        CALL BUTTON(7,'[m]easure',0)
        CALL BUTTON(8,'[r]epeat',0)
        IF(.NOT.LBUFFER) CALL BUTTON(8,'[r]epeat',3)
        CALL BUTTON(9,'repeat[?]',0)
        IF(.NOT.LBUFFER) CALL BUTTON(9,'repeat[?]',3)
        CALL BUTTON(10,'pol.[d]eg.',0)
        IF(.NOT.LBUFFER) CALL BUTTON(10,'pol.[d]eg.',3)
        CALL BUTTON(11,'new [f]ile',0)
        CALL BUTTON(12,'E[x]it',0)
        CALL BUTTON(13,'err[.]bars',0)
        IF(CERR.EQ.'n') CALL BUTTON(13,'err[.]bars',3)
        CALL BUTTON(14,'plot lines',0)
        IF(LPLINES) CALL BUTTON(14,'plot lines',1)
        IF(LYLIMAUTO)THEN
          CALL BUTTON(18,'Ylim auto',0)
          CALL BUTTON(18,'Ylim auto',1)
        ELSE
          CALL BUTTON(18,'Ylim fixed',0)
          CALL BUTTON(18,'Ylim fixed',1)
        END IF
C------------------------------------------------------------------------------
        XMIN=1.
        XMAX=REAL(NCHAN)
        DX=XMAX-XMIN
        XMIN=XMIN-DX/20.
        XMAX=XMAX+DX/20.
C..............................................................................
10      DO J=1,NCHAN
          SS(J)=A(J,NS0)/FLUX(J)
        END DO
        IF(CERR.EQ.'y')THEN
          DO J=1,NCHAN
            ESS(J)=ERRA(J,NS0)/FLUX(J)
          END DO
        ELSE
          DO J=1,NCHAN
            ESS(J)=0.0
          END DO
        END IF
20      IF(LYLIMAUTO)THEN
          YMIN=1.E20
          YMAX=-1.E20
          IF((CERR.EQ.'y').AND.LPBARS)THEN
            DO J=1,NCHAN
              IF((X(J).GE.XMIN).AND.(X(J).LE.XMAX))THEN
                IF(SS(J)-ESS(J).LT.YMIN) YMIN=SS(J)-ESS(J)
                IF(SS(J)+ESS(J).GT.YMAX) YMAX=SS(J)+ESS(J)
              END IF
            END DO
          ELSE
            DO J=1,NCHAN
              IF((X(J).GE.XMIN).AND.(X(J).LE.XMAX))THEN
                IF(SS(J).LT.YMIN) YMIN=SS(J)
                IF(SS(J).GT.YMAX) YMAX=SS(J)
              END IF
            END DO
          END IF
          DY=YMAX-YMIN
          YMIN=YMIN-DY/20.
          YMAX=YMAX+DY/20.
        END IF
ccc30      XMINL=(XMIN-1.)*DISP+STWV
        XMINL=(XMIN-1.)*DISP+STWV
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
          CALL PGBIN(NCHAN,X,SS,.TRUE.)
          IF(LPLINES) CALL SUBPLINES(.TRUE.)
          IF(CSKY.EQ.'y')THEN
            DY=YMAX-YMIN
            DO J=1,NCHAN
              SKYEFF(J)=SKY(J)/FACTORSKY-DY/2.5
            END DO
            CALL PGSCI(2)
            CALL PGBIN(NCHAN,X,SKYEFF,.TRUE.)
            CALL PGSCI(1)
          END IF
          IF((CERR.EQ.'y').AND.LPBARS)THEN
            CALL PGSCI(5)
            DO J=1,NCHAN
              CALL PGERRY(1,X(J),SS(J)-ESS(J),SS(J)+ESS(J),1.0)
            END DO
            CALL PGSCI(1)
          END IF
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
          CALL PGLABEL('Channel','Flux units',CHAR(32))
          WRITE(CDUMMY,*) NS0
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          CALL PGMTEXT('T',2.5,1.0,1.0,'Scan #'//CDUMMY(1:L))
          L1=TRUELEN(INFILE)
          L2=TRUELEN(OBJECT)
          CALL PGMTEXT('T',2.5,0.,0.,'file: '//INFILE(1:L1)//' '//
     +     OBJECT(1:L2))
        END DO
C..............................................................................
50      CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
        CALL CHLOWER(CH)
        CALL IFBUTTON(XC,YC,NB)
C
        NBLOCAL=INDEX('zkwjpnmr?dfx.',CH)
        IF((NBLOCAL.NE.0).AND.(CH.NE.' '))THEN
          CALL BUTTQEX(NBLOCAL,LBEXIST)
          IF(LBEXIST) NB=NBLOCAL
        END IF
C..............................................................................
        IF(NB.EQ.0)THEN
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
C..............................................................................
        ELSEIF(NB.EQ.1)THEN
          CALL BUTTON(1,'[z]oom (m)',5)
          WRITE(*,100)'Press mouse (limit #1)...'
          IF(LCOLOR(1)) CALL PGSCI(2)
          CALL RPGBAND(6,0,0.,0.,XC,YC,CH)
          XMIN=XC
          IF(XMIN.LT.1.) XMIN=1.
          IF(XMIN.GT.REAL(NCHAN)) XMIN=REAL(NCHAN)
          WRITE(*,100)'   cursor at '
          WRITE(*,*)XMIN
          WRITE(*,100)'Press mouse (limit #2)...'
          CALL RPGBAND(4,0,XMIN,0.,XC,YC,CH)
          IF(LCOLOR(1)) CALL PGSCI(1)
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
          CALL BUTTON(4,'[j]ump',0)
          CALL RPGERASW(0.,1.,0.,0.80)
          GOTO 10
C..............................................................................
        ELSEIF(NB.EQ.5)THEN
          CALL BUTTON(5,'[p]revious',5)
          IF(NS0.GT.1) NS0=NS0-1
          CALL BUTTON(5,'[p]revious',0)
          CALL RPGERASW(0.,1.,0.,0.80)
          GOTO 10
C..............................................................................
        ELSEIF(NB.EQ.6)THEN
          CALL BUTTON(6,'[n]ext',5)
          IF(NS0.LT.NSCAN) NS0=NS0+1
          CALL BUTTON(6,'[n]ext',0)
          CALL RPGERASW(0.,1.,0.,0.80)
          GOTO 10
C..............................................................................
        ELSEIF(NB.EQ.7)THEN
          CALL BUTTON(7,'[m]easure',5)
          WRITE(*,101) '(1) 1 gaussian'
          WRITE(*,100) '(2) 2 gaussian (same width and same area)'
          WRITE(*,101) ' with fixed separation'
          WRITE(*,100) '(3) 2 gaussian (same width and different area)'
          WRITE(*,101) ' with fixed separation'
          WRITE(*,101) '(4) 1 gaussian with fixed x0 and sigma'
          WRITE(*,100) 'Option '
          CTYPEFIT(1:1)=READC('1','1234')
          CALL MIDEWIDTH(CTYPEFIT,NCHAN,SS,ESS,CERR,0,NITER,YRMSTOL)
          LBUFFER=.TRUE.
          CALL BUTTON(8,'[r]epeat',0)
          CALL BUTTON(9,'repeat[?]',0)
          CALL BUTTON(7,'[m]easure',0)
          CALL BUTTON(10,'pol.[d]eg.',0)
C..............................................................................
        ELSEIF(NB.EQ.8)THEN
          CALL BUTTON(8,'[r]epeat',5)
          CALL MIDEWIDTH(CTYPEFIT,NCHAN,SS,ESS,CERR,1,NITER,YRMSTOL)
          CALL BUTTON(8,'[r]epeat',0)
C..............................................................................
        ELSEIF(NB.EQ.9)THEN
          CALL BUTTON(9,'repeat[?]',5)
          CALL MIDEWIDTH(CTYPEFIT,NCHAN,SS,ESS,CERR,2,NITER,YRMSTOL)
          CALL BUTTON(9,'repeat[?]',0)
C..............................................................................
        ELSEIF(NB.EQ.10)THEN
          CALL BUTTON(10,'pol.[d]eg.',5)
          WRITE(*,100)'Pol. degree '
          NTERMS_CONT=READILIM('0',0,MIN(19,NFIT-1))
          IF(COUT.EQ.'y')THEN
            WRITE(30,100)'Polynomial degree: '
            WRITE(30,*) NTERMS_CONT
          END IF
          NTERMS_CONT=NTERMS_CONT+1
          CALL BUTTON(10,'pol.[d]eg.',0)
C..............................................................................
        ELSEIF(NB.EQ.11)THEN
          CALL BUTTON(11,'new [f]ile',5)
          CALL BUTTON(11,'new [f]ile',0)
          GOTO 5
C..............................................................................
        ELSEIF(NB.EQ.12)THEN
          CALL BUTTON(12,'E[x]it',5)
          WRITE(*,100)'Do you really want to exit (y/n) '
          CSURE(1:1)=READC('n','yn')
          IF(CSURE.EQ.'n')THEN
            CALL BUTTON(12,'E[x]it',0)
          ELSE
            CALL PGEND
            STOP
          END IF
C..............................................................................
        ELSEIF(NB.EQ.13)THEN
          IF(LPBARS)THEN
            CALL BUTTON(13,'err[.]bars',0)
            LPBARS=.FALSE.
          ELSE
            CALL BUTTON(13,'err[.]bars',1)
            LPBARS=.TRUE.
          END IF
          CALL RPGERASW(0.,1.,0.,0.80)
          GOTO 10
C..............................................................................
        ELSEIF(NB.EQ.14)THEN
          IF(LPLINES)THEN
            LPLINES=.FALSE.
            CALL BUTTON(14,'plot lines',0)
          ELSE
            CALL BUTTON(14,'plot lines',1)
            LPLINES=.TRUE.
          END IF
          CALL RPGERASW(0.,1.,0.,0.80)
          GOTO 20
C..............................................................................
        ELSEIF(NB.EQ.15)THEN
          CALL BUTTON(15,'Sky Factor',5)
          WRITE(CDUMMY,*) FACTORSKY
          WRITE(*,100) 'Factor to be applied to the sky signal'
          FACTORSKY=READF(CDUMMY)
          CALL BUTTON(15,'Sky Factor',0)
C..............................................................................
        ELSEIF(NB.EQ.18)THEN
          IF(LYLIMAUTO)THEN
            LYLIMAUTO=.FALSE.
            CALL BUTTON(18,'Ylim fixed',0)
            CALL BUTTON(18,'Ylim fixed',1)
            WRITE(CDUMMY,*) YMIN
            WRITE(*,100) 'Ymin'
            YMIN=READF(CDUMMY)
            WRITE(CDUMMY,*) YMAX
            WRITE(*,100) 'Ymax'
            YMAX=READF(CDUMMY)
          ELSE
            LYLIMAUTO=.TRUE.
            CALL BUTTON(18,'Ylim auto',0)
            CALL BUTTON(18,'Ylim auto',1)
          END IF
          CALL RPGERASW(0.,1.,0.,0.80)
          GOTO 20
C..............................................................................
        END IF
        GOTO 50
C------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
150     FORMAT(79('-'))
        END
C
C******************************************************************************
C Define que canales van a ser utilizados para normalizar los espectros.
C Si LADD=.TRUE. la subrutina no inicializa el array IFCHAN, pudiendose asi
C mantener regiones previamente seleccionadas
C
        SUBROUTINE SUBREGIONS(IFCHAN,NCHAN,LADD,NCOLOR)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'futils.inc'
C
        LOGICAL IFCHAN(NCHAN),LADD
C
        INTEGER J,L
        INTEGER NC1,NC2,NC0
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        INTEGER NCOLOR
        REAL XC,YC
        REAL YMIN,YMAX,DY
        CHARACTER*1 CMODE,CH,COUT
        CHARACTER*50 CDUMMY
        LOGICAL LCOLOR(MAX_ID_RED)
C
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
        COMMON/BLKGEN3/YMIN,YMAX
        COMMON/BLKGEN4/COUT
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
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
10      IF(CMODE.EQ.'m')THEN
          WRITE(*,101)'Press <q>/<x>/<mouse right button> to EXIT'
          WRITE(*,100)'Press mouse (limit #1)...'
          IF(LCOLOR(1)) CALL PGSCI(2)
          CALL RPGBAND(6,0,0.,0.,XC,YC,CH)
          CALL CHLOWER(CH)
          IF((CH.EQ.'x').OR.(CH.EQ.'q'))THEN
            WRITE(*,*)
            GOTO 20
          END IF
          NC1=NINT(XC)
          IF(NC1.LT.1) NC1=1
          IF(NC1.GT.NCHAN) NC1=NCHAN
          WRITE(*,100)'   cursor at '
          WRITE(*,*)NC1
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(NCOLOR)
            CALL PGMOVE(REAL(NC1),YMIN)
            CALL PGDRAW(REAL(NC1),YMAX)
          END DO
          WRITE(*,100)'Press mouse (limit #2)...'
          IF(LCOLOR(1)) CALL PGSCI(2)
          CALL RPGBAND(6,0,0.,0.,XC,YC,CH)
          CALL CHLOWER(CH)
          IF((CH.EQ.'x').OR.(CH.EQ.'q'))THEN
            WRITE(*,*)
            GOTO 20
          END IF
          NC2=NINT(XC)
          IF(NC2.LT.1) NC2=1
          IF(NC2.GT.NCHAN) NC2=NCHAN
          WRITE(*,100)'   cursor at '
          WRITE(*,*)NC2
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(NCOLOR)
            CALL PGMOVE(REAL(NC2),YMIN)
            CALL PGDRAW(REAL(NC2),YMAX)
          END DO
          IF(NC1.GT.NC2)THEN
            NC0=NC1
            NC1=NC2
            NC2=NC0
          END IF
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGRECT(REAL(NC1),REAL(NC2),YMIN,YMIN+DY/100.)
            CALL PGRECT(REAL(NC1),REAL(NC2),YMAX,YMAX-DY/100.)
            IF(LCOLOR(1)) CALL PGSCI(ITERM)
          END DO
        ELSE
12        WRITE(*,100)'Channel region (0,0=EXIT) '
          CALL READ2I('0,0',NC1,NC2)
          IF((NC1.EQ.0).AND.(NC2.EQ.0)) GOTO 20
          IF((NC1.LT.1).OR.(NC2.GT.NCHAN).OR.(NC1.GT.NC2))THEN
            WRITE(*,101)'ERROR: invalid numbers. Try again.'
            GOTO 12
          END IF
        END IF
        DO J=NC1,NC2
          IFCHAN(J)=.TRUE.
        END DO
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM)) CALL PGSCI(NCOLOR)
          CALL PGMOVE(REAL(NC1),YMIN)
          CALL PGDRAW(REAL(NC1),YMAX)
          CALL PGMOVE(REAL(NC2),YMIN)
          CALL PGDRAW(REAL(NC2),YMAX)
          CALL PGRECT(REAL(NC1),REAL(NC2),YMIN,YMIN+DY/100.)
          CALL PGRECT(REAL(NC1),REAL(NC2),YMAX,YMAX-DY/100.)
          IF(LCOLOR(1)) CALL PGSCI(ITERM)
        END DO
        IF(COUT.EQ.'y')THEN
          WRITE(CDUMMY,'(A1,I6,A1,I6,A1)')'(',NC1,',',NC2,')'
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          WRITE(30,100)CDUMMY(1:L)
        END IF
        GOTO 10
C
20      NC0=0
        DO J=1,NCHAN
          IF(IFCHAN(J)) NC0=NC0+1
        END DO
        IF(NC0.EQ.0)THEN
          WRITE(*,101) 'ERROR: number of channels to be used = 0! '//
     +     'Try again.'
          GOTO 10
        END IF
C
100     FORMAT(A,$)
101     FORMAT(A)
        END
C
C******************************************************************************
C Mide las lineas
C IREPEAT=0: mide desde el principio (no es repeticion)
C IREPEAT=1: es repeticion (no ha de preguntar nada)
C IREPEAT=2: es repeticion (pero puede preguntar algo)
        SUBROUTINE MIDEWIDTH(CTYPEFIT,NCHAN,Y,YE,CERR,IREPEAT,
     +   NITER,YRMSTOL)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'futils.inc'
        INTEGER READILIM
        REAL READF
C
        CHARACTER*1 CTYPEFIT
        REAL Y(NCHAN),YE(NCHAN)
        CHARACTER*1 CERR
        INTEGER IREPEAT,NITER
        REAL YRMSTOL
C
        INTEGER NMAXSIMUL
        PARAMETER (NMAXSIMUL=1000)
        INTEGER NDEGMAX
        PARAMETER (NDEGMAX=16)
        REAL PI
        PARAMETER(PI= 3.141592654)
        REAL PI2
        PARAMETER(PI2=6.283185307)
C
        INTEGER J,JMIN,JMAX,K,L,NSIM
        INTEGER NFIT,NPLOT,NTERMS_CONT
        INTEGER NS0
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        INTEGER NSEED
        INTEGER NCHAN_CONTINUUM
        REAL X(NCMAX)
        REAL XFIT(NCMAX),YFIT(NCMAX),EYFIT(NCMAX)
        REAL XPLOT(NCMAX),YPLOT(NCMAX)
        REAL EYPLOT1(NCMAX),EYPLOT2(NCMAX)
        REAL XX,YY,EYY,YY1,EYY1,YY2,EYY2,POL(NCMAX),EPOL(NCMAX)
        REAL ERRYY
        REAL A(NDEGMAX+1),CHISQR,A_(NDEGMAX+1)
        REAL POL_SIMUL(NCMAX,NMAXSIMUL),POL_SIMUL_YY(NMAXSIMUL)
        REAL POL_SIMUL_YY1(NMAXSIMUL),POL_SIMUL_YY2(NMAXSIMUL)
        REAL WIDTH,ERRWIDTH
        REAL FLUX1,FLUX2,ERRFLUX1,ERRFLUX2
        REAL FLUXC,ERRFLUXC !flux in the continuum and error
        REAL FLUXG,EFLUXG
        REAL EW,ERREW,EW1,ERREW1,EW2,ERREW2
        REAL FWHM
        REAL AMP,AMP1,AMP2,X0,SIGMA
        REAL EAMP,EAMP1,EAMP2,EX0,ESIGMA
        REAL EEAMP,EEAMP1,EEAMP2,EEX0,EESIGMA
        REAL MINSIGMA
        REAL RCVEL1,WLMIN,WLMAX
        REAL DELTAX
        REAL R1,R2
        REAL RANRED
        CHARACTER*1 CRCONT,CBEFORE,COUT
        CHARACTER*50 CDUMMY
        LOGICAL IFCHAN_CONT(NCMAX)
        LOGICAL IFCHAN_LINE(NCMAX)
        LOGICAL LCOLOR(MAX_ID_RED)
C
        COMMON/BLKNFIT/NFIT
        COMMON/BLKGEN1/STWV,DISP
        COMMON/BLKGEN2/RCVEL1,WLMIN,WLMAX
        COMMON/BLKGEN4/COUT
        COMMON/BLKMWIDTH1/CRCONT
        COMMON/BLKMWIDTH2/IFCHAN_CONT,IFCHAN_LINE
        COMMON/BLKMWIDTH3/NTERMS_CONT
        COMMON/BLKMWIDTH4/NS0
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
        JMIN=0 !avoid compilation warning
        JMAX=0 !avoid compilation warning
C
        IF(COUT.EQ.'y')THEN
          WRITE(30,100)'>>> Current SCAN is #'
          WRITE(CDUMMY,*)NS0
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          WRITE(30,101)CDUMMY(1:L)
        END IF
C
        IF(IREPEAT.NE.0)THEN
          CBEFORE='y'
        ELSE
          CBEFORE='n'
        END IF
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
        IF(CBEFORE.EQ.'n')THEN
          WRITE(*,100) 'Remove continuum (y/n) '
          CRCONT(1:1)=READC(CRCONT,'yn')
          IF(CRCONT.EQ.'y')THEN
            WRITE(*,*)
            WRITE(*,101)'* Select channels to be employed to fit '//
     +       'continuum:'
            IF(COUT.EQ.'y')THEN
              WRITE(30,100)'Continuum regions: '
            END IF
            CALL SUBREGIONS(IFCHAN_CONT,NCHAN,.FALSE.,4)
            IF(COUT.EQ.'y') WRITE(30,*)
          END IF
        ELSE
          IF(CRCONT.EQ.'y')THEN
            IF(COUT.EQ.'y')THEN
              WRITE(30,101)'Continuum regions: (previous)'
            END IF
          END IF
        END IF
C
        IF(CRCONT.EQ.'y')THEN
          NFIT=0
          DO J=1,NCHAN
            IF(IFCHAN_CONT(J))THEN
              NFIT=NFIT+1
              XFIT(NFIT)=X(J)
              YFIT(NFIT)=Y(J)
            END IF
          END DO
          NCHAN_CONTINUUM=NFIT
          IF(CBEFORE.EQ.'n')THEN
            WRITE(*,100) 'Pol. degree '
            NTERMS_CONT=READILIM('0',0,MIN(16,NFIT-1))
            IF(COUT.EQ.'y')THEN
              WRITE(30,100) 'Polynomial degree: '
              WRITE(30,*) NTERMS_CONT
            END IF
            NTERMS_CONT=NTERMS_CONT+1
          ELSE
            IF(COUT.EQ.'y')THEN
              WRITE(30,101) 'Polynomial degree: (previous)'
            END IF
          END IF
          CALL POLFIT(XFIT,YFIT,YFIT,NFIT,NTERMS_CONT,0,A,CHISQR)
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
        ELSE
          JMIN=1
          JMAX=NCHAN
          NTERMS_CONT=1
          A(1)=0.
          NCHAN_CONTINUUM=0
        END IF
C calculamos continuo y errores en el continuo
        DO J=1,NCHAN
          XX=REAL(J)
          POL(J)=A(NTERMS_CONT)
          IF(NTERMS_CONT.GT.1)THEN
            DO K=NTERMS_CONT-1,1,-1
              POL(J)=POL(J)*XX+A(K)
            END DO
          END IF
        END DO
        IF(CERR.EQ.'y')THEN
          IF(IREPEAT.NE.1)THEN
            WRITE(CDUMMY,*) NITER
            WRITE(*,100) 'No of iterations to evaluate errors '//
     +       'in continuum'
            NITER=READILIM(CDUMMY,2,NMAXSIMUL)
          END IF
          NSEED=-1
          DO NSIM=1,NITER
            NFIT=0
            DO J=1,NCHAN
              IF(IFCHAN_CONT(J))THEN
                NFIT=NFIT+1
                XFIT(NFIT)=X(J)
                R1=RANRED(NSEED)
                R2=RANRED(NSEED)
                ERRYY=1.414213562*YE(J)*SQRT(-1.*LOG(1.-R1))*
     +           COS(PI2*R2)
                YFIT(NFIT)=Y(J)+ERRYY
              END IF
            END DO
            CALL POLFIT(XFIT,YFIT,YFIT,NFIT,NTERMS_CONT,0,A_,CHISQR)
            DO J=1,NCHAN
              XX=REAL(J)
              POL_SIMUL(J,NSIM)=A_(NTERMS_CONT)
              IF(NTERMS_CONT.GT.1)THEN
                DO K=NTERMS_CONT-1,1,-1
                  POL_SIMUL(J,NSIM)=POL_SIMUL(J,NSIM)*XX+A_(K)
                END DO
              END IF
            END DO
          END DO
          DO J=1,NCHAN
            EPOL(J)=0.0
            DO NSIM=1,NITER
              EPOL(J)=EPOL(J)+
     +         (POL(J)-POL_SIMUL(J,NSIM))*(POL(J)-POL_SIMUL(J,NSIM))
            END DO
            EPOL(J)=SQRT(EPOL(J)/REAL(NSIM-1))
          END DO
        ELSE
          DO J=1,NCHAN
            EPOL(J)=0.0
          END DO
        END IF
C dibujamos el continuo y su error
        NPLOT=0
        DO J=JMIN,JMAX
          NPLOT=NPLOT+1
          XPLOT(NPLOT)=REAL(J)
          YPLOT(NPLOT)=POL(J)
          IF(CERR.EQ.'y')THEN
            EYPLOT1(NPLOT)=POL(J)-EPOL(J)
            EYPLOT2(NPLOT)=POL(J)+EPOL(J)
          END IF
        END DO
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM)) CALL PGSCI(5)
          CALL PGBIN(NPLOT,XPLOT,YPLOT,.TRUE.)
          IF(CERR.EQ.'y')THEN
            IF(LCOLOR(ITERM)) CALL PGSCI(4)
            CALL PGBIN(NPLOT,XPLOT,EYPLOT1,.TRUE.)
            CALL PGBIN(NPLOT,XPLOT,EYPLOT2,.TRUE.)
          END IF
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
        END DO
C
        IF(CBEFORE.EQ.'n')THEN
          WRITE(*,*)
          WRITE(*,101) '* Select channels to be employed to measure '//
     +     'the line:'
          IF(COUT.EQ.'y')THEN
            WRITE(30,100) 'Line regions.....: '
          END IF
          CALL SUBREGIONS(IFCHAN_LINE,NCHAN,.FALSE.,3)
          IF(COUT.EQ.'y') WRITE(30,*)
        ELSE
          IF(COUT.EQ.'y')THEN
            WRITE(30,101) 'Line regions.....: (previous)'
          END IF
        END IF
        NFIT=0
        DO J=1,NCHAN
          IF(IFCHAN_LINE(J))THEN
            NFIT=NFIT+1
            XFIT(NFIT)=X(J)
            YFIT(NFIT)=Y(J)
            EYFIT(NFIT)=YE(J)
          END IF
        END DO
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM)) CALL PGSCI(7)
          CALL PGBIN(NFIT,XFIT,YFIT,.TRUE.)
        END DO
        FLUX1=0.
        FLUXC=0.
        FLUX2=0.
        WIDTH=0.
        IF(CERR.EQ.'y')THEN
          ERRFLUX1=0.
          ERRFLUXC=0.
          ERRFLUX2=0.
          ERRWIDTH=0.
        END IF
        NFIT=0
        DO J=1,NCHAN
          IF(IFCHAN_LINE(J))THEN
            XX=REAL(J)
            POL(J)=A(NTERMS_CONT)
            DO K=NTERMS_CONT-1,1,-1
              POL(J)=POL(J)*XX+A(K)
            END DO
            NFIT=NFIT+1
            YFIT(NFIT)=Y(J)-POL(J)
            FLUX1=FLUX1+(POL(J)-Y(J))*DISP
            FLUXC=FLUXC+POL(J)*DISP
            FLUX2=FLUX2+Y(J)*DISP
            IF(CRCONT.EQ.'y') WIDTH=WIDTH+DISP*(POL(J)-Y(J))/POL(J)
            IF(CERR.EQ.'y')THEN
              ERRFLUX1=ERRFLUX1+(YE(J)*YE(J)+EPOL(J)*EPOL(J))
              ERRFLUXC=ERRFLUXC+EPOL(J)*EPOL(J)
              ERRFLUX2=ERRFLUX2+YE(J)*YE(J)
              IF(CRCONT.EQ.'y')THEN
                ERRWIDTH=ERRWIDTH+(YE(J)*YE(J)+
     +           EPOL(J)*EPOL(J))/(POL(J)*POL(J))+
     +           ((POL(J)-Y(J))/(POL(J)*POL(J)))*
     +           ((POL(J)-Y(J))/(POL(J)*POL(J)))*EPOL(J)*EPOL(J)
              END IF
            END IF
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              CALL PGMOVE(XX,POL(J))
              CALL PGDRAW(XX,Y(J))
            END DO
          END IF
        END DO
        IF(CERR.EQ.'y')THEN
          ERRFLUX1=DISP*SQRT(ERRFLUX1)
          ERRFLUXC=DISP*SQRT(ERRFLUXC)
          ERRFLUX2=DISP*SQRT(ERRFLUX2)
          IF(CRCONT.EQ.'y') ERRWIDTH=DISP*SQRT(ERRWIDTH)
        ELSE
          ERRFLUX1=0.
          ERRFLUXC=0.
          ERRFLUX2=0.
          ERRWIDTH=0.
        END IF
        IF(IREPEAT.NE.1)THEN
          IF((CTYPEFIT.EQ.'1').OR.(CTYPEFIT.EQ.'2').OR.
     +     (CTYPEFIT.EQ.'3'))THEN
            WRITE(CDUMMY,*)YRMSTOL
            WRITE(*,100)'YRMSTOL for DOWNHILL '
            YRMSTOL=READF(CDUMMY)
          ELSEIF(CTYPEFIT.EQ.'4')THEN
            WRITE(*,100) 'X0 (pixels)'
            X0=READF('@')
            WRITE(*,100) 'Error in X0 '
            EX0=READF('0.0')
            WRITE(*,100) 'SIGMA (pixels)'
            SIGMA=READF('@')
            WRITE(*,100) 'Error in SIGMA '
            ESIGMA=READF('0.0')
            WRITE(*,100) 'Minimum SIGMA (spectral resolution)'
            MINSIGMA=READF('@')
          END IF
        END IF
        IF((CTYPEFIT.EQ.'2').OR.(CTYPEFIT.EQ.'3'))THEN
          WRITE(*,100) 'Separation between gaussians (pixels)'
          DELTAX=READF('@')
        ELSE
          DELTAX=0.0
        END IF
        IF(CERR.EQ.'y')THEN
          IF(IREPEAT.NE.1)THEN
            WRITE(CDUMMY,*) NITER
            WRITE(*,100)'No of iterations to evaluate errors '
            NITER=READILIM(CDUMMY,2,NMAXSIMUL)
          END IF
          IF(CTYPEFIT.EQ.'1')THEN
            CALL GAUSSFIT(NFIT,XFIT,YFIT,EYFIT,X0,SIGMA,AMP,
     +       EX0,ESIGMA,EAMP,EEX0,EESIGMA,EEAMP,YRMSTOL,NITER)
          ELSEIF(CTYPEFIT.EQ.'2')THEN
            CALL GAUSS2AFIT(NFIT,XFIT,YFIT,EYFIT,DELTAX,X0,SIGMA,AMP,
     +       EX0,ESIGMA,EAMP,EEX0,EESIGMA,EEAMP,YRMSTOL,NITER)
            AMP1=AMP
            AMP2=AMP
            EAMP1=EAMP
            EAMP2=EAMP
            EEAMP1=EEAMP
            EEAMP2=EEAMP
          ELSEIF(CTYPEFIT.EQ.'3')THEN
            CALL GAUSS2BFIT(NFIT,XFIT,YFIT,EYFIT,DELTAX,
     +       X0,SIGMA,AMP1,AMP2,
     +       EX0,ESIGMA,EAMP1,EAMP2,EEX0,EESIGMA,EEAMP1,EEAMP2,
     +       YRMSTOL,NITER)
          ELSEIF(CTYPEFIT.EQ.'4')THEN
            CALL GAUSSFITAMP(NFIT,XFIT,YFIT,EYFIT,X0,SIGMA,
     +       EX0,ESIGMA,MINSIGMA,AMP,EAMP,EEAMP,NITER)
            EEX0=0.0
            EESIGMA=0.0
          END IF
        ELSE
          IF(CTYPEFIT.EQ.'1')THEN
            CALL GAUSSFIT(NFIT,XFIT,YFIT,EYFIT,X0,SIGMA,AMP,
     +       EX0,ESIGMA,EAMP,EEX0,EESIGMA,EEAMP,YRMSTOL,0)
          ELSEIF(CTYPEFIT.EQ.'2')THEN
            CALL GAUSS2AFIT(NFIT,XFIT,YFIT,EYFIT,DELTAX,X0,SIGMA,AMP,
     +       EX0,ESIGMA,EAMP,EEX0,EESIGMA,EEAMP,YRMSTOL,0)
          ELSEIF(CTYPEFIT.EQ.'3')THEN
            CALL GAUSS2BFIT(NFIT,XFIT,YFIT,EYFIT,DELTAX,
     +       X0,SIGMA,AMP1,AMP2,
     +       EX0,ESIGMA,EAMP1,EAMP2,EEX0,EESIGMA,EEAMP1,EEAMP2,
     +       YRMSTOL,0)
          ELSEIF(CTYPEFIT.EQ.'4')THEN
            DO J=1,NFIT
              EYFIT(J)=1.0 !ponemos todos los errores a uno para ajuste pesado
            END DO
            CALL GAUSSFITAMP(NFIT,XFIT,YFIT,EYFIT,X0,SIGMA,
     +       EX0,ESIGMA,MINSIGMA,AMP,EAMP,EEAMP,0)
            EEX0=0.0
            EESIGMA=0.0
          END IF
        END IF
C
        IF((CTYPEFIT.EQ.'1').OR.(CTYPEFIT.EQ.'4'))THEN
          XX=X0
          YY=A(NTERMS_CONT)    !valor del continuo en el centro de la gaussiana
          DO K=NTERMS_CONT-1,1,-1
            YY=YY*XX+A(K)
          END DO
        ELSE
          XX=X0
          YY1=A(NTERMS_CONT) !valor del continuo en el centro de la gaussiana#1
          DO K=NTERMS_CONT-1,1,-1
            YY1=YY1*XX+A(K)
          END DO
          XX=X0+DELTAX
          YY2=A(NTERMS_CONT) !valor del continuo en el centro de la gaussiana#2
          DO K=NTERMS_CONT-1,1,-1
            YY2=YY2*XX+A(K)
          END DO
        END IF
        IF(CERR.EQ.'y')THEN
          DO NSIM=1,NITER
            NFIT=0
            DO J=1,NCHAN
              IF(IFCHAN_CONT(J))THEN
                NFIT=NFIT+1
                XFIT(NFIT)=X(J)
                R1=RANRED(NSEED)
                R2=RANRED(NSEED)
                ERRYY=1.414213562*YE(J)*SQRT(-1.*LOG(1.-R1))*
     +           COS(PI2*R2)
                YFIT(NFIT)=Y(J)+ERRYY
              END IF
            END DO
            CALL POLFIT(XFIT,YFIT,YFIT,NFIT,NTERMS_CONT,0,A_,CHISQR)
            IF((CTYPEFIT.EQ.'1').OR.(CTYPEFIT.EQ.'4'))THEN
              XX=X0
              POL_SIMUL_YY(NSIM)=A_(NTERMS_CONT)
              IF(NTERMS_CONT.GT.1)THEN
                DO K=NTERMS_CONT-1,1,-1
                  POL_SIMUL_YY(NSIM)=POL_SIMUL_YY(NSIM)*XX+A_(K)
                END DO
              END IF
            ELSE
              XX=X0
              POL_SIMUL_YY1(NSIM)=A_(NTERMS_CONT)
              IF(NTERMS_CONT.GT.1)THEN
                DO K=NTERMS_CONT-1,1,-1
                  POL_SIMUL_YY1(NSIM)=POL_SIMUL_YY1(NSIM)*XX+A_(K)
                END DO
              END IF
              XX=X0+DELTAX
              POL_SIMUL_YY2(NSIM)=A_(NTERMS_CONT)
              IF(NTERMS_CONT.GT.1)THEN
                DO K=NTERMS_CONT-1,1,-1
                  POL_SIMUL_YY2(NSIM)=POL_SIMUL_YY2(NSIM)*XX+A_(K)
                END DO
              END IF
            END IF
          END DO
          IF((CTYPEFIT.EQ.'1').OR.(CTYPEFIT.EQ.'4'))THEN
            EYY=0.0
            DO NSIM=1,NITER
              EYY=EYY+(YY-POL_SIMUL_YY(NSIM))*(YY-POL_SIMUL_YY(NSIM))
            END DO
            EYY=SQRT(EYY/REAL(NSIM-1))
          ELSE
            EYY1=0.0
            DO NSIM=1,NITER
              EYY1=EYY1+
     +         (YY1-POL_SIMUL_YY1(NSIM))*(YY1-POL_SIMUL_YY1(NSIM))
            END DO
            EYY1=SQRT(EYY1/REAL(NSIM-1))
            EYY2=0.0
            DO NSIM=1,NITER
              EYY2=EYY2+
     +         (YY2-POL_SIMUL_YY2(NSIM))*(YY2-POL_SIMUL_YY2(NSIM))
            END DO
            EYY2=SQRT(EYY2/REAL(NSIM-1))
          END IF
        END IF
C
        WRITE(*,152)
        WRITE(*,101) '>>> Type of fit'//
     +   '................................: '//CTYPEFIT
        IF(CRCONT.EQ.'y')THEN
          IF((CTYPEFIT.EQ.'1').OR.(CTYPEFIT.EQ.'4'))THEN
            WRITE(*,100) '>>> Continuum flux ..............'//
     +       '[obser. frame]: '
            WRITE(*,*) YY,EYY
            EW=-SQRT(2.*PI)*AMP*SIGMA*DISP/YY
            ERREW=SQRT(2.*PI)*
     +       SQRT(SIGMA*SIGMA*DISP*DISP*EAMP*EAMP/(YY*YY)+
     +       AMP*AMP*DISP*DISP*ESIGMA*ESIGMA/(YY*YY)+
     +       AMP*AMP*SIGMA*SIGMA*DISP*DISP*
     +       EYY*EYY/(YY*YY*YY*YY))
            WRITE(*,100) '>>> EW under gauss.(Angstroms)...'//
     +       '[obser. frame]: '
            WRITE(*,*) EW,ERREW
            WRITE(*,100) '>>> EW under gauss.(Angstroms)...'//
     +       '[obser. frame]: '
            WRITE(*,*) EW,ERREW
            WRITE(*,100) '>>> EW under spectrum (Angstroms)'//
     +       '[obser. frame]: '
            WRITE(*,*) WIDTH,ERRWIDTH
            WRITE(*,100) '>>> EW under gauss.(Angstroms).....'//
     +       '[rest frame]: '
            WRITE(*,*) EW/RCVEL1,ERREW/RCVEL1
            WRITE(*,100) '>>> EW under spectrum (Angstroms)..'//
     +       '[rest frame]: '
            WRITE(*,*) WIDTH/RCVEL1,ERRWIDTH/RCVEL1
          ELSE
            WRITE(*,100) '>>> Continuum flux #1............'//
     +       '[obser. frame]: '
            WRITE(*,*) YY1,EYY1
            EW1=-SQRT(2.*PI)*AMP1*SIGMA*DISP/YY1
            ERREW1=SQRT(2.*PI)*
     +       SQRT(SIGMA*SIGMA*DISP*DISP*EAMP1*EAMP1/(YY1*YY1)+
     +       AMP1*AMP1*DISP*DISP*ESIGMA*ESIGMA/(YY1*YY1)+
     +       AMP1*AMP1*SIGMA*SIGMA*DISP*DISP*
     +       EYY1*EYY1/(YY1*YY1*YY1*YY1))
            WRITE(*,100) '>>> EW under gauss.(Angstroms)#1.'//
     +       '[obser. frame]: '
            WRITE(*,*) EW1,ERREW1
            WRITE(*,100) '>>> EW under spectrum (Angstroms)'//
     +       '[obser. frame]: '
            WRITE(*,*) WIDTH,ERRWIDTH
            WRITE(*,100) '>>> EW under gauss.(Angstroms)#1...'//
     +       '[rest frame]: '
            WRITE(*,*) EW1/RCVEL1,ERREW1/RCVEL1
            WRITE(*,100) '>>> EW under spectrum (Angstroms)#1'//
     +       '[rest frame]: '
            WRITE(*,*) WIDTH/RCVEL1,ERRWIDTH/RCVEL1
            WRITE(*,100) '>>> Continuum flux #2............'//
     +       '[obser. frame]: '
            WRITE(*,*) YY2,EYY2
            EW2=-SQRT(2.*PI)*AMP2*SIGMA*DISP/YY2
            ERREW2=SQRT(2.*PI)*
     +       SQRT(SIGMA*SIGMA*DISP*DISP*EAMP2*EAMP2/(YY2*YY2)+
     +       AMP2*AMP2*DISP*DISP*ESIGMA*ESIGMA/(YY2*YY2)+
     +       AMP2*AMP2*SIGMA*SIGMA*DISP*DISP*
     +       EYY2*EYY2/(YY2*YY2*YY2*YY2))
            WRITE(*,100) '>>> EW under gauss.(Angstroms)#2.'//
     +       '[obser. frame]: '
            WRITE(*,*) EW2,ERREW2
            WRITE(*,100) '>>> EW under spectrum (Angstroms)'//
     +       '[obser. frame]: '
            WRITE(*,*) WIDTH,ERRWIDTH
            WRITE(*,100) '>>> EW under gauss.(Angstroms)#2...'//
     +       '[rest frame]: '
            WRITE(*,*) EW2/RCVEL1,ERREW2/RCVEL1
            WRITE(*,100) '>>> EW under spectrum (Angstroms)#2'//
     +       '[rest frame]: '
            WRITE(*,*) WIDTH/RCVEL1,ERRWIDTH/RCVEL1
          END IF
        END IF
        WRITE(*,100) '>>> Integrated Flux (spec.-cont.)...'//
     +   '[any frame]: '
        WRITE(*,*) FLUX1,ERRFLUX1
        WRITE(*,100) '>>> Integrated Flux (continuum).....'//
     +   '[any frame]: '
        WRITE(*,*) FLUXC,ERRFLUXC
        WRITE(*,100) '>>> Integrated Flux (spectrum)......'//
     +   '[any frame]: '
        WRITE(*,*) FLUX2,ERRFLUX2
        IF((CTYPEFIT.EQ.'1').OR.(CTYPEFIT.EQ.'4'))THEN
          FLUXG=-AMP*SQRT(2.*PI)*SIGMA*DISP
          EFLUXG=-SQRT(2.*PI)*DISP*SQRT(AMP*AMP*ESIGMA*ESIGMA+
     +     SIGMA*SIGMA*EAMP*EAMP)
          WRITE(*,100) '>>> Integrated Flux (gaussian)......'//
     +     '[any frame]: '
          WRITE(*,*) FLUXG,EFLUXG
        ELSE
          FLUXG=-AMP1*SQRT(2.*PI)*SIGMA*DISP
          EFLUXG=-SQRT(2.*PI)*DISP*SQRT(AMP1*AMP1*ESIGMA*ESIGMA+
     +     SIGMA*SIGMA*EAMP1*EAMP1)
          WRITE(*,100) '>>> Integrated Flux (gaussian)#1....'//
     +     '[any frame]: '
          WRITE(*,*) FLUXG,EFLUXG
          FLUXG=-AMP2*SQRT(2.*PI)*SIGMA*DISP
          EFLUXG=-SQRT(2.*PI)*DISP*SQRT(AMP2*AMP2*ESIGMA*ESIGMA+
     +     SIGMA*SIGMA*EAMP2*EAMP2)
          WRITE(*,100) '>>> Integrated Flux (gaussian)#2....'//
     +     '[any frame]: '
          WRITE(*,*) FLUXG,EFLUXG
        END IF
        WRITE(*,*)
        FWHM=SIGMA*DISP*SQRT(-8.*ALOG(0.5))
        WRITE(*,100) 'FWHM (Angstroms): '
        WRITE(*,*) FWHM,SQRT(-8.*ALOG(0.5))*DISP*ESIGMA
        IF((CTYPEFIT.EQ.'1').OR.(CTYPEFIT.EQ.'4'))THEN
          WRITE(*,101) 'Coefficients from fit: '//
     +     'y=a*exp[-((x-x0)**2)/(2 sig**2)]'
          WRITE(*,100) '> a   (+ error, rmsDOWNHILL) [flux]: '
          WRITE(*,*)AMP,EAMP,EEAMP
          WRITE(*,100) '> x0  (+ error, rmsDOWNHILL) [pix.]: '
          WRITE(*,*)X0,EX0,EEX0
          WRITE(*,100) '> sig (+ error, rmsDOWNHILL) [pix.]: '
          WRITE(*,*)SIGMA,ESIGMA,EESIGMA
          WRITE(*,100) '> Central wavelength [obser. frame]: '
          WRITE(*,*)STWV+(X0-1.)*DISP
          WRITE(*,100) '> Central wavelength...[rest frame]: '
          WRITE(*,*)(STWV+(X0-1.)*DISP)/RCVEL1
        ELSE
          WRITE(*,101) 'Coefficients from fit: '//
     +     'y=a1*exp[-((x-x0)**2)/(2 sig**2)] +'
          WRITE(*,101) '                       '//
     +     '+ a2*exp[-((x+deltax-x0)**2)/(2 sig**2)]'
          WRITE(*,100) '> a #1(+ error, rmsDOWNHILL) [flux]: '
          WRITE(*,*)AMP1,EAMP1,EEAMP1
          WRITE(*,100) '> a #2(+ error, rmsDOWNHILL) [flux]: '
          WRITE(*,*)AMP2,EAMP2,EEAMP2
          WRITE(*,100) '> x0#1(+ error, rmsDOWNHILL) [pix.]: '
          WRITE(*,*)X0,EX0,EEX0
          WRITE(*,100) '> x0#2(+ error, rmsDOWNHILL) [pix.]: '
          WRITE(*,*)X0+DELTAX,EX0,EEX0
          WRITE(*,100) '> sig (+ error, rmsDOWNHILL) [pix.]: '
          WRITE(*,*)SIGMA,ESIGMA,EESIGMA
          WRITE(*,100) '> Central wavelength#1[obser.frame]: '
          WRITE(*,*)STWV+(X0-1.)*DISP
          WRITE(*,100) '> Central wavelength#2[obser.frame]: '
          WRITE(*,*)STWV+(X0+DELTAX-1.)*DISP
          WRITE(*,100) '> Central wavelength#1.[rest frame]: '
          WRITE(*,*)(STWV+(X0-1.)*DISP)/RCVEL1
          WRITE(*,100) '> Central wavelength#2.[rest frame]: '
          WRITE(*,*)(STWV+(X0+DELTAX-1.)*DISP)/RCVEL1
        END IF
        WRITE(*,152)
C
        IF(COUT.EQ.'y')THEN
          WRITE(30,152)
          WRITE(30,101) '>>> Type of fit'//
     +     '................................: '//CTYPEFIT
          IF(CRCONT.EQ.'y')THEN
            IF((CTYPEFIT.EQ.'1').OR.(CTYPEFIT.EQ.'4'))THEN
              WRITE(30,100) '>>> Continuum flux ..............'//
     +         '[obser. frame]: '
              WRITE(30,*) YY,EYY
              WRITE(30,100) '>>> EW under gauss.(Angstroms)...'//
     +         '[obser. frame]: '
              WRITE(30,*) EW,ERREW
ccc cotas inferiores en EW usando como continuo el valor medido del continuo
ccc mas una vez el r.m.s.
ccc            write(30,*) NCHAN_CONTINUUM,-SQRT(2.*PI)*AMP*SIGMA*DISP/
ccc     +       (YY+SQRT(REAL(NCHAN_CONTINUUM))*EYY)
              WRITE(30,100) '>>> EW under spectrum (Angstroms)'//
     +         '[obser. frame]: '
              WRITE(30,*) WIDTH,ERRWIDTH
              WRITE(30,100) '>>> EW under gauss.(Angstroms).....'//
     +         '[rest frame]: '
              WRITE(30,*) EW/RCVEL1,ERREW/RCVEL1
ccc cotas inferiores en EW usando como continuo el valor medido del continuo
ccc mas una vez el r.m.s.
ccc            write(30,*) NCHAN_CONTINUUM,-SQRT(2.*PI)*AMP*SIGMA*DISP/
ccc     +       (YY+SQRT(REAL(NCHAN_CONTINUUM))*EYY)/RCVEL1
              WRITE(30,100) '>>> EW under spectrum (Angstroms)..'//
     +         '[rest frame]: '
              WRITE(30,*) WIDTH/RCVEL1,ERRWIDTH/RCVEL1
            ELSE
              WRITE(30,100) '>>> Continuum flux #1............'//
     +         '[obser. frame]: '
              WRITE(30,*) YY1,EYY1
              WRITE(30,100) '>>> EW under gauss.(Angstroms)...'//
     +         '[obser. frame]: '
              WRITE(30,*) EW1,ERREW1
              WRITE(30,100) '>>> EW under spectrum (Angstroms)'//
     +         '[obser. frame]: '
              WRITE(30,*) WIDTH,ERRWIDTH
              WRITE(30,100) '>>> EW under gauss.(Angstroms).....'//
     +         '[rest frame]: '
              WRITE(30,*) EW1/RCVEL1,ERREW1/RCVEL1
              WRITE(30,100) '>>> EW under spectrum (Angstroms)..'//
     +         '[rest frame]: '
              WRITE(30,*) WIDTH/RCVEL1,ERRWIDTH/RCVEL1
              WRITE(30,100) '>>> Continuum flux #2............'//
     +         '[obser. frame]: '
              WRITE(30,*) YY2,EYY2
              WRITE(30,100) '>>> EW under gauss.(Angstroms)...'//
     +         '[obser. frame]: '
              WRITE(30,*) EW2,ERREW2
              WRITE(30,100) '>>> EW under spectrum (Angstroms)'//
     +         '[obser. frame]: '
              WRITE(30,*) WIDTH,ERRWIDTH
              WRITE(30,100) '>>> EW under gauss.(Angstroms).....'//
     +         '[rest frame]: '
              WRITE(30,*) EW2/RCVEL1,ERREW2/RCVEL1
              WRITE(30,100) '>>> EW under spectrum (Angstroms)..'//
     +         '[rest frame]: '
              WRITE(30,*) WIDTH/RCVEL1,ERRWIDTH/RCVEL1
            END IF
          END IF
          WRITE(30,100) '>>> Integrated Flux (spec.-cont.)...'//
     +     '[any frame]: '
          WRITE(30,*) FLUX1,ERRFLUX1
          WRITE(30,100) '>>> Integrated Flux (continuum).....'//
     +     '[any frame]: '
          WRITE(30,*) FLUXC,ERRFLUXC
          WRITE(30,100) '>>> Integrated Flux (spectrum)......'//
     +     '[any frame]: '
          WRITE(30,*) FLUX2,ERRFLUX2
          IF((CTYPEFIT.EQ.'1').OR.(CTYPEFIT.EQ.'4'))THEN
            FLUXG=-AMP*SQRT(2.*PI)*SIGMA*DISP
            EFLUXG=-SQRT(2.*PI)*DISP*SQRT(AMP*AMP*ESIGMA*ESIGMA+
     +       SIGMA*SIGMA*EAMP*EAMP)
            WRITE(30,100) '>>> Integrated Flux (gaussian)......'//
     +       '[any frame]: '
            WRITE(30,*) FLUXG,EFLUXG
          ELSE
            FLUXG=-AMP1*SQRT(2.*PI)*SIGMA*DISP
            EFLUXG=-SQRT(2.*PI)*DISP*SQRT(AMP1*AMP1*ESIGMA*ESIGMA+
     +       SIGMA*SIGMA*EAMP1*EAMP1)
            WRITE(30,100) '>>> Integrated Flux (gaussian)#1....'//
     +       '[any frame]: '
            WRITE(30,*) FLUXG,EFLUXG
            FLUXG=-AMP2*SQRT(2.*PI)*SIGMA*DISP
            EFLUXG=-SQRT(2.*PI)*DISP*SQRT(AMP2*AMP2*ESIGMA*ESIGMA+
     +       SIGMA*SIGMA*EAMP2*EAMP2)
            WRITE(30,100) '>>> Integrated Flux (gaussian)#2....'//
     +       '[any frame]: '
            WRITE(30,*) FLUXG,EFLUXG
          END IF
          WRITE(30,*)
          FWHM=SIGMA*DISP*SQRT(-8.*ALOG(0.5))
          WRITE(30,100) 'FWHM (Angstroms): '
          WRITE(30,*) FWHM,SQRT(-8.*ALOG(0.5))*DISP*ESIGMA
          IF((CTYPEFIT.EQ.'1').OR.(CTYPEFIT.EQ.'4'))THEN
            WRITE(30,101) 'Coefficients from fit: '//
     +       'y=a*exp[-((x-x0)**2)/(2 sig**2)]'
            WRITE(30,100) '> a   (+ error, rmsDOWNHILL) [flux]: '
            WRITE(30,*)AMP,EAMP,EEAMP
            WRITE(30,100) '> x0  (+ error, rmsDOWNHILL) [pix.]: '
            WRITE(30,*)X0,EX0,EEX0
            WRITE(30,100) '> sig (+ error, rmsDOWNHILL) [pix.]: '
            WRITE(30,*)SIGMA,ESIGMA,EESIGMA
            WRITE(30,100) '> Central wavelength [obser. frame]: '
            WRITE(30,*)STWV+(X0-1.)*DISP
            WRITE(30,100) '> Central wavelength...[rest frame]: '
            WRITE(30,*)(STWV+(X0-1.)*DISP)/RCVEL1
          ELSE
            WRITE(30,101) 'Coefficients from fit: '//
     +       'y=a1*exp[-((x-x0)**2)/(2 sig**2)] +'
            WRITE(30,101) '                       '//
     +       '+ a2*exp[-((x+deltax-x0)**2)/(2 sig**2)]'
            WRITE(30,100) '> a #1(+ error, rmsDOWNHILL) [flux]: '
            WRITE(30,*)AMP1,EAMP1,EEAMP1
            WRITE(30,100) '> a #2(+ error, rmsDOWNHILL) [flux]: '
            WRITE(30,*)AMP2,EAMP2,EEAMP2
            WRITE(30,100) '> x0#1(+ error, rmsDOWNHILL) [pix.]: '
            WRITE(30,*)X0,EX0,EEX0
            WRITE(30,100) '> x0#2(+ error, rmsDOWNHILL) [pix.]: '
            WRITE(30,*)X0+DELTAX,EX0,EEX0
            WRITE(30,100) '> sig (+ error, rmsDOWNHILL) [pix.]: '
            WRITE(30,*)SIGMA,ESIGMA,EESIGMA
            WRITE(30,100) '> Central wavelength#1[obser.frame]: '
            WRITE(30,*)STWV+(X0-1.)*DISP
            WRITE(30,100) '> Central wavelength#2[obser.frame]: '
            WRITE(30,*)STWV+(X0+DELTAX-1.)*DISP
            WRITE(30,100) '> Central wavelength#1.[rest frame]: '
            WRITE(30,*)(STWV+(X0-1.)*DISP)/RCVEL1
            WRITE(30,100) '> Central wavelength#2.[rest frame]: '
            WRITE(30,*)(STWV+(X0+DELTAX-1.)*DISP)/RCVEL1
          END IF
          WRITE(30,152)
        END IF
C
        NPLOT=0
        DO J=JMIN,JMAX
          NPLOT=NPLOT+1
          XPLOT(NPLOT)=REAL(J)
          POL(J)=A(NTERMS_CONT)
          DO K=NTERMS_CONT-1,1,-1
            POL(J)=POL(J)*XPLOT(NPLOT)+A(K)
          END DO
        END DO
        IF((CTYPEFIT.EQ.'1').OR.(CTYPEFIT.EQ.'4'))THEN
          NPLOT=0
          DO J=JMIN,JMAX
            NPLOT=NPLOT+1
            YPLOT(NPLOT)=POL(J)+AMP*EXP(-(XPLOT(NPLOT)-X0)*
     +       (XPLOT(NPLOT)-X0)/(2.*SIGMA*SIGMA))
          END DO
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(6)
            CALL PGLINE(NPLOT,XPLOT,YPLOT)
              IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END DO
        ELSE
          NPLOT=0
          DO J=JMIN,JMAX
            NPLOT=NPLOT+1
            YPLOT(NPLOT)=POL(J)+AMP1*EXP(-(XPLOT(NPLOT)-X0)*
     +       (XPLOT(NPLOT)-X0)/(2.*SIGMA*SIGMA))
          END DO
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(8)
            CALL PGLINE(NPLOT,XPLOT,YPLOT)
              IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END DO
          NPLOT=0
          DO J=JMIN,JMAX
            NPLOT=NPLOT+1
            YPLOT(NPLOT)=POL(J)+AMP2*EXP(-(XPLOT(NPLOT)-X0-DELTAX)*
     +       (XPLOT(NPLOT)-X0-DELTAX)/(2.*SIGMA*SIGMA))
          END DO
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(8)
            CALL PGLINE(NPLOT,XPLOT,YPLOT)
              IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END DO
          NPLOT=0
          DO J=JMIN,JMAX
            NPLOT=NPLOT+1
            YPLOT(NPLOT)=POL(J)+AMP1*EXP(-(XPLOT(NPLOT)-X0)*
     +       (XPLOT(NPLOT)-X0)/(2.*SIGMA*SIGMA))+
     +       AMP2*EXP(-(XPLOT(NPLOT)-X0-DELTAX)*
     +       (XPLOT(NPLOT)-X0-DELTAX)/(2.*SIGMA*SIGMA))
          END DO
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(6)
            CALL PGLINE(NPLOT,XPLOT,YPLOT)
              IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END DO
        END IF
C
100     FORMAT(A,$)
101     FORMAT(A)
152     FORMAT(79('*'))
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
