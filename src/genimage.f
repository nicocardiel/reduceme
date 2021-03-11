C------------------------------------------------------------------------------
C Version 28-April-2015                                        file: genimage.f
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
C Program: genimage
C Classification: arithmetic & manipulations
C Description: Creates an artificial image using, if necessary, regions from
C other images and/or tabulated data from ASCII files.
C
Comment
C
C Crea una imagen sintetica a gusto del usuario
C
        PROGRAM GENIMAGE
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER TRUELEN
        INTEGER READI
        INTEGER READILIM
        REAL READF
C
        INTEGER NMAXCOLUMN
        PARAMETER (NMAXCOLUMN=100)        !numero maximo de columnas a leer en 
                                          !la lectura de espectros tabulados en 
                                          !ficheros ASCII
        INTEGER NMAXROWTABLE
        PARAMETER (NMAXROWTABLE=1000000) !número máximo de filas que pueden
                                         !leerse de un fichero ASCII
        INTEGER NMAXTERM
        PARAMETER (NMAXTERM=20) !numero maximo de terminos en polinomio
C
        REAL PI
        PARAMETER (PI=3.141592654)
        REAL PI2
        PARAMETER (PI2=6.283185307)
        REAL SQRT2
        PARAMETER (SQRT2=1.414213562)
        REAL C
        PARAMETER (C=299792.458)
C
        REAL LININTERP
        REAL INTEGTAB
        REAL FINTGAUSS
        REAL RANRED
C
        INTEGER I,J,L
        INTEGER J1,J2
        INTEGER K,KTOT,N1,N2
        INTEGER IFLAG,IFLAG1,IFLAG2
        INTEGER II,JJ
        INTEGER NS1,NS2,NS0
        INTEGER NC1,NC2,NC0
        INTEGER NSCAN2,NCHAN2
        INTEGER NP
        INTEGER NSMAX_LOCAL,NCMAX_LOCAL
        INTEGER NSKIP,NCOLUMN,NCOLUMNW,NCOLUMNF
        INTEGER NTERM
        INTEGER I0
        INTEGER NLINES,NLINEST
        INTEGER NEXTINFO
        INTEGER NSEED
        REAL A(NCMAX,NSMAX)
        REAL B(NCMAX,NSMAX)
        REAL S(NCMAX)                     !dimensionado al mayor de NCMAX/NSMAX
        REAL S_SPL(NCMAX),A_SPL(NCMAX),B_SPL(NCMAX),C_SPL(NCMAX)
        REAL NEWVAL
        REAL X(NCMAX),Y(NCMAX),XP,YP
        REAL XTABLE(NMAXROWTABLE),YTABLE(NMAXROWTABLE)
        REAL FSTWV,FDISP
        REAL FACTOR_WAVE,FACTOR_FLUX
        REAL RVEL,RCVEL,RCVEL1
        REAL STWV2,DISP2
        REAL AIRMASS0,TIMEXPOS0,STWV0,DISP0
        REAL DUMCOLUMN(NMAXCOLUMN)
        REAL COEF(NMAXTERM)
        REAL FACTOR1
        REAL G_AMP,G_SIGMA,G_X0
        REAL G_TOTAL_FLUX,G_SIGMA_KMS,G_X0_LDO,TG_SIGMA,LDO_CENTRAL
        REAL GX0,GWIDTH,GFACTOR
        REAL WL,WL1,WL2
        REAL R1,R2,STDERR
        CHARACTER*1 COPC,CINTER,CFILE,CXTYPE,CCHANGE,CCC,CEXISTING
        CHARACTER*1 CCONT,CSURE
        CHARACTER*50 CDUMMY
        CHARACTER*255 INFILE,DATAFILE,OUTFILE,DUMLINE
        CHARACTER*255 FILELINES
        CHARACTER*255 OBJECT0,FITSFILE0,COMMENT0
        LOGICAL LOUTLIMITS
        LOGICAL LOGFILE
C------------------------------------------------------------------------------
        THISPROGRAM='genimage'
        CALL WELCOME('08-March-2005')
C------------------------------------------------------------------------------
C protecciones
        NSMAX_LOCAL=NSMAX
        NCMAX_LOCAL=NCMAX
        !las matrices están dimensionadas típicamente a NCMAX
        !independientemente de si se utilizan para introducir espectros o 
        !cortes en la dirección espacial:
        IF(NSMAX_LOCAL.GT.NCMAX_LOCAL)STOP 'FATAL ERROR: NSMAX.GT.NCMAX'
        NSEED=-1
C------------------------------------------------------------------------------
        WRITE(*,100)'Initialize frame to existing file (y/n) '
        CEXISTING(1:1)=READC('n','yn')
        IF(CEXISTING.EQ.'n')THEN
          STWV0=0.
          DISP0=0.
          AIRMASS0=0.
          TIMEXPOS0=0.
          OBJECT0='[from genimage]'
          FITSFILE0='[none]'
          COMMENT0='[from genimage]'
          WRITE(*,101)'Enter image dimensions:'
          WRITE(*,100)'NSCAN '
          NSCAN=READILIM('@',1,NSMAX)
          WRITE(*,100)'NCHAN '
          NCHAN=READILIM('@',1,NCMAX)
C inicializamos imagen de salida a cero
          WRITE(*,100)'Initializing whole image to zero...'
          DO I=1,NSCAN
            DO J=1,NCHAN
              A(J,I)=0.
            END DO
          END DO
          WRITE(*,101)'...OK!'
        ELSE
          WRITE(*,100)'Existing file name'
          INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
          DO I=1,NSCAN
            READ(20) (A(J,I),J=1,NCHAN)
          END DO
          CLOSE(20)
          STWV0=STWV
          DISP0=DISP
          AIRMASS0=AIRMASS
          TIMEXPOS0=TIMEXPOS
          OBJECT0=OBJECT
          FITSFILE0=FITSFILE
          COMMENT0=COMMENT
        END IF
C------------------------------------------------------------------------------
C menu principal del programa
15      WRITE(*,*)
        WRITE(*,100) '>>> NSCAN: '
        WRITE(*,*) NSCAN
        WRITE(*,100) '>>> NCHAN: '
        WRITE(*,*) NCHAN
        WRITE(*,100) '>>> STWV: '
        WRITE(*,*) STWV0
        WRITE(*,100) '>>> DISP: '
        WRITE(*,*) DISP0
        WRITE(*,*)
        WRITE(*,101)'1 - Replace spectra by spectrum file '//
     +   '(with same NCHAN)'
        WRITE(*,101)'2 - Replace spatial directions by spatial '//
     +   'direction from file (with same NSCAN)'
        WRITE(*,101)'3 - Replace image regions by constants'
        WRITE(*,101)'4 - Replace spectra by fit to splines'
        WRITE(*,101)'5 - Replace spectra by tabulated spectrum'
        WRITE(*,101)'6 - Replace WHOLE image by tabulated data'
        WRITE(*,101)'7 - Replace image region by other image region '//
     +   '(with arbitrary dimensions)'
        WRITE(*,101)'8 - Replace spectra by polynomial'
        WRITE(*,101)'a - Add Gaussian noise (with constant STD) '//
     +   'to image region'
        WRITE(*,101)'b - Add Gaussian noise (with STD given by '//
     +   'another image with same dimensions)'
        WRITE(*,101)'m - Multiply spatial profile by a Gaussian'
        WRITE(*,101)'g - Add gaussian'
        WRITE(*,101)'h - Add gaussians from line list in file'
        WRITE(*,101)'0 - EXIT (save image)'
        WRITE(*,101)'q - QUIT (do NOT save image)'
        WRITE(*,*)
        WRITE(*,100)'Option'
        COPC(1:1)=READC('@','qQ012345678aAbBmMgGhH')
        IF((COPC.EQ.'q').OR.(COPC.EQ.'Q'))THEN
          WRITE(*,100)'Are you sure (y/n) '
          CSURE(1:1)=READC('n','yn')
          IF(CSURE.EQ.'y') STOP
          GOTO 15
        ELSEIF(COPC.EQ.'0')THEN
          GOTO 900
        ELSEIF(COPC.EQ.'1')THEN
          GOTO 20
        ELSEIF(COPC.EQ.'2')THEN
          GOTO 30
        ELSEIF(COPC.EQ.'3')THEN
          GOTO 40
        ELSEIF(COPC.EQ.'4')THEN
          GOTO 50
        ELSEIF(COPC.EQ.'5')THEN
          GOTO 60
        ELSEIF(COPC.EQ.'6')THEN
          GOTO 70
        ELSEIF(COPC.EQ.'7')THEN
          GOTO 80
        ELSEIF(COPC.EQ.'8')THEN
          GOTO 85
        ELSEIF((COPC.EQ.'a').OR.(COPC.EQ.'A'))THEN
          GOTO 300
        ELSEIF((COPC.EQ.'b').OR.(COPC.EQ.'B'))THEN
          GOTO 310
        ELSEIF((COPC.EQ.'m').OR.(COPC.EQ.'M'))THEN
          GOTO 87
        ELSEIF((COPC.EQ.'g').OR.(COPC.EQ.'G'))THEN
          GOTO 90
        ELSEIF((COPC.EQ.'h').OR.(COPC.EQ.'H'))THEN
          GOTO 94
        END IF
C------------------------------------------------------------------------------
C introducimos forma espectral en los scans solicitados
20      WRITE(*,*)
        WRITE(*,100)'Spectrum file name'
        INFILE=INFILEX(20,'@',NSCAN2,NCHAN2,STWV2,DISP2,1,.FALSE.)
        IF(NCHAN2.NE.NCHAN)THEN
          WRITE(*,101)'ERROR: NCHAN in last image is different.'
          CLOSE(20)
          GOTO 15
        END IF
        IF(NSCAN2.NE.1)THEN
          WRITE(*,101)'ERROR: NSCAN is greater than 1.'
          CLOSE(20)
          GOTO 15
        END IF
        IF(STWV2.NE.STWV0)THEN
          WRITE(*,101)'WARNING: STWV in last image is different.'
          WRITE(*,101)'Do you want to continue anyway (y/n) '
          CCONT(1:1)=READC('n','yn')
          IF(CCONT.EQ.'n') GOTO 15
        END IF
        IF(DISP2.NE.DISP0)THEN
          WRITE(*,101)'WARNING: DISP in last image is different.'
          WRITE(*,101)'Do you want to continue anyway (y/n) '
          CCONT(1:1)=READC('n','yn')
          IF(CCONT.EQ.'n') GOTO 15
        END IF
        READ(20) (S(J),J=1,NCHAN)
        CLOSE(20)
22      WRITE(*,100)'Enter scan region to introduce spectrum (0,0=EXIT)'
        CALL READ2I('0,0',NS1,NS2)
        IF((NS1.EQ.0).AND.(NS2.EQ.0)) GOTO 15
        IF((NS1.LT.1).OR.(NS2.GT.NSCAN).OR.(NS1.GT.NS2))THEN
          WRITE(*,101)'ERROR: invalid entry. Try again.'
          GOTO 22
        END IF
        DO I=NS1,NS2
          DO J=1,NCHAN
            A(J,I)=S(J)
          END DO
        END DO
        GOTO 22
C------------------------------------------------------------------------------
C introducimos forma espacial en los canales solicitados
30      WRITE(*,*)
        WRITE(*,100)'Spatial form file name'
        INFILE=INFILEX(20,'@',NSCAN2,NCHAN2,STWV2,DISP2,1,.FALSE.)
        IF(NSCAN2.NE.NSCAN)THEN
          WRITE(*,101)'ERROR: NSCAN in last image is different.'
          CLOSE(20)
          GOTO 15
        END IF
        IF(NCHAN2.NE.1)THEN
          WRITE(*,101)'ERROR: NCHAN is greater than 1.'
          CLOSE(20)
          GOTO 15
        END IF
        IF(STWV2.NE.STWV0)THEN
          WRITE(*,101)'WARNING: STWV in last image is different.'
          WRITE(*,101)'Do you want to continue anyway (y/n) '
          CCONT(1:1)=READC('n','yn')
          IF(CCONT.EQ.'n') GOTO 15
        END IF
        IF(DISP2.NE.DISP0)THEN
          WRITE(*,101)'WARNING: DISP in last image is different.'
          WRITE(*,101)'Do you want to continue anyway (y/n) '
          CCONT(1:1)=READC('n','yn')
          IF(CCONT.EQ.'n') GOTO 15
        END IF
        DO I=1,NSCAN
          READ(20) S(I)
        END DO
        CLOSE(20)  
32      WRITE(*,100)'Enter channel region to introduce spatial '//
     +   'direction form (0,0=EXIT) '
        CALL READ2I('0,0',NC1,NC2)
        IF((NC1.EQ.0).AND.(NC2.EQ.0)) GOTO 15
        IF((NC1.LT.1).OR.(NC2.GT.NCHAN).OR.(NC1.GT.NC2))THEN
          WRITE(*,101)'ERROR: invalid entry. Try again.'
          GOTO 32
        END IF
        DO I=1,NSCAN
          DO J=NC1,NC2
            A(J,I)=S(I)
          END DO
        END DO
        GOTO 32
C------------------------------------------------------------------------------
C introducimos constante en las regiones solicitadas
40      WRITE(*,*)
        WRITE(*,101)'Enter region and pixel value:'
41      WRITE(*,100)'Scan region (0,0=EXIT)...'
        CALL READ2I('@',NS1,NS2)
        IF((NS1.EQ.0).AND.(NS2.EQ.0)) GOTO 15
        IF((NS1.LT.1).OR.(NS2.GT.NSCAN).OR.(NS1.GT.NS2))THEN
          WRITE(*,101)'ERROR: invalid entry. Try again.'
          GOTO 41
        END IF
42      WRITE(*,100)'Channel region (0,0=EXIT)'
        CALL READ2I('@',NC1,NC2)
        IF((NC1.EQ.0).AND.(NC2.EQ.0)) GOTO 15
        IF((NC1.LT.1).OR.(NC2.GT.NCHAN).OR.(NC1.GT.NC2))THEN
          WRITE(*,101)'ERROR: invalid entry. Try again.'
          GOTO 42
        END IF
        WRITE(*,100)'Pixel(s) value'
        NEWVAL=READF('@')
        DO I=NS1,NS2
          DO J=NC1,NC2
            A(J,I)=NEWVAL
          END DO
        END DO
        GOTO 40
C------------------------------------------------------------------------------
C introducimos ajuste a splines en los scans solicitados
50      WRITE(*,*)
        WRITE(*,100)'Data through [k]eyboard or [f]ile.........'//
     +   '.......(k/f) '
        CFILE(1:1)=READC('k','kf')
        WRITE(*,100)'X-coordinate in [w]avelength or [c]hannel '//
     +   'number (w/c) '
        CXTYPE(1:1)=READC('w','wc')
C
        NP=0
        IF(CFILE.EQ.'f')THEN
          WRITE(*,100)'Data file name'
          DATAFILE=INFILEX(25,'@',0,0,.0,.0,3,.FALSE.)           !fichero ASCII
51        READ(25,*,END=52)X(NP+1),Y(NP+1)
          NP=NP+1
          IF(NP.GT.100)THEN
            WRITE(*,101)'ERROR: no. of points too large.'
            CLOSE(25)
            GOTO 15
          END IF
          GOTO 51
52        CLOSE(25)
          WRITE(*,110)'No. of points read: ',NP
        ELSE
ccc53        WRITE(*,100)'No. of points (max. 100) '
          WRITE(*,100)'No. of points (max. 100) '
          NP=READILIM('@',2,100)
          WRITE(*,*)
          DO I=1,NP
            WRITE(*,'(A2,I3.3,A3,$)')'X(',I,')'
            X(I)=READF('@')
            WRITE(*,'(A2,I3.3,A3,$)')'Y(',I,')'
            Y(I)=READF('@')
            WRITE(*,*)
          END DO
        END IF
C
        IF(CXTYPE.EQ.'w')THEN
          WRITE(*,100)'STWV'
          FSTWV=READF('@')
          WRITE(*,100)'DISP'
          FDISP=READF('@')
          STWV=FSTWV
          DISP=FDISP
        ELSE
          FSTWV=1.
          FDISP=1.
          STWV=0.
          DISP=0.
        END IF
C
        WRITE(*,100) 'Fitting polynomial...'
        IF(CXTYPE.EQ.'w')THEN
          CALL CUBSPL(X,Y,NP,1,S_SPL,A_SPL,B_SPL,C_SPL)
          I0=1
          DO J=1,NCHAN
            XP=REAL(J-1)*FDISP+FSTWV
            CALL CUBSPLX(X,Y,A_SPL,B_SPL,C_SPL,NP,I0,XP,YP)
            S(J)=YP
          END DO
        ELSE
          CALL CUBSPL(X,Y,NP,1,S_SPL,A_SPL,B_SPL,C_SPL)
          I0=1
          DO J=1,NCHAN
            XP=REAL(J)
            CALL CUBSPLX(X,Y,A_SPL,B_SPL,C_SPL,NP,I0,XP,YP)
            S(J)=YP
          END DO
        END IF
        WRITE(*,101) '  ...OK!'
C
        GOTO 22
C------------------------------------------------------------------------------
C introducimos tabla en las regiones solicitadas
60      WRITE(*,*)
        IF(NSCAN.EQ.1)THEN
          NS1=1
          NS2=1
          GOTO 62
        END IF
        WRITE(*,101)'Enter region:'
61      WRITE(*,100)'Scan region (0,0=EXIT)...'
        CALL READ2I('@',NS1,NS2)
        IF((NS1.EQ.0).AND.(NS2.EQ.0)) GOTO 15
        IF((NS1.LT.1).OR.(NS2.GT.NSCAN).OR.(NS1.GT.NS2))THEN
          WRITE(*,101)'ERROR: invalid entry. Try again.'
          GOTO 61
        END IF
C En la opcion 1 sólo se lee una columna, que debe contener el flujo en cada
C pixel. Se asume que el incremento en l.d.o. es constante y conocido.
C En las opciones 2 y 3 se leen dos columnas, una es la longitud de onda y
C otra el flujo. La diferencia entre estas dos opciones estriba en que en la
C opcion 2 se sustituye el flujo en el pixel por el resultado de la
C interpolación lineal en la tabla X,Y, mientras que en la opción 3 se calcula
C la integral de la función tabulada para los límites en l.d.o. de cada pixel.
62      WRITE(*,101) 'Select type of ASCII file:'
        WRITE(*,101) '1) Table n,Y: flux in a single column, FIXED '//
     +   'increment in wavelength assumed'
        WRITE(*,101) '2) Table X,Y: X=wavelength, Y=flux, linear '//
     +   'interpolation'
        WRITE(*,101) '3) Table X,Y: X=wavelength, Y=flux, computing '//
     +   'area'
        WRITE(*,101) '0) EXIT'
        WRITE(*,100) 'Option (0/1/2/3) '
        CINTER(1:1)=READC('0','0123')
        IF(CINTER.EQ.'0') GOTO 15
C
        WRITE(*,100)'File with tabulated spectrum'
        DATAFILE=INFILEX(25,'@',0,0,.0,.0,3,.FALSE.)  !abrimos el fichero ASCII
C
        IF(CINTER.EQ.'1')THEN !.......................................Table n,Y
          WRITE(*,100)'No. of lines to be skipped '
          NSKIP=READI('0')
          IF(NSKIP.GT.0)THEN
            DO J=1,NSKIP
              READ(25,*)
            END DO
          END IF
          WRITE(*,100)'Column number where spectrum is located '
          NCOLUMN=READILIM('1',1,NMAXCOLUMN)
          WRITE(*,100)'Factor to be applied to the flux values '
          FACTOR_FLUX=READF('1.0')
          DO J=1,NCHAN
            READ(25,*)(DUMCOLUMN(I),I=1,NCOLUMN)
            S(J)=DUMCOLUMN(NCOLUMN)
            IF(J.EQ.1)THEN
              WRITE(*,100)'>>> FIRST LINE: ... '
              WRITE(DUMLINE,*) S(J)
              WRITE(*,101)DUMLINE(1:TRUELEN(DUMLINE))
            END IF
            IF(J.EQ.NCHAN)THEN
              WRITE(*,100)'>>> LAST  LINE: ... '
              WRITE(DUMLINE,*) S(J)
              WRITE(*,101)DUMLINE(1:TRUELEN(DUMLINE))
            END IF
          END DO
          CLOSE(25)
          DO I=NS1,NS2
            DO J=1,NCHAN
              A(J,I)=S(J)*FACTOR_FLUX
            END DO
          END DO
          GOTO 15
        END IF
C seguimos con las opciones 2 y 3
        WRITE(*,101) 'Please, confirm current wavelength calibration:'
        WRITE(CDUMMY,*) STWV0
        WRITE(*,100) 'STWV '
        STWV0=READF(CDUMMY)
        WRITE(CDUMMY,*) DISP0
        WRITE(*,100) 'DISP '
        DISP0=READF(CDUMMY)
        WRITE(*,100)'Column number where wavelength is located '
        NCOLUMNW=READILIM('1',1,NMAXCOLUMN)
        WRITE(*,100)'Column number where flux is located...... '
        NCOLUMNF=READILIM('2',1,NMAXCOLUMN)
        WRITE(*,100)'Factor to be applied to the wavelength values '
        FACTOR_WAVE=READF('1.0')
        WRITE(*,100)'Radial velocity (km/sec) '
        RVEL=READF('0.0')
        RCVEL=RVEL/C                                     !z
        RCVEL1=1.+RCVEL                                  !(1+z)
        RCVEL1=RCVEL1/SQRT(1.-RCVEL*RCVEL)               !correcion relativista
        WRITE(*,100)'Factor to be applied to the flux values '
        FACTOR_FLUX=READF('1.0')
C leemos la tabla completa y cerramos el fichero ASCII
        NCOLUMN=MAX(NCOLUMNW,NCOLUMNF)
        K=0
67      READ(25,*,END=68) (DUMCOLUMN(I),I=1,NCOLUMN)
        IF(K+1.GT.NMAXROWTABLE)THEN
          CLOSE(25)
          WRITE(*,101) 'FATAL ERROR: ASCII table is too large!'
          STOP
        END IF
        XTABLE(K+1)=DUMCOLUMN(NCOLUMNW)*FACTOR_WAVE*RCVEL1
        YTABLE(K+1)=DUMCOLUMN(NCOLUMNF)*FACTOR_FLUX*RCVEL1
        K=K+1
        GOTO 67
68      CLOSE(25)
        KTOT=K
C evaluamos el valor de cada píxel en función de la l.d.o.
        IF(CINTER.EQ.'2')THEN !.................Table X,Y, linear interpolation
          DO J=1,NCHAN
            WL=REAL(DBLE(STWV0)+DBLE(J-1)*DBLE(DISP0))
            S(J)=LININTERP(KTOT,XTABLE,YTABLE,WL,IFLAG,N1,N2)
            IF(IFLAG.NE.0)THEN
              WRITE(*,100) 'IFLAG='
              WRITE(*,*) IFLAG
              WRITE(*,101) 'FATAL ERROR: while executing LININTERP.'
              STOP
            END IF
          END DO
        ELSE !........................................Table X,Y, computing area
          DO J=1,NCHAN
            WL=REAL(DBLE(STWV0)+DBLE(J-1)*DBLE(DISP0))
            WL1=WL-DISP0/2.0
            WL2=WL+DISP0/2.0
            S(J)=INTEGTAB(KTOT,XTABLE,YTABLE,WL1,WL2,IFLAG1,IFLAG2)
!           IF((IFLAG1.NE.0).OR.(IFLAG2.NE.0))THEN
!             WRITE(*,100) 'IFLAG1,IFLAG2='
!             WRITE(*,*) IFLAG1,IFLAG2
!             WRITE(*,101) 'FATAL ERROR: while executing INTEGTAB.'
!             STOP
!           END IF
          END DO
        END IF
C
        DO I=NS1,NS2
          DO J=1,NCHAN
            A(J,I)=S(J)
          END DO
        END DO
        GOTO 15
C------------------------------------------------------------------------------
70      WRITE(*,*)
        WRITE(*,101)'* WARNING: File with tabulated data must contain'
        WRITE(*,101)'           1 column with NSCAN x NCHAN elements'
        WRITE(*,*)
        WRITE(*,100)'Do you want to continue (y/n) '
        CCC(1:1)=READC('y','yn')
        IF(CCC.EQ.'y')THEN
          WRITE(*,100)'File with tabulated data'
          DATAFILE=INFILEX(25,'@',0,0,.0,.0,3,.FALSE.)           !fichero ASCII
          WRITE(*,100)'Reading scan #00000'
          DO I=1,NSCAN
            WRITE(*,'(A,I5,$)')'\b\b\b\b\b',I
            DO J=1,NCHAN
              READ(25,*) A(J,I)
            END DO
          END DO
          CLOSE(25)
          WRITE(*,*)
        END IF
C
        GOTO 15
C------------------------------------------------------------------------------
80      WRITE(*,*)
        WRITE(*,100)'File name'
        INFILE=INFILEX(20,'@',NSCAN2,NCHAN2,STWV2,DISP2,1,.FALSE.)
        IF((NSCAN2.NE.NSCAN).OR.(NCHAN2.NE.NCHAN))THEN
          WRITE(*,101)'WARNING: dimensions in last image are different.'
          WRITE(*,100)'Do you want to continue anyway (y/n) '
          CCONT(1:1)=READC('n','yn')
          IF(CCONT.EQ.'n') GOTO 15
        END IF
        IF(STWV2.NE.STWV0)THEN
          WRITE(*,101)'WARNING: STWV in last image is different.'
          WRITE(*,101)'Do you want to continue anyway (y/n) '
          CCONT(1:1)=READC('n','yn')
          IF(CCONT.EQ.'n') GOTO 15
        END IF
        IF(DISP2.NE.DISP0)THEN
          WRITE(*,101)'WARNING: DISP in last image is different.'
          WRITE(*,101)'Do you want to continue anyway (y/n) '
          CCONT(1:1)=READC('n','yn')
          IF(CCONT.EQ.'n') GOTO 15
        END IF
        DO I=1,NSCAN2
          READ(20) (B(J,I),J=1,NCHAN2)
        END DO
        CLOSE(20)
83      WRITE(*,101)'Enter region to be selected from new frame:'
81      WRITE(*,100)'Scan region (0,0=EXIT)...'
        CALL READ2I('@',NS1,NS2)
        IF((NS1.EQ.0).AND.(NS2.EQ.0)) GOTO 15
        IF((NS1.LT.1).OR.(NS2.GT.NSCAN2).OR.(NS1.GT.NS2))THEN
          WRITE(*,101)'ERROR: invalid entry. Try again.'
          GOTO 81
        END IF
82      WRITE(*,100)'Channel region (0,0=EXIT)'
        CALL READ2I('@',NC1,NC2)
        IF((NC1.EQ.0).AND.(NC2.EQ.0)) GOTO 15
        IF((NC1.LT.1).OR.(NC2.GT.NCHAN2).OR.(NC1.GT.NC2))THEN
          WRITE(*,101)'ERROR: invalid entry. Try again.'
          GOTO 82
        END IF
        WRITE(*,101)'Enter origin in original frame'//
     +   ' where new region will be inserted:'
        WRITE(*,100)'Scan '
        WRITE(CDUMMY,*) NS1
        NS0=READILIM(CDUMMY,1,NSCAN)
        WRITE(*,100)'Channel '
        WRITE(CDUMMY,*) NC1
        NC0=READILIM(CDUMMY,1,NCHAN)
        LOUTLIMITS=.FALSE.               !vigila si insertamos fuera de limites
        DO I=NS1,NS2
          II=NS0+(I-NS1)
          IF((II.GE.1).AND.(II.LE.NSCAN))THEN
            DO J=NC1,NC2
              JJ=NC0+(J-NC1)
              IF((JJ.GE.1).AND.(JJ.LE.NCHAN))THEN
                A(JJ,II)=B(J,I)
              ELSE
                LOUTLIMITS=.TRUE.
              END IF
            END DO
          ELSE
            LOUTLIMITS=.TRUE.
          END IF
        END DO
        IF(LOUTLIMITS)THEN
          WRITE(*,101) 'WARNING: some pixels of the selected region '//
     +     'are out of the original frame.' 
          WRITE(*,101) '-------> selected region has been clipped!'
        END IF
        GOTO 83
C------------------------------------------------------------------------------
C introducimos polinomio en los scans solicitados
85      WRITE(*,*)
        WRITE(*,100)'X-coordinate in [w]avelength or [c]hannel '//
     +   'number (w/c) '
        CXTYPE(1:1)=READC('w','wc')
C
        WRITE(*,100) 'Polynomial degree '
        NTERM=READILIM('@',0,NMAXTERM)
        NTERM=NTERM+1
        DO I=1,NTERM
          WRITE(*,'(A2,I2.2,A2,$)') 'a(',I-1,') '
          COEF(I)=READF('@')
        END DO
C
        IF(CXTYPE.EQ.'w')THEN
          WRITE(*,100)'STWV'
          FSTWV=READF('@')
          WRITE(*,100)'DISP'
          FDISP=READF('@')
          STWV=FSTWV
          DISP=FDISP
        ELSE
          FSTWV=1.
          FDISP=1.
          STWV=0.
          DISP=0.
        END IF
C
        WRITE(*,100) 'Computing polynomial...'
        IF(CXTYPE.EQ.'w')THEN
          DO J=1,NCHAN
            XP=REAL(J-1)*FDISP+FSTWV
            YP=COEF(NTERM)
            DO I=NTERM-1,1,-1
              YP=XP*YP+COEF(I)
            END DO
            S(J)=YP
          END DO
        ELSE
          DO J=1,NCHAN
            XP=REAL(J)
            YP=COEF(NTERM)
            DO I=NTERM-1,1,-1
              YP=XP*YP+COEF(I)
            END DO
            S(J)=YP
          END DO
        END IF
        WRITE(*,101) '  ...OK!'
C
        GOTO 22
C------------------------------------------------------------------------------
C multiplicamos el perfil espacial por una gaussiana
87      WRITE(*,*)
        WRITE(*,100) 'X0 (center of the Gaussian ---spatial dir.---)'
        GX0=READF('@')
        WRITE(*,100) 'Gaussian width in pixels' 
        GWIDTH=READF('@')
        DO I=1,NSCAN
          GFACTOR=EXP(-(REAL(I)-GX0)*(REAL(I)-GX0)/(2.*GWIDTH*GWIDTH))
          GFACTOR=GFACTOR/(SQRT(2.0*PI)*GWIDTH)
          DO J=1,NCHAN
            A(J,I)=A(J,I)*GFACTOR
          END DO
        END DO
C
        GOTO 15
C------------------------------------------------------------------------------
90      CONTINUE
        WRITE(*,101) '* Enter gaussian parameters:'
        WRITE(*,101) 'Flux=Amp*exp[-(x-x0)^2/(2*sigma^2)]'
        WRITE(*,100) 'Amp (counts) '
        G_AMP=READF('1.0')
        WRITE(*,100) 'Sigma (pixels) '
        G_SIGMA=READF('2.0')
        WRITE(*,100) 'x0 (pixel)'
        G_X0=READF('@')
C generamos un espectro con la gaussiana anterior
        FACTOR1=-1./(2.*G_SIGMA*G_SIGMA)
        DO J=1,NCHAN
          S(J)=FINTGAUSS(REAL(J)-0.5,REAL(J)+0.5,20,G_X0,FACTOR1)
          S(J)=G_AMP*S(J)
        END DO
92      WRITE(*,100)'Enter scan region to introduce spectrum (0,0=EXIT)'
        CALL READ2I('0,0',NS1,NS2)
        IF((NS1.EQ.0).AND.(NS2.EQ.0)) GOTO 15
        IF((NS1.LT.1).OR.(NS2.GT.NSCAN).OR.(NS1.GT.NS2))THEN
          WRITE(*,101)'ERROR: invalid entry. Try again.'
          GOTO 92
        END IF
        DO I=NS1,NS2
          DO J=1,NCHAN
            A(J,I)=A(J,I)+S(J)
          END DO
        END DO
        GOTO 92
C------------------------------------------------------------------------------
94      CONTINUE
        WRITE(*,101) '* Please, confirm wavelength calibration:'
        WRITE(CDUMMY,*) STWV0
        WRITE(*,100) 'STWV '
        STWV0=READF(CDUMMY)
        WRITE(CDUMMY,*) DISP0
        WRITE(*,100) 'DISP '
        DISP0=READF(CDUMMY)
        WRITE(*,100) 'Name of file with line list '//
     +   '(2 columns: wavel. & flux)? '
        READ(*,101) FILELINES
        INQUIRE(FILE=FILELINES,EXIST=LOGFILE)
        IF(.NOT.LOGFILE)THEN
          WRITE(*,101) 'ERROR: this file does not exist.'
          WRITE(*,100) 'Press <CR> to continue...'
          READ(*,*)
          GOTO 15
        END IF
        LDO_CENTRAL=STWV0+DISP0*REAL(NCHAN)/2.0
        WRITE(CDUMMY,*) DISP0/(LDO_CENTRAL)*C !aprox. 1 pixel
        WRITE(*,100) 'Sigma (km/s) '
        G_SIGMA_KMS=READF(CDUMMY)
        WRITE(*,100) 'Times sigma to extend the tails of the lines '//
     +   'at each side '
        TG_SIGMA=READF('9.0')
C       !contamos el numero de lineas que hay en el fichero
        NLINEST=0
        OPEN(17,FILE=FILELINES,STATUS='OLD',FORM='FORMATTED')
95      READ(17,*,END=96) G_X0_LDO,G_AMP
        NLINEST=NLINEST+1
        GOTO 95
96      CLOSE(17)
        IF(NLINEST.EQ.0)THEN
          WRITE(*,101) 'ERROR: this file is empty!'
          WRITE(*,100) 'press <CR> to continue...'
          READ(*,*)
          GOTO 15
        ELSE
          WRITE(*,100) 'No. of lines read: '
          WRITE(*,*) NLINEST
        END IF
C
97      WRITE(*,100)'Enter scan region to introduce spectrum (0,0=EXIT)'
        CALL READ2I('0,0',NS1,NS2)
        IF((NS1.EQ.0).AND.(NS2.EQ.0)) GOTO 15
        IF((NS1.LT.1).OR.(NS2.GT.NSCAN).OR.(NS1.GT.NS2))THEN
          WRITE(*,101)'ERROR: invalid entry. Try again.'
          GOTO 97
        END IF
C
        OPEN(17,FILE=FILELINES,STATUS='OLD',FORM='FORMATTED')
        NLINES=0
98      READ(17,*,END=99) G_X0_LDO,G_TOTAL_FLUX
        NLINES=NLINES+1
        CALL SHOWPERC(1,NLINEST,1,NLINES,NEXTINFO)
        !generamos un espectro con cada gaussiana
        G_X0=1.0+(G_X0_LDO-STWV0)/DISP0 !.........centro de la linea en pixels
        G_SIGMA=(G_SIGMA_KMS*G_X0_LDO/C)/DISP0    !sigma de la linea en pixels 
          !(calculada en la l.d.o. central); es una aproximacion porque cambia
             !con la l.d.o. y usamos el mismo valor en todo el recorrido de la
                               !linea, pero deberia ser una buena aproximacion
        !calculamos el valor de la amplitud de la gaussiana para que el area
        !total coincida con el flujo leido en el fichero
        G_AMP=G_TOTAL_FLUX/(SQRT(2.0*PI)*G_SIGMA)
        FACTOR1=-1./(2.*G_SIGMA*G_SIGMA)
        DO J=1,NCHAN
          S(J)=0.0
        END DO
        J1=G_X0-TG_SIGMA*G_SIGMA
        J2=G_X0+TG_SIGMA*G_SIGMA
        DO J=J1,J2
          IF((J.GE.1).AND.(J.LE.NCHAN))THEN
            S(J)=FINTGAUSS(REAL(J)-0.5,REAL(J)+0.5,20,G_X0,FACTOR1)
            S(J)=G_AMP*S(J)
          END IF
        END DO
        DO I=NS1,NS2
          DO J=1,NCHAN
            A(J,I)=A(J,I)+S(J)
          END DO
        END DO
        GOTO 98
99      CLOSE(17)
        GOTO 15
C------------------------------------------------------------------------------
C introducimos constante en las regiones solicitadas
300     WRITE(*,*)
        WRITE(*,101)'Enter region and pixel value:'
301     WRITE(*,100)'Scan region (0,0=EXIT)...'
        CALL READ2I('@',NS1,NS2)
        IF((NS1.EQ.0).AND.(NS2.EQ.0)) GOTO 15
        IF((NS1.LT.1).OR.(NS2.GT.NSCAN).OR.(NS1.GT.NS2))THEN
          WRITE(*,101)'ERROR: invalid entry. Try again.'
          GOTO 301
        END IF
302     WRITE(*,100)'Channel region (0,0=EXIT)'
        CALL READ2I('@',NC1,NC2)
        IF((NC1.EQ.0).AND.(NC2.EQ.0)) GOTO 15
        IF((NC1.LT.1).OR.(NC2.GT.NCHAN).OR.(NC1.GT.NC2))THEN
          WRITE(*,101)'ERROR: invalid entry. Try again.'
          GOTO 302
        END IF
        WRITE(*,100)'Seed for random number generator '//
     +   '(-1: internal clock, > 0: fixed) '
        NSEED=READI('-1')
        WRITE(*,100)'Standard deviation'
        STDERR=READF('@')
        DO I=NS1,NS2
          DO J=NC1,NC2
            R1=RANRED(NSEED)
            R2=RANRED(NSEED)
            A(J,I)=A(J,I)+
     +       SQRT2*STDERR*SQRT(-1.*ALOG(1.-R1))*COS(PI2*R2)
          END DO
        END DO
        GOTO 15
C------------------------------------------------------------------------------
310     WRITE(*,*)
        WRITE(*,100)'File name to be used for STD'
        INFILE=INFILEX(20,'@',NSCAN2,NCHAN2,STWV2,DISP2,1,.FALSE.)
        IF((NSCAN2.NE.NSCAN).OR.(NCHAN2.NE.NCHAN))THEN
          WRITE(*,101)'ERROR: dimensions in last image are different.'
          WRITE(*,100) 'Press <CR> to continue...'
          READ(*,*)
          GOTO 15
        END IF
        IF(STWV2.NE.STWV0)THEN
          WRITE(*,101)'WARNING: STWV in last image is different.'
          WRITE(*,101)'Do you want to continue anyway (y/n) '
          CCONT(1:1)=READC('n','yn')
          IF(CCONT.EQ.'n') GOTO 15
        END IF
        IF(DISP2.NE.DISP0)THEN
          WRITE(*,101)'WARNING: DISP in last image is different.'
          WRITE(*,101)'Do you want to continue anyway (y/n) '
          CCONT(1:1)=READC('n','yn')
          IF(CCONT.EQ.'n') GOTO 15
        END IF
        DO I=1,NSCAN2
          READ(20) (B(J,I),J=1,NCHAN2)
        END DO
        CLOSE(20)
C
        WRITE(*,100)'Seed for random number generator '//
     +   '(-1: internal clock, > 0: fixed) '
        NSEED=READI('-1')
        DO I=1,NSCAN
          DO J=1,NCHAN
            R1=RANRED(NSEED)
            R2=RANRED(NSEED)
            A(J,I)=A(J,I)+
     +       SQRT2*B(J,I)*SQRT(-1.*ALOG(1.-R1))*COS(PI2*R2)
          END DO
        END DO
        GOTO 15
C
C------------------------------------------------------------------------------

C establecemos unos parametros de cabecera por defecto
900     STWV=STWV0
        DISP=DISP0
        AIRMASS=AIRMASS0
        TIMEXPOS=TIMEXPOS0
        OBJECT=OBJECT0
        FITSFILE=FITSFILE0
        COMMENT=COMMENT0
C------------------------------------------------------------------------------
C mostramos valores actuales de la cabecera
        WRITE(*,*)
        WRITE(CDUMMY,'(I10,A1,I10)')NSCAN,',',NCHAN
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,101)'Image size (NSCAN,NCHAN): '//CDUMMY(1:L)
        WRITE(CDUMMY,*)STWV
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,101)'STWV    : '//CDUMMY(1:L)
        WRITE(CDUMMY,*)DISP
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,101)'DISP    : '//CDUMMY(1:L)
        WRITE(CDUMMY,*)AIRMASS
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,101)'Airmass : '//CDUMMY(1:L)
        WRITE(CDUMMY,*)TIMEXPOS
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,101)'Timexpos: '//CDUMMY(1:L)
        L=TRUELEN(OBJECT)
        IF(L.GT.0)THEN
          WRITE(*,101)'Object  : '//OBJECT(1:TRUELEN(OBJECT))
          WRITE(*,110)'No. of characters: ',L
        ELSE
          OBJECT=CHAR(32)
          WRITE(*,101)'Object  : [not found]'
        END IF
        L=TRUELEN(FITSFILE)
        IF(L.GT.0)THEN
          WRITE(*,101)'FITSfile: '//FITSFILE(1:TRUELEN(FITSFILE))
          WRITE(*,110)'No. of characters: ',L
        ELSE
          FITSFILE=CHAR(32)
          WRITE(*,101)'FITSfile: [not found]'
        END IF
        L=TRUELEN(COMMENT)
        IF(L.GT.0)THEN
          WRITE(*,101)'Comment : '//COMMENT(1:TRUELEN(COMMENT))
          WRITE(*,110)'No. of characters: ',L
        ELSE
          COMMENT=CHAR(32)
          WRITE(*,101)'Comment : [not found]'
        END IF
        WRITE(*,*)
C------------------------------------------------------------------------------
C pedimos confirmacion de la cabecera
        WRITE(*,100)'Change head information (y/n) '
        CCHANGE(1:1)=READC('n','yn')
C
        IF(CCHANGE.EQ.'y')THEN
          WRITE(CDUMMY,*)STWV
          WRITE(*,100)'STWV     '
          STWV=READF(CDUMMY)
c
          WRITE(CDUMMY,*)DISP
          WRITE(*,100)'DISP     '
          DISP=READF(CDUMMY)
c
          WRITE(CDUMMY,*)AIRMASS
          WRITE(*,100)'AIRMASS  '
          AIRMASS=READF(CDUMMY)
c
          WRITE(CDUMMY,*)TIMEXPOS
          WRITE(*,100)'TIMEXPOS '
          TIMEXPOS=READF(CDUMMY)
c
          L=TRUELEN(OBJECT)
          IF(L.EQ.0)THEN
            WRITE(*,101)'OBJECT  : [not found] '
          ELSE
            WRITE(*,101)'OBJECT  : '//OBJECT(1:L)
          END IF
          WRITE(*,100)'Change OBJECT value.....(y/n) '
          CCC(1:1)=READC('n','yn')
          IF(CCC.EQ.'y')THEN
            WRITE(*,100)'OBJECT   '
            OBJECT=READC('@','@')
          END IF
c
          L=TRUELEN(FITSFILE)
          IF(L.EQ.0)THEN
            WRITE(*,101)'FITSFILE: [not found] '
          ELSE
            WRITE(*,101)'FITSFILE: '//FITSFILE(1:L)
          END IF
          WRITE(*,100)'Change FITSFILE value...(y/n) '
          CCC(1:1)=READC('n','yn')
          IF(CCC.EQ.'y')THEN
            WRITE(*,100)'FITSFILE '
            FITSFILE=READC('@','@')
          END IF
c
          L=TRUELEN(COMMENT)
          IF(L.EQ.0)THEN
            WRITE(*,100)'COMMENT : [not found] '
          ELSE
            WRITE(*,101)'COMMENT : '
            WRITE(*,101)COMMENT(1:L)
          END IF
          WRITE(*,100)'Change COMMENT value....(y/n) '
          CCC(1:1)=READC('n','yn')
          IF(CCC.EQ.'y')THEN
            WRITE(*,101)'COMMENT  '
            COMMENT=READC('@','@')
          END IF
        END IF
C------------------------------------------------------------------------------
C salvamos finalmente el fichero creado
        WRITE(*,*)
        WRITE(*,100)'Output file name'
        OUTFILE=OUTFILEX(30,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        DO I=1,NSCAN
          WRITE(30) (A(J,I),J=1,NCHAN)
        END DO
        CLOSE(30)
C
        STOP
100     FORMAT(A,$)
101     FORMAT(A)
110     FORMAT(A,I6)
C
        END
