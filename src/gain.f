C------------------------------------------------------------------------------
C Version 07-September-2007                                         file:gain.f
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
C Program: gain
C Classification: arithmetic & manipulations
C Description: Measure the gain using several flatfield images.
C
Comment
C
C Descripcion: el programa utiliza varias imagenes de flatfield (obtenidas con
C el mismo tiempo de exposicion y, por tanto, con la misma se\~{n}al) para 
C determinar la ganancia. Para ello, el programa mide la varianza en muchas
C regiones peque\~{n}as (cuadradas) de las imagenes, en las que se supone que
C la se\~{n}al no varia demasiado. En cada uno de estos cuadrados, el programa
C determina el cociente entre la imagen considerada y la imagen promedio
C (obtenida como suma de todas las imagenes individuales disponibles). Dicho
C cociente permite escalar la imagen suma y sustraerla a la imagen individual.
C De esta forma, en el cuadrado considerado obtenemos "algo" que debe ser solo
C el ruido de la imagen (varianza=1/Gain*Ncuentas+ReadoutNoise), dado que la
C sustraccion anterior debe haber eliminado el efecto de la respuesta pixel a
C pixel. Sin embargo, como el numero de imagenes utilizado no es tremendamente
C grande, la imagen promedio todavia se ve afectada del ruido (el efecto
C disminuye con 1/SQRT[N]). Por eso, lo que representamos graficamente es el
C ruido medido dividido por (1-1/N). Una vez realizado esto para cuadrados de
C diferente tama\~{n}o, el programa permite repetir el proceso utilizando N-1
C imagenes (es decir, dejando fuera una de ellas). El resultado obtenido debe
C coincidir dentro de los errores estadisticos.
C
        PROGRAM GAIN
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER TRUELEN
        INTEGER READILIM
C
        INTEGER NMAXFILES
        PARAMETER (NMAXFILES=80)
        INTEGER NMAXBOXES
        PARAMETER (NMAXBOXES=500)
        INTEGER NMAXSIZE
        PARAMETER (NMAXSIZE=101)
C
        REAL RANRED
        REAL FMEAN0
        REAL FMEDIAN1
C
        INTEGER I,J,K,L,N
        INTEGER NPIX,NF
        INTEGER IBOX(NMAXBOXES),JBOX(NMAXBOXES)
        INTEGER I1,I2,J1,J2
        INTEGER NFILES,NFILESEFF,NOUT,NOUTLIMIT
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        INTEGER NSCAN_,NCHAN_
        INTEGER NC1BIAS,NC2BIAS,NS1BIAS,NS2BIAS
        INTEGER NC1,NC2,NS1,NS2
        INTEGER NSIZE,DSIZE,NSIZE1,NSIZE2
        INTEGER NBOXES,NGAIN
        INTEGER NSEED
        REAL STWV_,DISP_
        REAL A(NCMAX,NSMAX)
        REAL PIXEL(NMAXSIZE*NMAXSIZE)
        REAL SUMA(NCMAX,NSMAX)
        REAL BOXMEDIAN(NMAXBOXES),BOXMEDIAN_
        REAL FACTOR
        REAL XP(NMAXFILES*NMAXBOXES),YP(NMAXFILES*NMAXBOXES)
        REAL XDUM,EXDUM,YDUM,EYDUM
        REAL XXDUM(1),YYDUM(1)
        REAL XMIN,XMAX,YMIN,YMAX,DX,DY
        REAL XMIN_,XMAX_,YMIN_,YMAX_,DX_,DY_
        REAL G(NMAXSIZE/2+1),ERRG(NMAXSIZE/2+1),XG(NMAXSIZE/2+1)
        REAL GAINFIN(NMAXSIZE/2+1,NMAXFILES+1)
        REAL ERRGAINFIN(NMAXSIZE/2+1,NMAXFILES+1)
        REAL MEANW,ERRMEAN,SUMERR,ERR_RMS,ERR_EXP
        DOUBLE PRECISION BIAS(NMAXFILES)
        DOUBLE PRECISION XF(NMAXBOXES),EXF(NMAXBOXES)
        DOUBLE PRECISION YF(NMAXBOXES),EYF(NMAXBOXES)
        DOUBLE PRECISION COEFA,COEFB,VCOEFA,VCOEFB
        CHARACTER*1 CNOUT
        CHARACTER*50 CDUMMY
        CHARACTER*75 LISTFILES,FILEIN(NMAXFILES),FILEDUM
        LOGICAL LCOLOR(MAX_ID_RED)
        LOGICAL LOGFILE
        LOGICAL IFCHAN_BIAS(NCMAX),IFSCAN_BIAS(NSMAX)
C------------------------------------------------------------------------------
        OUTFILEX=OUTFILEX
C
        THISPROGRAM='gain'
        CALL WELCOME('07-September-2007')
C------------------------------------------------------------------------------
        IF(NCMAX.LT.NMAXFILES)THEN
          WRITE(*,101) 'FATAL ERROR: NCMAX.LT.NMAXFILES'
          STOP
        END IF
C------------------------------------------------------------------------------
C leemos lista de ficheros
        WRITE(*,100)'Name of file with list of images to examine'
        LISTFILES=INFILEX(10,'@',0,0,0.,0.,3,.FALSE.)
        N=1
10      READ(10,101,END=12) FILEIN(N)
        N=N+1
        GOTO 10
12      CLOSE(10)
        NFILES=N-1
        WRITE(*,100) '>>> No. of files to be read: '
        WRITE(*,*) NFILES
        IF(NFILES.GT.NMAXFILES)THEN
          WRITE(*,101) 'FATAL ERROR: no. of files > NMAXFILES'
          STOP
        ELSEIF(NFILES.LT.4)THEN
          WRITE(*,101) 'FATAL ERROR: no. of files < 4'
          STOP
        END IF
C------------------------------------------------------------------------------
C chequeamos que los ficheros existen
        DO N=1,NFILES
          INQUIRE(FILE=FILEIN(N),EXIST=LOGFILE)
          IF(.NOT.LOGFILE)THEN
            WRITE(*,101) 'FATAL ERROR: the following file does '//
     >       'not exist: '
            WRITE(*,101) FILEIN(N)
            STOP
          ELSE
            L=TRUELEN(FILEIN(N))
            WRITE(*,100) FILEIN(N)(1:L)
            WRITE(*,101) '...The file exists. OK!'
          END IF
        END DO
        WRITE(*,100) 'Press <CR> to continue...'
        READ(*,*)
C------------------------------------------------------------------------------
C comprobamos que el tama\~{n}o de las imagenes es constante
        WRITE(*,101) 'Checking image size:'
        FILEDUM=
     +   INFILEX(20,FILEIN(1),NSCAN,NCHAN,STWV,DISP,11,.FALSE.)
        CLOSE(20)
        DO N=2,NFILES
          FILEDUM=
     +     INFILEX(20,FILEIN(N),NSCAN_,NCHAN_,STWV_,DISP_,11,.FALSE.)
          CLOSE(20)
          IF(NSCAN_.NE.NSCAN)THEN
            WRITE(*,101) 'FATAL ERROR: NSCAN is not constant.'
            STOP
          END IF
          IF(NCHAN_.NE.NCHAN)THEN
            WRITE(*,101) 'FATAL ERROR: NCHAN is not constant.'
            STOP
          END IF
          IF(STWV_.NE.STWV)THEN
            WRITE(*,101) 'FATAL ERROR: STWV is not constant.'
            STOP
          END IF
          IF(DISP_.NE.DISP)THEN
            WRITE(*,101) 'FATAL ERROR: DISP is not constant.'
            STOP
          END IF
        END DO
C------------------------------------------------------------------------------
C Definimos region para medir bias
        WRITE(*,*)
        WRITE(*,101) '* Define region to measure BIAS'
C
        DO J=1,NCHAN
          IFCHAN_BIAS(J)=.FALSE.
        END DO
        NC1BIAS=1
        NC2BIAS=1
        DO WHILE((NC1BIAS.NE.0).AND.(NC2BIAS.NE.0))
          WRITE(*,100) 'Channel region (0,0=EXIT) '
          CALL READ2I('0,0',NC1BIAS,NC2BIAS)
          IF((NC1BIAS.EQ.0).AND.(NC2BIAS.EQ.0))THEN
          ELSEIF((NC1BIAS.GE.1).AND.(NC2BIAS.LE.NCHAN))THEN
            DO J=NC1BIAS,NC2BIAS
              IFCHAN_BIAS(J)=.TRUE.
            END DO
          ELSE
            WRITE(*,101) 'ERROR: numbers out of range. Try again.'
            WRITE(*,100) 'Press <CR> to continue...'
            READ(*,*)
          END IF
        END DO
C
        DO I=1,NSCAN
          IFSCAN_BIAS(I)=.FALSE.
        END DO
        NS1BIAS=1
        NS2BIAS=1
        DO WHILE((NS1BIAS.NE.0).AND.(NS2BIAS.NE.0))
          WRITE(*,100) 'Scan region (0,0=EXIT) '
          CALL READ2I('0,0',NS1BIAS,NS2BIAS)
          IF((NS1BIAS.EQ.0).AND.(NS2BIAS.EQ.0))THEN
          ELSEIF((NS1BIAS.GE.1).AND.(NS2BIAS.LE.NSCAN))THEN
            DO I=NS1BIAS,NS2BIAS
              IFSCAN_BIAS(I)=.TRUE.
          END DO
            ELSE
            WRITE(*,101) 'ERROR: numbers out of range. Try again.'
            WRITE(*,100) 'Press <CR> to continue...'
            READ(*,*)
          END IF
        END DO
C------------------------------------------------------------------------------
C calculamos el bias de cada imagen
        DO N=1,NFILES
          WRITE(*,'(A,I2.2,A1,I2.2)') 'Measuring BIAS in frame #',
     +     N,'/',NFILES
          FILEDUM=
     +     INFILEX(20,FILEIN(N),NSCAN,NCHAN,STWV,DISP,11,.FALSE.)
          DO I=1,NSCAN
            READ(20) (A(J,I),J=1,NCHAN)
          END DO
          CLOSE(20)
          K=0
          BIAS(N)=0.D0
          DO I=1,NSCAN
            IF(IFSCAN_BIAS(I))THEN
              DO J=1,NCHAN
                IF(IFCHAN_BIAS(J))THEN
                  K=K+1
                  BIAS(N)=BIAS(N)+DBLE(A(J,I))
                END IF
              END DO
            END IF
          END DO
          IF(K.GT.0)THEN
            BIAS(N)=BIAS(N)/DBLE(K)
          ELSE
            BIAS(N)=0.0D0
          END IF
          DO I=1,NSCAN
            DO J=1,NCHAN
              SUMA(J,I)=SUMA(J,I)+A(J,I)-REAL(BIAS(N))
            END DO
          END DO
        END DO
C mostramos el valor de bias medido en cada imagen
        DO N=1,NFILES
          WRITE(*,'(A,I2.2,A,$)') 'File#, bias: ',N,': '
          WRITE(*,*) BIAS(N)
        END DO
        WRITE(*,100) 'Press <CR> to continue...'
        READ(*,*)
C------------------------------------------------------------------------------
C Definimos region util de las imagenes
        WRITE(*,*)
        WRITE(*,101) '* Define useful image region (rectangle)'
C
        WRITE(CDUMMY,'(A,I10)') '1,',NCHAN
        CALL RMBLANK(CDUMMY,CDUMMY,L)
20      WRITE(*,100) 'Channel region '
        CALL READ2I(CDUMMY(1:L),NC1,NC2)
        IF((NC1.LT.1).OR.(NC2.GT.NCHAN))THEN
          WRITE(*,101) 'ERROR: numbers out of range. Try again.'
          WRITE(*,100) 'Press <CR> to continue...'
          READ(*,*)
          GOTO 20
        END IF
C
        WRITE(CDUMMY,'(A,I10)') '1,',NSCAN
        CALL RMBLANK(CDUMMY,CDUMMY,L)
22      WRITE(*,100) 'Scan region '
        CALL READ2I(CDUMMY(1:L),NS1,NS2)
        IF((NS1.LT.1).OR.(NS2.GT.NSCAN))THEN
          WRITE(*,101) 'ERROR: numbers out of range. Try again.'
          WRITE(*,100) 'Press <CR> to continue...'
          READ(*,*)
          GOTO 22
        END IF
C------------------------------------------------------------------------------
C definimos numero de cajas para realizar la estadistica
        WRITE(*,100) 'Number of boxes '
        NBOXES=READILIM('100',10,NMAXBOXES)
C------------------------------------------------------------------------------
C definimos tama\~{n}o de la caja
30      WRITE(*,100) 'Minimum box size (odd) '
        NSIZE1=READILIM('3',3,NMAXSIZE-2)
        IF(MOD(NSIZE1,2).EQ.0)THEN
          NSIZE1=NSIZE1+1
          WRITE(*,100) 'WARNING: NSIZE1 set to: '
          WRITE(*,*) NSIZE1
          WRITE(*,100) 'Press <CR> to continue...'
          READ(*,*)
        END IF
        IF(NSIZE1.GT.(NS2-NS1+1))THEN
          WRITE(*,100) 'Box size is larger than useful scan region. '
          WRITE(*,101) 'Try again.'
          WRITE(*,100) 'Press <CR> to continue...'
          READ(*,*)
          GOTO 30
        END IF
        IF(NSIZE1.GT.(NC2-NC1+1))THEN
          WRITE(*,100) 'Box size is larger than useful channel region. '
          WRITE(*,101) 'Try again.'
          WRITE(*,100) 'Press <CR> to continue...'
          READ(*,*)
          GOTO 30
        END IF
C
32      WRITE(*,100) 'Maximum box size (odd) '
        WRITE(CDUMMY,*) NMAXSIZE
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        NSIZE2=READILIM(CDUMMY(1:L),NSIZE1+2,NMAXSIZE)
        IF(MOD(NSIZE2,2).EQ.0)THEN
          NSIZE2=NSIZE2+1
          WRITE(*,100) 'WARNING: NSIZE2 set to: '
          WRITE(*,*) NSIZE2
          WRITE(*,100) 'Press <CR> to continue...'
          READ(*,*)
        END IF
        IF(NSIZE2.GT.(NS2-NS1+1))THEN
          WRITE(*,100) 'Box size is larger than useful scan region. '
          WRITE(*,101) 'Try again.'
          WRITE(*,100) 'Press <CR> to continue...'
          READ(*,*)
          GOTO 32
        END IF
        IF(NSIZE2.GT.(NC2-NC1+1))THEN
          WRITE(*,100) 'Box size is larger than useful channel region. '
          WRITE(*,101) 'Try again.'
          WRITE(*,100) 'Press <CR> to continue...'
          READ(*,*)
          GOTO 32
        END IF
C
        WRITE(*,100) 'Step in box size (even) '
        DSIZE=READILIM('2',2,98)
        IF(MOD(DSIZE,2).NE.0)THEN
          DSIZE=DSIZE+1
          WRITE(*,100) 'WARNING: DSIZE set to: '
          WRITE(*,*) DSIZE
          WRITE(*,100) 'Press <CR> to continue...'
          READ(*,*)
        END IF
        XMIN_=REAL(NSIZE1)
        XMAX_=REAL(NSIZE2)
        DX_=XMAX_-XMIN_
        XMIN_=XMIN_-DX_/20.
        XMAX_=XMAX_+DX_/20.
C------------------------------------------------------------------------------
C abrimos la salida grafica
        CALL PIDEGTER(NTERM,IDN,LCOLOR)
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGSCH(1.0)
          CALL PGSLW(2)
          CALL PGSCF(2)
        END DO
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
        WRITE(*,100) 'Are you measuring the gain with N-1 images (y/n) '
        CNOUT(1:1)=READC('y','yn')
        IF(CNOUT.EQ.'n')THEN
          NOUTLIMIT=0
        ELSE
          WRITE(*,100) 'No. of last image to be removed (0=NONE) '
          NOUTLIMIT=READILIM('0',1,NFILES)
        END IF
        NSEED=-1
C BUCLE: en numero de imagen
        NOUT=0
C------------------------------------------------------------------------------
C BUCLE: en tama\~{n}o de la caja
40      NGAIN=0
        NSIZE=NSIZE1
C------------------------------------------------------------------------------
C calculamos la imagen suma
        IF(NOUT.EQ.0)THEN
          NFILESEFF=NFILES
        ELSE
          NFILESEFF=NFILES-1
        END IF
C
        DO I=1,NSCAN
          DO J=1,NCHAN
            SUMA(J,I)=0.0
          END DO
        END DO
        DO N=1,NFILES
          IF(N.NE.NOUT)THEN
            WRITE(*,'(A,I2.2,A1,I2.2,A,I2.2,A)') 
     +       'Computing added image with frame #',N,'/',NFILES,
     +       ' (file #',NOUT,' removed)'
            FILEDUM=
     +       INFILEX(20,FILEIN(N),NSCAN,NCHAN,STWV,DISP,11,.FALSE.)
            DO I=1,NSCAN
              READ(20) (A(J,I),J=1,NCHAN)
            END DO
            CLOSE(20)
            DO I=1,NSCAN
              DO J=1,NCHAN
                SUMA(J,I)=SUMA(J,I)+A(J,I)-REAL(BIAS(N))
              END DO
            END DO
          END IF
        END DO
        DO I=1,NSCAN
          DO J=1,NCHAN
            SUMA(J,I)=SUMA(J,I)/REAL(NFILESEFF)
          END DO
        END DO
C------------------------------------------------------------------------------
C calculamos los centros de las cajas
50      DO K=1,NBOXES
          IBOX(K)=NS1+INT(RANRED(NSEED)*REAL(NS2-NS1)+1)
          IF(IBOX(K)-NS1.LT.NSIZE/2) IBOX(K)=NS1+NSIZE/2
          IF(NS2-IBOX(K).LT.NSIZE/2) IBOX(K)=NS2-NSIZE/2
          JBOX(K)=NC1+INT(RANRED(NSEED)*REAL(NC2-NC1)+1)
          IF(JBOX(K)-NC1.LT.NSIZE/2) JBOX(K)=NC1+NSIZE/2
          IF(NC2-JBOX(K).LT.NSIZE/2) JBOX(K)=NC2-NSIZE/2
          I1=IBOX(K)-NSIZE/2
          I2=IBOX(K)+NSIZE/2
          J1=JBOX(K)-NSIZE/2
          J2=JBOX(K)+NSIZE/2
          NPIX=0
          DO I=I1,I2
            DO J=J1,J2
              NPIX=NPIX+1
            END DO
          END DO
          IF(NPIX.NE.NSIZE*NSIZE)THEN
            WRITE(*,101) 'ERROR: unexpected no. of pixels.'
            WRITE(*,100) 'Box center (I,J): '
            WRITE(*,*) IBOX(K),JBOX(K)
            WRITE(*,100) 'Box limits (I1,I2,J1,J2): '
            WRITE(*,*) I1,I2,J1,J2
            WRITE(*,100) 'Press <CR> to continue...'
            READ(*,*)
          END IF
        END DO
C------------------------------------------------------------------------------
C en la imagen suma, calculamos el numero promedio de cuentas en cada caja
        DO K=1,NBOXES
          I1=IBOX(K)-NSIZE/2
          I2=IBOX(K)+NSIZE/2
          J1=JBOX(K)-NSIZE/2
          J2=JBOX(K)+NSIZE/2
          NPIX=0
          DO I=I1,I2
            DO J=J1,J2
              NPIX=NPIX+1
              PIXEL(NPIX)=SUMA(J,I)
            END DO
          END DO
          BOXMEDIAN(K)=FMEDIAN1(NPIX,PIXEL)
        END DO
C------------------------------------------------------------------------------
C abrimos de nuevo cada imagen individual, calculamos el promedio en cada
C caja, escalamos la imagen suma en dicha caja al valor promedio anterior,
C sustraemos la imagen suma escalada y medimos varianza y valor promedio
        NF=0
        DO N=1,NFILES
          IF(N.NE.NOUT)THEN
            WRITE(*,'(A,I2.2,A1,I2.2,A,I2.2,A)') 
     +       'Measuring in frame #',N,'/',NFILES,' (file #',NOUT,
     +       ' removed)'
C..............................................................................
            FILEDUM=
     +       INFILEX(20,FILEIN(N),NSCAN,NCHAN,STWV,DISP,11,.FALSE.)
            DO I=1,NSCAN
              READ(20) (A(J,I),J=1,NCHAN)
            END DO
            CLOSE(20)
            DO I=1,NSCAN
              DO J=1,NCHAN
                A(J,I)=A(J,I)-REAL(BIAS(N))
              END DO
            END DO
C..............................................................................
            DO K=1,NBOXES
              I1=IBOX(K)-NSIZE/2
              I2=IBOX(K)+NSIZE/2
              J1=JBOX(K)-NSIZE/2
              J2=JBOX(K)+NSIZE/2
              NPIX=0
              DO I=I1,I2
                DO J=J1,J2
                  NPIX=NPIX+1
                  PIXEL(NPIX)=A(J,I)
                END DO
              END DO
              BOXMEDIAN_=FMEDIAN1(NPIX,PIXEL)
              FACTOR=BOXMEDIAN(K)/BOXMEDIAN_
              NPIX=0
              DO I=I1,I2
                DO J=J1,J2
                  NPIX=NPIX+1
                  PIXEL(NPIX)=A(J,I)-SUMA(J,I)*FACTOR
                END DO
              END DO
              XDUM=FMEAN0(NPIX,PIXEL,YDUM)
              NF=NF+1
              XP(NF)=BOXMEDIAN_
C NOTA: dividimos la varianza por (1-1/N) para tener en cuenta el hecho de
C que la imagen promedio se ve afectada por el propio ruido que pretendemos
C medir.
              YP(NF)=YDUM*YDUM/(1.-1./REAL(NFILESEFF))
            END DO
C..............................................................................
          END IF
        END DO
C------------------------------------------------------------------------------
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGPAGE
          CALL PGIDEN_RED
          CALL PGSVP(0.07,0.47,0.6,0.9)
          CALL PGSWIN(-0.6,REAL(NCHAN)+0.6,
     +                -0.6,REAL(NSCAN)+0.6)
          CALL PGBOX('BCNTSI',0.0,0,'BCNTSI',0.0,0)
          CALL PGLABEL('channel','scan',' ')
          CALL PGSCI(5)
          CALL PGSFS(2)
          CALL PGRECT(REAL(NC1),REAL(NC2),REAL(NS1),REAL(NS2))
          CALL PGSCI(2)
          DO K=1,NBOXES
            I1=IBOX(K)-NSIZE/2
            I2=IBOX(K)+NSIZE/2
            J1=JBOX(K)-NSIZE/2
            J2=JBOX(K)+NSIZE/2
            CALL PGRECT(REAL(J1),REAL(J2),REAL(I1),REAL(I2))
          END DO
          CALL PGSFS(1)
          CALL PGSCI(1)
        END DO
C------------------------------------------------------------------------------
        CALL FINDMM(NF,XP,XMIN,XMAX)
        CALL FINDMM(NF,YP,YMIN,YMAX)
        DX=XMAX-XMIN
        DY=YMAX-YMIN
        XMIN=0.
        XMAX=XMAX+DX/20.
        YMIN=0.
        YMAX=YMAX+DY/20.
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGSVP(0.58,0.98,0.6,0.9)
          CALL PGSWIN(XMIN,XMAX,YMIN,YMAX)
          CALL PGBOX('BCNTSI',0.0,0,'BCNTSI',0.0,0)
          CALL PGIDEN_RED
          CALL PGPOINT(NF,XP,YP,1)
          CALL PGLABEL('number of counts','modified variance',' ')
          WRITE(CDUMMY,*) NFILESEFF
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          CALL PGMTEXT('T',-2.0,0.05,0.0,'NFILESEFF: '//CDUMMY(1:L))
          WRITE(CDUMMY,*) NSIZE
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          CALL PGMTEXT('T',-4.0,0.05,0.0,'NSIZE: '//CDUMMY(1:L))
          WRITE(CDUMMY,*) NOUT
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          CALL PGMTEXT('T',-6.0,0.05,0.0,'NOUT: '//CDUMMY(1:L))
        END DO
C
        DO K=1,NBOXES
          XDUM=0.
          YDUM=0.
          DO N=1,NFILESEFF
            NF=K+(N-1)*NBOXES
            XDUM=XDUM+XP(NF)
            YDUM=YDUM+YP(NF)
          END DO
          XDUM=XDUM/REAL(NFILESEFF)
          YDUM=YDUM/REAL(NFILESEFF)
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGSCI(3)
            XXDUM(1)=XDUM
            YYDUM(1)=YDUM
            CALL PGPOINT(1,XXDUM,YYDUM,17)
            CALL PGSCI(1)
          END DO
          XF(K)=DBLE(XDUM)
          YF(K)=DBLE(YDUM)
          EXDUM=0.
          EYDUM=0.
          DO N=1,NFILESEFF
            NF=K+(N-1)*NBOXES
            EXDUM=EXDUM+(XP(NF)-XDUM)*(XP(NF)-XDUM)
            EYDUM=EYDUM+(YP(NF)-YDUM)*(YP(NF)-YDUM)
          END DO
          EXDUM=SQRT(EXDUM/REAL(NFILESEFF-1))
          EYDUM=SQRT(EYDUM/REAL(NFILESEFF-1))
          EXF(K)=DBLE(EXDUM)
          EYF(K)=DBLE(EYDUM)
        END DO
C
        CALL LINREGEY(NBOXES,XF,YF,EYF,COEFA,COEFB,VCOEFA,VCOEFB)
C
        WRITE(*,101) '* Line fit y=a+bx: OLS(X|Y) weighting with EY'
        WRITE(*,100) 'a, sigma(a): '
        WRITE(*,*) COEFA,DSQRT(VCOEFA)
        WRITE(*,100) 'b, sigma(b): '
        WRITE(*,*) COEFB,DSQRT(VCOEFB)
C
        NGAIN=NGAIN+1
        XG(NGAIN)=REAL(NSIZE)
        G(NGAIN)=1.0D0/COEFB
        ERRG(NGAIN)=DSQRT(VCOEFB)/(COEFB*COEFB)
        WRITE(*,100) 'gain (1/b), sigma(gain): '
        WRITE(*,*) G(NGAIN),ERRG(NGAIN)
C
        GAINFIN(NGAIN,NOUT+1)=G(NGAIN)
        ERRGAINFIN(NGAIN,NOUT+1)=ERRG(NGAIN)
C
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGSCI(2)
          CALL PGMOVE(XMIN,REAL(COEFA)+REAL(COEFB)*XMIN)
          CALL PGDRAW(XMAX,REAL(COEFA)+REAL(COEFB)*XMAX)
          CALL PGSCI(1)
        END DO
C
        YMIN_=G(1)-ERRG(1)
        YMAX_=G(1)+ERRG(1)
        IF(NGAIN.GT.1)THEN
          DO I=2,NGAIN
            IF(G(I)-ERRG(I).LT.YMIN_) YMIN_=G(I)-ERRG(I)
            IF(G(I)+ERRG(I).GT.YMAX_) YMAX_=G(I)+ERRG(I)
          END DO
        END IF
        DY_=YMAX_-YMIN_
        YMIN_=YMIN_-DY_/20.
        YMAX_=YMAX_+DY_/20.
C
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGSVP(0.07,0.75,0.15,0.45)
          CALL PGSWIN(XMIN_,XMAX_,YMIN_,YMAX_)
          CALL PGBOX('BCNTS',0.0,0,'BCNTS',0.0,0)
          CALL PGLABEL('box size','gain',' ')
          CALL PGSCI(5)
          DO I=1,NGAIN
            CALL PGPOINT(1,XG(I),G(I),17)
            CALL PGERRY(1,XG(I),G(I)-ERRG(I),G(I)+ERRG(I),1.0)
          END DO
          CALL PGSCI(1)
        END DO
C
        MEANW=0.
        SUMERR=0.
        DO I=1,NGAIN
          SUMERR=SUMERR+1./(ERRG(I)*ERRG(I))
          MEANW=MEANW+G(I)/(ERRG(I)*ERRG(I))
        END DO
        MEANW=MEANW/SUMERR
        ERRMEAN=1./SQRT(SUMERR)
        ERR_EXP=SQRT(REAL(NGAIN))*ERRMEAN
        ERR_RMS=0.
        IF(NGAIN.GT.1)THEN
          DO I=1,NGAIN
            ERR_RMS=ERR_RMS+(G(I)-MEANW)*(G(I)-MEANW)/(ERRG(I)*ERRG(I))
          END DO
          ERR_RMS=REAL(NGAIN)/REAL(NGAIN-1)*ERR_RMS/SUMERR
          ERR_RMS=SQRT(ERR_RMS)
        END IF
C
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
c
          CALL PGMTEXT('T',-2.0,1.10,1.,'<g>')
          CALL PGMTEXT('T',-2.0,1.13,.5,'=')
          WRITE(CDUMMY,*) MEANW
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          CALL PGMTEXT('T',-2.0,1.15,0.,CDUMMY(1:L))
c
          CALL PGMTEXT('T',-4.0,1.10,1.,'\\gs\\d<g>\\u')
          CALL PGMTEXT('T',-4.0,1.13,.5,'=')
          WRITE(CDUMMY,*) ERRMEAN
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          CALL PGMTEXT('T',-4.0,1.15,0.,CDUMMY(1:L))
c
          CALL PGMTEXT('T',-6.0,1.10,1.,'\\gs\\drms\\u')
          CALL PGMTEXT('T',-6.0,1.13,.5,'=')
          WRITE(CDUMMY,*) ERR_RMS
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          CALL PGMTEXT('T',-6.0,1.15,0.,CDUMMY(1:L))
c
          CALL PGMTEXT('T',-8.0,1.10,1.,'\\gs\\dexp\\u')
          CALL PGMTEXT('T',-8.0,1.13,.5,'=')
          WRITE(CDUMMY,*) ERR_EXP
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          CALL PGMTEXT('T',-8.0,1.15,0.,CDUMMY(1:L))
c
          CALL PGMTEXT('T',-10.,1.10,1.,'N\\ddata\\u')
          CALL PGMTEXT('T',-10.,1.13,.5,'=')
          WRITE(CDUMMY,*) NGAIN
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          CALL PGMTEXT('T',-10.,1.15,0.,CDUMMY(1:L))
c          
          CALL PGSCI(7)
          CALL PGSLS(4)
          CALL PGMOVE(XMIN_,MEANW)
          CALL PGDRAW(XMAX_,MEANW)
          CALL PGSLS(1)
          CALL PGSCI(1)
        END DO
C
        NSIZE=NSIZE+DSIZE
        IF(NSIZE.GT.NSIZE2)THEN
          NOUT=NOUT+1
          IF(NOUT.GT.NOUTLIMIT)THEN
            GOTO 90
          ELSE
            GOTO 40
          END IF
        ELSE
          GOTO 50
        END IF
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
90      IF(NOUTLIMIT.EQ.0) GOTO 95
C Dibujamos todas las estimaciones de la ganancia
        YMIN_=GAINFIN(1,1)-ERRGAINFIN(1,1)
        YMAX_=GAINFIN(1,1)+ERRGAINFIN(1,1)
        DO N=0,NOUTLIMIT
          DO I=1,NGAIN
            IF(GAINFIN(I,N+1)-ERRGAINFIN(I,N+1).LT.YMIN_) 
     +       YMIN_=GAINFIN(I,N+1)-ERRGAINFIN(I,N+1)
            IF(GAINFIN(I,N+1)+ERRGAINFIN(I,N+1).GT.YMAX_) 
     +       YMAX_=GAINFIN(I,N+1)+ERRGAINFIN(I,N+1)
          END DO
        END DO
        DY_=YMAX_-YMIN_
        YMIN_=YMIN_-DY_/20.
        YMAX_=YMAX_+DY_/20.
C
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGPAGE
          CALL PGIDEN_RED
          CALL PGSCH(1.5)
          CALL PGENV(XMIN_,XMAX_,YMIN_,YMAX_,0,0)
          CALL PGLABEL('box size','gain',' ')
          DO N=0,NOUTLIMIT
            IF(N.EQ.0)THEN
              CALL PGSCI(2)
            ELSE
              CALL PGSCI(3)
            END IF
            DO I=1,NGAIN
              IF(N.EQ.0)THEN
                CALL PGPOINT(1,XG(I),GAINFIN(I,N+1),21)
              ELSE
                CALL PGPOINT(1,XG(I),GAINFIN(I,N+1),17)
              END IF
              CALL PGERRY(1,XG(I),
     +         GAINFIN(I,N+1)-ERRGAINFIN(I,N+1),
     +         GAINFIN(I,N+1)+ERRGAINFIN(I,N+1),1.0)
            END DO
          END DO
          CALL PGSCI(1)
        END DO
C
        MEANW=0.
        SUMERR=0.
        DO I=1,NGAIN
          SUMERR=SUMERR+1./(ERRGAINFIN(I,1)*ERRGAINFIN(I,1))
          MEANW=MEANW+GAINFIN(I,1)/(ERRGAINFIN(I,1)*ERRGAINFIN(I,1))
        END DO
        MEANW=MEANW/SUMERR
        ERRMEAN=1./SQRT(SUMERR)
        ERR_EXP=SQRT(REAL(NGAIN))*ERRMEAN
        ERR_RMS=0.
        IF(NGAIN.GT.1)THEN
          DO I=1,NGAIN
            ERR_RMS=ERR_RMS+
     +       (GAINFIN(I,1)-MEANW)*(GAINFIN(I,1)-MEANW)/
     +       (ERRGAINFIN(I,1)*ERRGAINFIN(I,1))
          END DO
          ERR_RMS=REAL(NGAIN)/REAL(NGAIN-1)*ERR_RMS/SUMERR
          ERR_RMS=SQRT(ERR_RMS)
        END IF
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGSCI(2)
          CALL PGSLS(4)
          CALL PGMOVE(XMIN_,MEANW)
          CALL PGDRAW(XMAX_,MEANW)
          CALL PGSLS(1)
          CALL PGSCI(1)
        END DO
        WRITE(*,*)
        WRITE(*,'(A,I2.2,A)') '* Results using ',NFILES,' frames:'
        WRITE(*,100) 'mean..........: '
        WRITE(*,*) MEANW
        WRITE(*,100) 'mean error....: '
        WRITE(*,*) ERRMEAN
        WRITE(*,100) 'expected error: '
        WRITE(*,*) ERR_EXP
        WRITE(*,100) 'r.m.s. error..: '
        WRITE(*,*) ERR_RMS
        WRITE(*,100) 'Ndata.........: '
        WRITE(*,*) NGAIN
C
        MEANW=0.
        SUMERR=0.
        DO N=1,NOUTLIMIT
          DO I=1,NGAIN
            SUMERR=SUMERR+1./(ERRGAINFIN(I,N+1)*ERRGAINFIN(I,N+1))
            MEANW=MEANW+GAINFIN(I,N+1)/
     +       (ERRGAINFIN(I,N+1)*ERRGAINFIN(I,N+1))
          END DO
        END DO
        MEANW=MEANW/SUMERR
        ERRMEAN=1./SQRT(SUMERR)
        ERR_EXP=SQRT(REAL(NGAIN*NOUTLIMIT))*ERRMEAN
        ERR_RMS=0.
        IF(NGAIN*NOUTLIMIT.GT.1)THEN
          DO N=1,NOUTLIMIT
            DO I=1,NGAIN
              ERR_RMS=ERR_RMS+
     +         (GAINFIN(I,N+1)-MEANW)*(GAINFIN(I,N+1)-MEANW)/
     +         (ERRGAINFIN(I,N+1)*ERRGAINFIN(I,N+1))
            END DO
          END DO
          ERR_RMS=REAL(NGAIN*NOUTLIMIT)/REAL(NGAIN*NOUTLIMIT-1)*
     +     ERR_RMS/SUMERR
          ERR_RMS=SQRT(ERR_RMS)
        END IF
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGSCI(3)
          CALL PGSLS(4)
          CALL PGMOVE(XMIN_,MEANW)
          CALL PGDRAW(XMAX_,MEANW)
          CALL PGSLS(1)
          CALL PGSCI(1)
        END DO
        WRITE(*,*)
        WRITE(*,'(A,I2.2,A,I2.2,A)') 
     +   '* Results using ',NFILES-1,' frames (',NOUTLIMIT,' times): '
        WRITE(*,100) 'mean..........: '
        WRITE(*,*) MEANW
        WRITE(*,100) 'mean error....: '
        WRITE(*,*) ERRMEAN
        WRITE(*,100) 'expected error: '
        WRITE(*,*) ERR_EXP
        WRITE(*,100) 'r.m.s. error..: '
        WRITE(*,*) ERR_RMS
        WRITE(*,100) 'Ndata.........: '
        WRITE(*,*) NGAIN*NOUTLIMIT
C------------------------------------------------------------------------------
C cerramos la salida grafica
95      CALL PGEND
C------------------------------------------------------------------------------
C fin del programa
        STOP
100     FORMAT(A,$)
101     FORMAT(A)
        END
C
C******************************************************************************
C Version 20-September-1999 
C------------------------------------------------------------------------------
C Copyright: N. Cardiel, Departamento de Astrofisica
C            Facultad de Ciencias Fisicas
C            Universidad Complutense de Madrid, 28040-Madrid, Spain
C            E-mail: ncl@astrax.fis.ucm.es
C------------------------------------------------------------------------------
C This program is free software; you can redistribute it and/or modify it
C under the terms of the GNU General Public License as published by the Free
C Software Foundation; either version 2 of the License, or (at your option) any
C later version. See the file gnu-public-license.txt for details.
C------------------------------------------------------------------------------
C
C     SUBROUTINE LINREGEY(N,X,Y,EY,A,B,VARA,VARB)
C
C
C Linear regresion fit weighting with the errors EY. The formulae are identical
C to those shown in "Data reduction and error analysis for Physical Sciences",
C Philip R. Bevington, chapter 6.
C 
C
C******************************************************************************
C Input:    N        : number of data pairs
C           X(N),Y(N): data pairs
C           EY(N)    : errors
C Output:   A        : intercept coefficient
C           B        : slope
C           VARA     : variance of the intercept coefficient
C           VARB     : variance of the slope
C******************************************************************************
C
        SUBROUTINE LINREGEY(N,X,Y,EY,A,B,VARA,VARB)
        IMPLICIT NONE
C routine parameters
        INTEGER N
        DOUBLE PRECISION X(N),Y(N),EY(N)
        DOUBLE PRECISION A,B
        DOUBLE PRECISION VARA,VARB
C local variables
        INTEGER I
        DOUBLE PRECISION DELTA,SERR
        DOUBLE PRECISION SX,SY,SXX,SXY
C------------------------------------------------------------------------------
        SERR=0.D0
        DO I=1,N
          SERR=SERR+1.D0/(EY(I)*EY(I))
        END DO
C
        SX=0.D0
        DO I=1,N
          SX=SX+X(I)/(EY(I)*EY(I))
        END DO
C
        SY=0.D0
        DO I=1,N
          SY=SY+Y(I)/(EY(I)*EY(I))
        END DO
C
        SXX=0.D0
        DO I=1,N
          SXX=SXX+X(I)*X(I)/(EY(I)*EY(I))
        END DO
C
        SXY=0.D0
        DO I=1,N
          SXY=SXY+X(I)*Y(I)/(EY(I)*EY(I))
        END DO
C
        DELTA=SERR*SXX-SX*SX
C------------------------------------------------------------------------------
        B=(SERR*SXY-SX*SY)/DELTA
        VARB=SERR/DELTA
        A=SY/SERR-B*SX/SERR
        VARA=SXX/DELTA
C------------------------------------------------------------------------------
        END
