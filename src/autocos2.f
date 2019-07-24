C------------------------------------------------------------------------------
C Version 14-October-1997                                      file:autocos2.f
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
C Program: autocos2
C Classification: cosmic rays
C Description: Automatic removal of cosmic rays in many similar images 
C simultaneously (maximum number of images=50, maximum number of cosmic
C rays=1000).
C
Comment
C
C
        PROGRAM AUTOCOS2
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER TRUELEN
        INTEGER READILIM
        REAL READF
C
        INTEGER NBINMAX
        PARAMETER (NBINMAX=1000)
C
        INTEGER NF,I,II,III,K,J,NPIXCR,NCOSF
        INTEGER NCOS(50),NC,NCR
        INTEGER XPOS(50,1000),YPOS(50,1000)
        INTEGER NBIN
        INTEGER NY(NBINMAX)
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        INTEGER NS1,NS2,NC1,NC2
        REAL TSIGMA,FILNUM,MEAN,SIG,MAX,TSIGMA2
        REAL A(NCMAX,NSMAX),B(NCMAX,NSMAX),SB(NCMAX,NSMAX)
        REAL TIMES(50,1000)
        REAL TSIGMA0,TSIGMAMIN,TSIGMAMAX,BMIN,BMAX,DB
        REAL YMIN,YMAX,DY
        REAL X(NBINMAX),Y(NBINMAX)
        REAL TSLIMIT
        CHARACTER*1 CERR,LOUT,CALL,CREPLOT,CLIST
        CHARACTER*50 CDUMMY
        CHARACTER*75 FILEN(50),LOGFILE,LISTFILES
        CHARACTER*75 FILE1,FILE2,FILE3,FILENAME
        LOGICAL LCOLOR(MAX_ID_RED)
        LOGICAL CRFOUND
        LOGICAL LEXIST
        LOGICAL LANY,IFSCAN(NSMAX),IFCHAN(NCMAX)
C
        COMMON/BLKAUTOCOS1/NSCAN,NCHAN
        COMMON/BLKAUTOCOS2/A,B,SB
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
        THISPROGRAM='autocos2'
        CALL WELCOME('14-October-1997')
C------------------------------------------------------------------------------
        NBIN=NBINMAX
C
        WRITE(*,100)'Work with error files (y/n) '
        CERR(1:1)=READC('n','yn')
C
        CALL PIDEGTER(NTERM,IDN,LCOLOR)
C
        WRITE(*,100)'Name of file with list of images to examine'
        LISTFILES=INFILEX(9,'@',0,0,0.,0.,3,.FALSE.)
C------------------------------------------------------------------------------
C Abrimos todos los ficheros y comprobamos que los ficheros de destino no
C existen todavia (en caso contrario el programa se detiene con un mensaje
C de error).
        WRITE(*,101)'Opening files...'
C
        DO I=1,50
          NCOS(I)=0
          READ(9,101,END=999) FILEN(I)
          FILE1=FILEN(I)
          FILE2=FILE1(1:TRUELEN(FILE1))//'c'
          IF(CERR.EQ.'y') CALL GUESSEF(FILE2,FILE3)
          INQUIRE(FILE=FILE1,EXIST=LEXIST)
          IF(.NOT.LEXIST)THEN
            WRITE(*,101)'ERROR: file '//FILE1(1:TRUELEN(FILE1))//
     +       'does not exist.'
            IF(I.GT.1)THEN
              DO II=1,I-1
                CLOSE(10+I)
              END DO
            END IF
            STOP
          END IF
          INQUIRE(FILE=FILE2,EXIST=LEXIST)
          IF(LEXIST)THEN
            WRITE(*,101)'ERROR: file '//FILE2(1:TRUELEN(FILE2))//
     +       'already exist.'
            IF(I.GT.1)THEN
              DO II=1,I-1
                CLOSE(10+I)
              END DO
            END IF
            STOP
          END IF
          IF(CERR.EQ.'y')THEN
            INQUIRE(FILE=FILE3,EXIST=LEXIST)
            IF(LEXIST)THEN
              WRITE(*,101)'ERROR: file '//FILE3(1:TRUELEN(FILE3))//
     +         'already exist.'
              IF(I.GT.1)THEN
                DO II=1,I-1
                  CLOSE(10+I)
                  END DO
              END IF
              STOP
            END IF
          END IF
          WRITE(*,100)'* Opening: '
          WRITE(*,101)FILEN(I)(1:TRUELEN(FILEN(I)))
          FILENAME=
     +     INFILEX(10+I,FILE1,NSCAN,NCHAN,STWV,DISP,11,.FALSE.)
        END DO
C Si el numero de ficheros es mayor que 50, el programa se detiene
        WRITE(*,101)'ERROR: MAX NUMBER OF FILES IS 50'
        DO I=1,50
          CLOSE(10+I)
        END DO
        STOP
C------------------------------------------------------------------------------
C exigimos que al menos tengamos 3 ficheros para trabajar
999     NF=I-1
        FILNUM=FLOAT(NF)
        WRITE(*,100)'files: '
        WRITE(*,*)NF
        IF(NF.LT.3)THEN
          DO I=1,NF
            CLOSE(10+I)
          END DO
          WRITE(*,101)'FATAL ERROR: number of files must be .GT.2'
          STOP
        END IF
        CLOSE(9)
C------------------------------------------------------------------------------
C Determinamos la region que va a ser limpiada
        WRITE(*,100)'Are you cleaning all the pixels (y/n) '
        CALL(1:1)=READC('y','yn')
        IF(CALL.EQ.'y')THEN
          DO I=1,NSCAN
            IFSCAN(I)=.TRUE.
          END DO
          DO J=1,NCHAN
            IFCHAN(J)=.TRUE.
          END DO
        ELSE
          WRITE(*,101)'Define region to be cleaned:'
C
          DO I=1,NSCAN
            IFSCAN(I)=.FALSE.
          END DO
          LANY=.FALSE.
          DO WHILE(.NOT.LANY)
            NS1=1
            NS2=1
            DO WHILE((NS1.NE.0).AND.(NS2.NE.0))
              WRITE(*,100)'Scan region (0,0=EXIT) '
              CALL READ2I('0,0',NS1,NS2)
              IF((NS1.EQ.0).AND.(NS2.EQ.0))THEN
                WRITE(*,101)'Thanks!'
              ELSEIF((NS1.LT.1).OR.(NS2.GT.NSCAN).OR.(NS1.GT.NS2))THEN
                WRITE(*,101)'ERROR: numbers out of range. Try again.'
              ELSE
                DO I=NS1,NS2
                  IFSCAN(I)=.TRUE.
                END DO
              END IF
            END DO
            DO I=1,NSCAN
              IF(IFSCAN(I)) LANY=.TRUE.
            END DO
            IF(.NOT.LANY)THEN
              WRITE(*,101)'ERROR: total number of scans = 0!'
              WRITE(*,101)'Try again.'
            END IF
          END DO
C
          DO J=1,NCHAN
            IFCHAN(J)=.FALSE.
          END DO
          LANY=.FALSE.
          DO WHILE(.NOT.LANY)
            NC1=1
            NC2=1
            DO WHILE((NC1.NE.0).AND.(NC2.NE.0))
              WRITE(*,100)'Channel region (0,0=EXIT) '
              CALL READ2I('0,0',NC1,NC2)
              IF((NC1.EQ.0).AND.(NC2.EQ.0))THEN
                WRITE(*,101)'Thanks!'
              ELSEIF((NC1.LT.1).OR.(NC2.GT.NCHAN).OR.(NC1.GT.NC2))THEN
                WRITE(*,101)'ERROR: numbers out of range. Try again.'
              ELSE
                DO J=NC1,NC2
                  IFCHAN(J)=.TRUE.
                END DO
              END IF
            END DO
            DO J=1,NCHAN
              IF(IFCHAN(J)) LANY=.TRUE.
            END DO
            IF(.NOT.LANY)THEN
              WRITE(*,101)'ERROR: total number of channels = 0!'
              WRITE(*,101)'Try again.'
            END IF
          END DO
C
        END IF
C------------------------------------------------------------------------------
C Calculamos MEAN y SIGMA en cada pixel
C
        TSIGMAMIN=+1.E10
        TSIGMAMAX=-TSIGMAMIN
C
        DO I=1,NSCAN
          WRITE(*,'(A,I4,$)')'\b\b\b\b',I
          DO II=1,NF
            READ(II+10) (A(J,II),J=1,NCHAN)
          END DO
          IF(IFSCAN(I))THEN
            DO J=1,NCHAN
              IF(IFCHAN(J))THEN
                MAX=-1.E10
                MEAN=0.0
                SIG=0.0
                DO II=1,NF
                  MEAN=MEAN+A(J,II)
                  IF(A(J,II).GT.MAX) MAX=A(J,II)
                END DO
                MEAN=(MEAN-MAX)/(FILNUM-1.)
                DO II=1,NF
                  SIG=SIG+(A(J,II)-MEAN)*(A(J,II)-MEAN)
                END DO
                SIG=SIG-(MAX-MEAN)*(MAX-MEAN)
                IF(SIG.LE.0.)THEN
                  SIG=0.
                ELSE
                  SIG=SQRT(SIG/(FILNUM-2.))
                END IF
                DO II=1,NF
                  IF(SIG.LE.0.0)THEN
                    TSIGMA0=0.
                  ELSE
                    TSIGMA0=(A(J,II)-MEAN)/SIG
                  END IF
                  IF(TSIGMA0.LT.TSIGMAMIN) TSIGMAMIN=TSIGMA0
                  IF(TSIGMA0.GT.TSIGMAMAX) TSIGMAMAX=TSIGMA0
                END DO
              END IF                                                 !IFCHAN(J)
            END DO
          END IF                                                     !IFSCAN(I)
        END DO
C------------------------------------------------------------------------------
C Ya hemos leido todas las imagenes, asi que cerramos los ficheros
        DO II=1,NF
          CLOSE(II+10)
        END DO
C------------------------------------------------------------------------------
C Hacemos un histograma con los valores de sigma
20      WRITE(*,101)'* Histogram Npixels vs Sigma(N files)'
        WRITE(CDUMMY,*)TSIGMAMIN
        WRITE(*,100)'Minimum value (sigma units) '
        BMIN=READF(CDUMMY)
        WRITE(CDUMMY,*)TSIGMAMAX
        WRITE(*,100)'Maximum value (sigma units) '
        BMAX=READF(CDUMMY)
        WRITE(CDUMMY,*)NBIN
        WRITE(*,100)'No. of bins '
        NBIN=READILIM(CDUMMY,1,NBINMAX)
        WRITE(*,100)'Times sigma to estimate number of C.R. '
        TSLIMIT=READF('10.0')
        WRITE(*,100)'List possible C.R. (y/n) '
        CLIST(1:1)=READC('n','yn')
        NCR=0
        DO II=1,NF
          NCOS(II)=0
        END DO
C
        DO K=1,NBIN
          NY(K)=0
        END DO
        DB=(BMAX-BMIN)/REAL(NBIN)
C
C abrimos todos los ficheros otra vez
        DO II=1,NF
          FILE1=FILEN(II)
          WRITE(*,100)'* Opening: '
          WRITE(*,101)FILEN(I)(1:TRUELEN(FILEN(I)))
          FILENAME=
     +     INFILEX(II+10,FILE1,NSCAN,NCHAN,STWV,DISP,11,.FALSE.)
        END DO
C
        WRITE(*,*)
C
        DO I=1,NSCAN
          WRITE(*,'(A,I4,$)')'\b\b\b\b',I
          DO II=1,NF
            READ(II+10) (A(J,II),J=1,NCHAN)
          END DO
          IF(IFSCAN(I))THEN
            DO J=1,NCHAN
              IF(IFCHAN(J))THEN
                MEAN=0.
                MAX=0.
                SIG=0.
                DO II=1,NF
                  MEAN=MEAN+A(J,II)
                  IF(A(J,II).GT.MAX) MAX=A(J,II)
                END DO
                MEAN=(MEAN-MAX)/(FILNUM-1.)
                DO II=1,NF
                  SIG=SIG+(A(J,II)-MEAN)*(A(J,II)-MEAN)
                END DO
                SIG=SIG-(MAX-MEAN)*(MAX-MEAN)
                IF(SIG.LE.0.)THEN
                  SIG=0.
                ELSE
                  SIG=SQRT(SIG/(FILNUM-2.))
                END IF
                B(J,I)=MEAN
                SB(J,I)=SIG
                DO II=1,NF
                  IF(SIG.LE.0.0)THEN
                    TSIGMA0=0.
                  ELSE
                    TSIGMA0=(A(J,II)-MEAN)/SIG
                  END IF
                  IF(TSIGMA0.GT.TSLIMIT)THEN
                    NCR=NCR+1
                    IF(CLIST.EQ.'y')THEN
                      WRITE(*,113) J,I,TSIGMA0,II
                    END IF
                    NCOS(II)=NCOS(II)+1
                    IF(NCOS(II).GT.1000)THEN
                      WRITE(*,101)'No. cosmics > 1000!'
                      WRITE(*,*)
                      WRITE(*,100)'press <CR> to continue...'
                      READ(*,*)
                      DO III=1,NF
                        CLOSE(III+10)
                      END DO
                      GOTO 20
                    END IF
                    TIMES(II,NCOS(II))=TSIGMA0
                    XPOS(II,NCOS(II))=J
                    YPOS(II,NCOS(II))=I
                  END IF
                  K=NINT((TSIGMA0-BMIN)/DB)+1
                  IF((K.GE.1).AND.(K.LE.NBIN)) NY(K)=NY(K)+1
                END DO
              END IF
            END DO
          END IF
        END DO
        WRITE(*,*)
        DO II=1,NF
          CLOSE(II+10)
        END DO
C
        YMAX=0.
        YMIN=-0.3
        DO K=1,NBIN
          IF(NY(K).GT.0)THEN
            Y(K)=ALOG10(REAL(NY(K)))
            IF(Y(K).GT.YMAX) YMAX=Y(K)
          ELSE
            Y(K)=-1.
          END IF
        END DO
        DY=YMAX-YMIN
        YMAX=YMAX+DY/50.
C
        DO K=1,NBIN
          X(K)=REAL(K)
        END DO
C
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGENV(BMIN,BMAX,YMIN,YMAX,0,-2)
          IF(LCOLOR(ITERM)) CALL PGSCI(3)
          CALL PGSLS(2)
          CALL PGMOVE(TSLIMIT,YMIN)
          CALL PGDRAW(TSLIMIT,YMAX)
          CALL PGSLS(1)
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
          CALL PGBOX('BCNITS',0.0,0,'BCNITS',0.0,0)
          CALL PGWINDOW(1.,REAL(NBIN),YMIN,YMAX)
          CALL PGIDEN_RED
          IF(LCOLOR(ITERM)) CALL PGSCI(2)
          CALL PGBIN(NBIN,X,Y,.TRUE.)
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
          CALL PGLABEL('Deviations of the mean (in \\gs units)',
     +       'Log(N\\dpixels\\u)',LISTFILES)
        END DO
C
        WRITE(*,100)'>>> Estimated number of C.R.: '
        WRITE(*,*) NCR
        WRITE(*,100)'Replot (y/n) '
        CREPLOT(1:1)=READC('n','yn')
        IF(CREPLOT.EQ.'y') GOTO 20
        CALL PGEND
C------------------------------------------------------------------------------
        WRITE(*,100)'Output to a log file (y,n) '
        LOUT(1:1)=READC('y','yn')
        IF(LOUT.EQ.'y') THEN
          WRITE(*,100)'Log filename'
          LOGFILE=OUTFILEX(10,'@',0,0,0.,0.,3,.FALSE.)
        END IF
C
50      WRITE(*,100)'Times sigma to detect cosmics '
        WRITE(CDUMMY,*)TSLIMIT
        TSIGMA=READF(CDUMMY)
        IF(TSIGMA.LE.0.0) GOTO 50
C
51      WRITE(*,100)'Times sigma to detect cosmic tails '
        WRITE(CDUMMY,*)TSLIMIT
        TSIGMA2=READF(CDUMMY)
        IF(TSIGMA2.LE.0.0) GOTO 51
C------------------------------------------------------------------------------
C Abrimos de nuevo los ficheros (tanto de entrada como de salida)
        DO I=1,NF
          FILE1=FILEN(I)
          FILE2=FILE1(1:TRUELEN(FILE1))//'c'
          WRITE(*,100)'* Opening: '
          WRITE(*,101)FILEN(I)(1:TRUELEN(FILEN(I)))
          FILENAME=
     +     INFILEX(11,FILE1,NSCAN,NCHAN,STWV,DISP,11,.FALSE.)
          FILENAME=
     +     OUTFILEX(12,FILE2,NSCAN,NCHAN,STWV,DISP,11,.FALSE.)
          IF(CERR.EQ.'y')THEN
            CALL GUESSEF(FILE2,FILE3)
            FILENAME=
     +       OUTFILEX(13,FILE3,NSCAN,NCHAN,STWV,DISP,11,.TRUE.)
          END IF

          DO K=1,NSCAN
            READ(11) (A(J,K),J=1,NCHAN)
          END DO
          CLOSE(11)

          WRITE(*,101) FILEN(I)
          IF(LOUT.EQ.'y') WRITE(10,101) FILEN(I)
          NCOSF=0
          IF(NCOS(I).GT.0) THEN
            DO NC=1,NCOS(I)

              CALL TAILS(XPOS(I,NC),YPOS(I,NC),TSIGMA2,CRFOUND,NPIXCR)

              IF(CRFOUND) THEN
                NCOSF=NCOSF+1 
                WRITE(*,111) XPOS(I,NC),YPOS(I,NC),TIMES(I,NC),NPIXCR
                IF(LOUT.EQ.'y') WRITE(10,111) XPOS(I,NC),YPOS(I,NC),
     +             TIMES(I,NC),NPIXCR
              END IF
            END DO
          END IF
          WRITE(*,112) NCOSF  
          IF(LOUT.EQ.'y') WRITE(10,112) NCOSF

          DO K=1,NSCAN
            WRITE(12) (A(J,K),J=1,NCHAN)
          END DO
          CLOSE(12)
          IF(CERR.EQ.'y') CLOSE(13)
        END DO

        CLOSE(10)

        STOP
100     FORMAT(A,$)
101     FORMAT(A)
111     FORMAT('CR at pixel ',I5,I5,' (',F6.1,' sigmas)  ',I4,
     +     ' pixels')
112     FORMAT('Total No. CR:',I4) 
113     FORMAT('CR at pixel ',I5,I5,' (',F6.1,' sigmas), file#',
     +    I2.2)
        END
C
C******************************************************************************
C Busca alrededor del pixel sospechoso para detectar todos los pixels
C afectados por el rayo cosmico.
        SUBROUTINE TAILS(XP,YP,SIGMA,CRFOUND,NPIXCR)
C
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
C
        INTEGER XP,YP
        INTEGER II,JJ,NPIXCR,KK,LL,NPIX
        LOGICAL CRFOUND
        LOGICAL NEIGHBOR
        REAL SIGMA
        REAL DEV
        REAL A(NCMAX,NSMAX),B(NCMAX,NSMAX),SB(NCMAX,NSMAX)
        LOGICAL MASK(NCMAX,NSMAX),VUELTA

        COMMON/BLKAUTOCOS1/NSCAN,NCHAN
        COMMON/BLKAUTOCOS2/A,B,SB
        COMMON/BLKAUTOCOS3/MASK
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        CRFOUND=.FALSE.
        NPIXCR=0
C Pixel central
        DEV=(A(XP,YP)-B(XP,YP))/SB(XP,YP)
        IF(DEV.GT.SIGMA)THEN
          NPIXCR=1
          A(XP,YP)=B(XP,YP)
          DO II=1,NSCAN
            DO JJ=1,NCHAN
              MASK(JJ,II)=.FALSE.
            END DO
          END DO
          MASK(XP,YP)=.TRUE.
          CRFOUND=.TRUE.
        ELSE
          RETURN
        END IF
C Movimiento en coronas cuadradas de tamanho creciente
        NPIX=0
        KK=0
        VUELTA=.TRUE.
        DO WHILE(NPIX.LT.(NCHAN*NSCAN-1).AND.VUELTA)
          VUELTA=.FALSE.
          KK=KK+2
C posicion de arranque (sobre la diagonal)
          JJ=XP-KK/2
          II=YP-KK/2
C nos movemos hacia arriba 
          DO LL=1,KK
            II=II+1
            IF(JJ.GE.1.AND.JJ.LE.NCHAN.AND.II.GE.1.AND.II.LE.NSCAN)THEN
              NPIX=NPIX+1
              DEV=(A(JJ,II)-B(JJ,II))/SB(JJ,II)
              IF(DEV.GT.SIGMA)THEN
                IF(NEIGHBOR(JJ,II))THEN
                  MASK(JJ,II)=.TRUE.
                  A(JJ,II)=B(JJ,II)
                  NPIXCR=NPIXCR+1
                  VUELTA=.TRUE.
                END IF
              END IF
            END IF
          END DO
C nos movemos hacia la derecha
          DO LL=1,KK
            JJ=JJ+1
            IF(JJ.GE.1.AND.JJ.LE.NCHAN.AND.II.GE.1.AND.II.LE.NSCAN)THEN
              NPIX=NPIX+1
              DEV=(A(JJ,II)-B(JJ,II))/SB(JJ,II)
              IF(DEV.GT.SIGMA)THEN
                IF(NEIGHBOR(JJ,II))THEN
                  MASK(JJ,II)=.TRUE.
                  A(JJ,II)=B(JJ,II)
                  NPIXCR=NPIXCR+1
                  VUELTA=.TRUE.
                END IF
              END IF
            END IF
          END DO
C nos movemos hacia abajo
          DO LL=1,KK
            II=II-1
            IF(JJ.GE.1.AND.JJ.LE.NCHAN.AND.II.GE.1.AND.II.LE.NSCAN)THEN
              NPIX=NPIX+1
              DEV=(A(JJ,II)-B(JJ,II))/SB(JJ,II)
              IF(DEV.GT.SIGMA)THEN
                IF(NEIGHBOR(JJ,II))THEN
                  MASK(JJ,II)=.TRUE.
                  A(JJ,II)=B(JJ,II)
                  NPIXCR=NPIXCR+1
                  VUELTA=.TRUE.
                END IF
              END IF
            END IF
          END DO
C nos movemos hacia la izquierda, hasta la posicion de arranque
          DO LL=1,KK
            JJ=JJ-1
            IF(JJ.GE.1.AND.JJ.LE.NCHAN.AND.II.GE.1.AND.II.LE.NSCAN)THEN
              NPIX=NPIX+1
              DEV=(A(JJ,II)-B(JJ,II))/SB(JJ,II)
              IF(DEV.GT.SIGMA)THEN
                IF(NEIGHBOR(JJ,II))THEN
                  MASK(JJ,II)=.TRUE.
                  A(JJ,II)=B(JJ,II)
                  NPIXCR=NPIXCR+1
                  VUELTA=.TRUE.
                END IF
              END IF
            END IF
          END DO
        END DO

        END
C
C******************************************************************************
C Determina si algun pixel de alrededor ya ha sido marcado como posible C.R.
        LOGICAL FUNCTION NEIGHBOR(J,I)

        IMPLICIT NONE
        INTEGER J,I
        INCLUDE 'redlib.inc'
C
        INTEGER II,JJ
        LOGICAL MASK(NCMAX,NSMAX)
        COMMON/BLKAUTOCOS1/NSCAN,NCHAN
        COMMON/BLKAUTOCOS3/MASK
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        NEIGHBOR=.FALSE.
        DO II=I-1,I+1
          DO JJ=J-1,J+1
            IF((JJ.GE.1).AND.(JJ.LE.NCHAN).AND.(II.GE.1).AND.
     +      (II.LE.NSCAN))THEN
              IF((II.EQ.I).AND.(JJ.EQ.J))THEN
              ELSE
                IF(MASK(JJ,II))THEN
                  NEIGHBOR=.TRUE.
                  RETURN
                END IF
              END IF
            END IF
          END DO
        END DO

        END
