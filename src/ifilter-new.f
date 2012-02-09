C------------------------------------------------------------------------------
C Version 13-May-2005                                           file: ifilter.f
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
C Program: ifilter
C Classification: arithmetic & manipulations
C Description: Applies a filter to an image.
C
Comment
C
C permite pasar un filtro a una imagen
C
        PROGRAM IFILTER
        IMPLICIT NONE
C Include NSMAX,NCMAX,NSCAN,NCHAN
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER READILIM
        REAL READF
C
        INTEGER NMAXBOX_X,NMAXBOX_Y
        PARAMETER(NMAXBOX_X=1025,NMAXBOX_Y=1025)
        INTEGER NMAXDEG
        PARAMETER (NMAXDEG=19) !maximum polynomial degree for POLFIT
C
        INTEGER I,J,K
        INTEGER II,JJ
        INTEGER IWIDTH,NDEGREE
        INTEGER IWIDTHMIN
        INTEGER NNS1,NNS2,NNC1,NNC2
        INTEGER NS1,NS2,NC1,NC2
        INTEGER NS1_,NS2_,NC1_,NC2_
        INTEGER NBOX_X,NBOX_Y,NBOX2_X,NBOX2_Y
        INTEGER NLIM1,NLIM2
        INTEGER NEXTINFO
        REAL A(NCMAX,NSMAX),B(NCMAX,NSMAX)
        REAL SPMASK(NCMAX)
        REAL PIXEL(NMAXBOX_X*NMAXBOX_Y)
        REAL FMEAN1,FMEAN2,FMEDIAN,FMEDIAN1
        REAL TSIGMA
        REAL XF(NMAXBOX_X),YF(NMAXBOX_X)
        REAL COEF(NMAXDEG+1),CHISQR
        REAL THRESHOLD
        CHARACTER*1 CMASK,CSAVE
        CHARACTER*1 CFILT,CMORE,CEXTEND,CWHOLE,CDIRECTION
        CHARACTER*75 INFILE,OUTFILE
        CHARACTER* 50 CDUMMY
        LOGICAL ANOMASK(NCMAX,NSMAX)
        LOGICAL LFIT(NMAXBOX_X)
        LOGICAL LOOP
C------------------------------------------------------------------------------
        THISPROGRAM='ifilter'
        CALL WELCOME('13-May-2005')
C
        THRESHOLD=0.0  !avoid compilation warning
        IWIDTH=3       !avoid compilation warning
        IWIDTHMIN=3    !avoid compilation warning
        CDIRECTION='@' !avoid compilation warning
C------------------------------------------------------------------------------
        WRITE(*,100) 'Input file name'
        INFILE=INFILEX(30,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        DO I=1,NSCAN
          READ(30) (A(J,I),J=1,NCHAN)
        END DO
        CLOSE(30)
C
        WRITE(*,100) 'Are you using a mask image (y/n) '
        CMASK(1:1)=READC('n','yn')
        IF(CMASK.EQ.'y')THEN
          WRITE(*,101) '> NOTE: mask pixels with signal != 0.0 '//
     +     'will be masked'
          WRITE(*,100) 'Input mask file name'
          INFILE=INFILEX(31,'@',NSCAN,NCHAN,STWV,DISP,21,.FALSE.)
          DO I=1,NSCAN
            READ(31) (SPMASK(J),J=1,NCHAN)
            DO J=1,NCHAN
              ANOMASK(J,I)=(SPMASK(J).EQ.0.0)
            END DO
          END DO
          CLOSE(31)
        ELSE
          DO I=1,NSCAN
            DO J=1,NCHAN
              ANOMASK(J,I)=.TRUE.
            END DO
          END DO
        END IF
C
C inicializamos la imagen final a la imagen inicial por si acaso solo filtramos
C un trozo de la imagen
        DO I=1,NSCAN
          DO J=1,NCHAN
            B(J,I)=A(J,I)
          END DO
        END DO
C
        CFILT='3' !opcion por defecto inicial
10      WRITE(*,101) '---------------'
        WRITE(*,101) '(1) MEAN filter'
        WRITE(*,101) '(2) MEAN filter excluding points'
        WRITE(*,101) '(3) MEDIAN filter'
        WRITE(*,101) '(4) 1-D Savitzky-Golay filter'
        WRITE(*,101) '(5) 1-D Savitzky-Golay filter (excluding points)'
        WRITE(*,101) '(6) Set to zero below threshold'
        WRITE(*,101) '(0) EXIT'
        WRITE(*,*)
        WRITE(*,100) 'Option '
        CFILT(1:1)=READC(CFILT,'0123456')
C
        IF(CFILT.EQ.'0')THEN
          STOP
        ELSEIF(CFILT.EQ.'2')THEN
          WRITE(*,100) 'Times sigma to exclude points '
          TSIGMA=READF('3.0')
        ELSEIF((CFILT.EQ.'4').OR.(CFILT.EQ.'5'))THEN
          WRITE(*,100) 'Choose direction (x/y)'
          CDIRECTION(1:1)=READC('@','xy')
          WRITE(*,100) 'Window width (pixels) '
          IWIDTH=READILIM('3',3,NMAXBOX_X)
          IF(CDIRECTION.EQ.'x')THEN
            IF(IWIDTH.GT.NCHAN)THEN
              WRITE(*,101) 'ERROR: window width larger than image size'
              WRITE(*,100) 'Press <CR> to continue...'
              READ(*,*)
              GOTO 10
            END IF
          ELSE
            IF(IWIDTH.GT.NSCAN)THEN
              WRITE(*,101) 'ERROR: window width larger than image size'
              WRITE(*,100) 'Press <CR> to continue...'
              READ(*,*)
              GOTO 10
            END IF
          END IF
          WRITE(*,100) 'Polynomial degree'
          NDEGREE=READILIM('@',0,NMAXDEG)
          IF(CFILT.EQ.'5')THEN
            WRITE(*,100) 'Times sigma to exclude points '
            TSIGMA=READF('3.0')
          END IF
          IF(CMASK.EQ.'y')THEN
            WRITE(*,*)
            WRITE(*,101) 'WARNING: Due to the use of the mask, '//
     +       'the running window will be stretched'
            WRITE(*,101) '         until reaching a mininum number '//
     +       'of pixels (if available).'
            WRITE(*,*)
            WRITE(*,100) 'Mininum number of unmasked pixels in '//
     +       'that case'
            WRITE(CDUMMY,*) IWIDTH
            IWIDTHMIN=READILIM(CDUMMY,NDEGREE+1,IWIDTH)
          END IF
        ELSEIF(CFILT.EQ.'6')THEN
          WRITE(*,100) 'Threshold'
          THRESHOLD=READF('@')
        END IF
C------------------------------------------------------------------------------
        WRITE(*,101) 'Define rectangle in image to aply filter:'
        CWHOLE='y'
        WRITE(*,100) 'first scan '
        NNS1=READILIM('1',1,NSCAN)
        WRITE(*,100) 'last scan '
        WRITE(CDUMMY,*) NSCAN
        NNS2=READILIM(CDUMMY,NNS1,NSCAN)
        IF((NNS1.NE.1).OR.(NNS2.NE.NSCAN)) CWHOLE='n'
C
        WRITE(*,100) 'first channel '
        NNC1=READILIM('1',1,NCHAN)
        WRITE(*,100) 'last channel '
        WRITE(CDUMMY,*) NCHAN
        NNC2=READILIM(CDUMMY,NNC1,NCHAN)
        IF((NNC1.NE.1).OR.(NNC2.NE.NCHAN)) CWHOLE='n'
C
        IF((CFILT.EQ.'4').OR.(CFILT.EQ.'5')) GOTO 40
        IF(CFILT.EQ.'6') GOTO 42
C
        IF(CWHOLE.EQ.'n')THEN
          WRITE(*,100)'Use data outside region to be filtered (y/n) '
          CEXTEND(1:1)=READC('y','yn')
        ELSE
          CEXTEND='y'
        END IF
C------------------------------------------------------------------------------
21      WRITE(*,100) 'Box  width: X-direction (must be odd)'
        NBOX_X=READILIM('@',1,MIN(NCHAN,NMAXBOX_X))
        IF(MOD(NBOX_X,2).EQ.0)THEN
          WRITE(*,101) 'ERROR: box width must be odd. Try again.'
          GOTO 21
        END IF
22      WRITE(*,100) 'Box height: Y-direction (must be odd)'
        NBOX_Y=READILIM('@',1,MIN(NSCAN,NMAXBOX_Y))
        IF(MOD(NBOX_Y,2).EQ.0)THEN
          WRITE(*,101) 'ERROR: box width must be odd. Try again.'
          GOTO 22
        END IF
C------------------------------------------------------------------------------
        NBOX2_X=NBOX_X/2
        NBOX2_Y=NBOX_Y/2
        DO I=NNS1,NNS2
          CALL SHOWPERC(NNS1,NNS2,1,I,NEXTINFO)
          NS1=I-NBOX2_Y
          NS2=I+NBOX2_Y
          IF(CEXTEND.EQ.'n')THEN
            NLIM1=NNS1
            NLIM2=NNS2
          ELSE
            NLIM1=1
            NLIM2=NSCAN
          END IF
          IF(NS1.LT.NLIM1) NS1=NLIM1
          IF(NS2.GT.NLIM2) NS2=NLIM2
          DO J=NNC1,NNC2
            NC1=J-NBOX2_X
            NC2=J+NBOX2_X
            IF(CEXTEND.EQ.'n')THEN
              NLIM1=NNC1
              NLIM2=NNC2
            ELSE
              NLIM1=1
              NLIM2=NCHAN
            END IF
            IF(NC1.LT.NLIM1) NC1=NLIM1
            IF(NC2.GT.NLIM2) NC2=NLIM2
            K=0
            DO II=NS1,NS2
              DO JJ=NC1,NC2
                IF(ANOMASK(JJ,II))THEN
                  K=K+1
                  PIXEL(K)=A(JJ,II)
                END IF
              END DO
            END DO
            IF(CFILT.EQ.'1')THEN                                   !mean filter
              B(J,I)=FMEAN1(K,PIXEL)
            ELSEIF(CFILT.EQ.'2')THEN          !mean filter rejecting with sigma
              B(J,I)=FMEAN2(K,PIXEL,TSIGMA)
            ELSEIF(CFILT.EQ.'3')THEN                             !median filter
              FMEDIAN=FMEDIAN1(K,PIXEL)
              B(J,I)=FMEDIAN
            ELSE
              STOP 'FATAL ERROR: invalid CFILT value.'
            END IF
          END DO
        END DO
        GOTO 50
C------------------------------------------------------------------------------
C Filtro de Savitzky-Golay
40      CONTINUE
        IF(CDIRECTION.EQ.'x')THEN
C..............................................................................
          DO I=NNS1,NNS2
            CALL SHOWPERC(NNS1,NNS2,1,I,NEXTINFO)
            DO J=NNC1,NNC2
              NC1=J-IWIDTH/2
              IF(NC1.LT.NNC1)THEN
                NC1=NNC1
                NC2=NC1+IWIDTH-1
                IF(NC2.GT.NNC2) STOP 'FATAL ERROR #1A'
              ELSE
                NC2=J+IWIDTH/2
                IF(NC2.GT.NNC2)THEN
                  NC2=NNC2
                  NC1=NC2-IWIDTH+1
                  IF(NC1.LT.NNC1) STOP 'FATAL ERROR #2A'
                END IF
              END IF
              K=0
              DO JJ=NC1,NC2
                IF(ANOMASK(JJ,I))THEN
                  K=K+1
                  XF(K)=REAL(JJ)
                  YF(K)=A(JJ,I)
                END IF
              END DO
              !codigo nuevo
              IF(K.LT.IWIDTHMIN)THEN
                !nos movemos un pixel a cada lado (si es posible) y vamos
                !incrementando el numero de pixeles utilizables
                NC1_=NC1
                NC2_=NC2
                LOOP=.TRUE.
                DO WHILE(LOOP)
                  IF(NC1_.GT.NNC1)THEN !podemos movernos 1 pixel a la izquierda
                    NC1_=NC1_-1
                    IF(ANOMASK(NC1_,I))THEN
                      K=K+1
                      XF(K)=REAL(NC1_)
                      YF(K)=A(NC1_,I)
                    END IF
                  END IF
                  IF(NC2_.LT.NNC2)THEN !podemos movernos 1 pixel a la derecha
                    NC2_=NC2_+1
                    IF(ANOMASK(NC2_,I))THEN
                      K=K+1
                      XF(K)=REAL(NC2_)
                      YF(K)=A(NC2_,I)
                    END IF
                  END IF
                  LOOP=(K.LT.IWIDTHMIN)
                  IF((NC1_.EQ.NNC1).AND.(NC2_.EQ.NNC2))THEN !ya no puede crecer
                    LOOP=.FALSE.
                  END IF
                END DO
              END IF
              !fin de codigo nuevo
              IF(K.GE.NDEGREE+1)THEN
                IF(CFILT.EQ.'4')THEN
                  CALL POLFIT(XF,YF,YF,K,NDEGREE+1,0,COEF,CHISQR)
                ELSE
                  CALL POLFITSIG(K,XF,YF,TSIGMA,NDEGREE,COEF,LFIT)
                END IF
                B(J,I)=COEF(NDEGREE+1)
                DO K=NDEGREE,1,-1
                  B(J,I)=B(J,I)*REAL(J)+COEF(K)
                END DO
              ELSE
                !STOP 'FATAL ERROR #3A'
              END IF
            END DO
          END DO
        ELSE
C..............................................................................
          DO J=NNC1,NNC2
            CALL SHOWPERC(NNC1,NNC2,1,I,NEXTINFO)
            DO I=NNS1,NNS2
              NS1=I-IWIDTH/2
              IF(NS1.LT.NNS1)THEN
                NS1=NNS1
                NS2=NS1+IWIDTH-1
                IF(NS2.GT.NNS2) STOP 'FATAL ERROR #1B'
              ELSE
                NS2=I+IWIDTH/2
                IF(NS2.GT.NNS2)THEN
                  NS2=NNS2
                  NS1=NS2-IWIDTH+1
                  IF(NS1.LT.NNS1) STOP 'FATAL ERROR #2B'
                END IF
              END IF
              K=0
              DO II=NS1,NS2
                IF(ANOMASK(J,II))THEN
                  K=K+1
                  XF(K)=REAL(II)
                  YF(K)=A(J,II)
                END IF
              END DO
              !codigo nuevo
              IF(K.LT.IWIDTHMIN)THEN
                !nos movemos un pixel a cada lado (si es posible) y vamos
                !incrementando el numero de pixeles utilizables
                NS1_=NS1
                NS2_=NS2
                LOOP=.TRUE.
                DO WHILE(LOOP)
                  IF(NS1_.GT.NNS1)THEN !podemos movernos 1 pixel hacia abajo
                    NS1_=NS1_-1
                    IF(ANOMASK(J,NS1_))THEN
                      K=K+1
                      XF(K)=REAL(NS1_)
                      YF(K)=A(J,NS1_)
                    END IF
                  END IF
                  IF(NS2_.LT.NNS2)THEN !podemos movernos 1 pixel hacia arriba
                    NS2_=NS2_+1
                    IF(ANOMASK(J,NS2_))THEN
                      K=K+1
                      XF(K)=REAL(NS2_)
                      YF(K)=A(J,NS2_)
                    END IF
                  END IF
                  LOOP=(K.LT.IWIDTHMIN)
                  IF((NS1_.EQ.NNS1).AND.(NS2_.EQ.NNS2))THEN !ya no puede crecer
                    LOOP=.FALSE.
                  END IF
                END DO
              END IF
              !fin de codigo nuevo
              IF(K.GE.NDEGREE+1)THEN
                IF(CFILT.EQ.'4')THEN
                  CALL POLFIT(XF,YF,YF,K,NDEGREE+1,0,COEF,CHISQR)
                ELSE
                  CALL POLFITSIG(K,XF,YF,TSIGMA,NDEGREE,COEF,LFIT)
                END IF
                B(J,I)=COEF(NDEGREE+1)
                DO K=NDEGREE,1,-1
                  B(J,I)=B(J,I)*REAL(I)+COEF(K)
                END DO
              ELSE
                !STOP 'FATAL ERROR #3A'
              END IF
            END DO
          END DO
        END IF
        GOTO 50
C------------------------------------------------------------------------------
C Set to zero below threshold
42      CONTINUE
        DO I=NNS1,NNS2
          CALL SHOWPERC(NNS1,NNS2,1,I,NEXTINFO)
          DO J=NNC1,NNC2
            IF(A(J,I).LT.THRESHOLD)THEN
              B(J,I)=0.0
            END IF
          END DO
        END DO
        GOTO 50
C------------------------------------------------------------------------------
50      CONTINUE
C
        WRITE(*,100) 'Are you saving the result (y/n) '
        CSAVE(1:1)=READC('y','yn')
        IF(CSAVE.EQ.'y')THEN
          WRITE(*,100) 'Output file name'
          OUTFILE=OUTFILEX(40,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
          DO I=1,NSCAN
            WRITE(40) (B(J,I),J=1,NCHAN)
          END DO
          CLOSE(40)
        END IF
C
        WRITE(*,100) 'More filters (y/n) '
        CMORE(1:1)=READC('y','yn')
        IF(CMORE.EQ.'y')THEN
          DO I=1,NSCAN
            DO J=1,NCHAN
              A(J,I)=B(J,I)
            END DO
          END DO
          GOTO 10
        END IF
C------------------------------------------------------------------------------
        STOP
100     FORMAT(A,$)
101     FORMAT(A)
C
        END
