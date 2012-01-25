C------------------------------------------------------------------------------
C Version 13-October-2007                                   file: fftcorrzoom.f
C------------------------------------------------------------------------------
C Copyright N. Cardiel & J. Gorgas, Departamento de Astrofisica
C Universidad Complutense de Madrid, 28040-Madrid, Spain
C E-mail: ncl@astrax.fis.ucm.es or fjg@astrax.fis.ucm.es
C------------------------------------------------------------------------------
C This routine is free software; you can redistribute it and/or modify it
C under the terms of the GNU General Public License as published by the Free
C Software Foundation; either version 2 of the License, or (at your option) any
C later version. See the file gnu-public-license.txt for details.
C------------------------------------------------------------------------------
Comment
C
C SUBROUTINE FFTCORRZOOM(X,Y,N,NL,NAME1,NAME2,NFIT,LPLOT,X0,Y0)
C
C Input: X,Y,N,NL,NAME1,NAME2,NFIT,LPLOT
C Output: X0,Y0
C
C This subroutine plots the correlation function given by X(N),Y(N), zooms in
C around the maximum and fit a second-order polynomial in order to obtain the
C exact location of the peak. The fitted value is returned through X0,Y0.
C
C REAL X(N) -> X-coordinates of the correlation function
C REAL Y(N) -> Y-coordinates of the correlation function
C INTEGER N -> dimension of X and N
C INTEGER NL -> number of pixels affected by zero padding
C CHARACTER*(*) NAME1 -> description of the first data set employed to compute
C                        cross correlation
C CHARACTER*(*) NAME2 -> description of the second data set employed to compute
C                        cross correlation
C INTEGER NFIT -> no. of points around peak to fit maximum (if NFIT=0, the
C                 routine asks for this number; otherwise the routine returns 
C                 without prompting)
C LOGICAL LPLOT -> if .TRUE., plot zoomed region
C REAL X0 -> X-offset of the peak of the correlation function
C REAL Y0 -> peak value of the correlation function
C
Comment
C------------------------------------------------------------------------------
        SUBROUTINE FFTCORRZOOM(X,Y,N,NL,NAME1,NAME2,NFIT,LPLOT,X0,Y0)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'futils.inc'
        REAL READF
        INTEGER READILIM
        INTEGER TRUELEN
C
        INTEGER N,NL
        REAL X(N),Y(N)
        CHARACTER*(*) NAME1,NAME2
        INTEGER NFIT
        LOGICAL LPLOT
        REAL X0,Y0
C
        INTEGER I,IMAX,L
        INTEGER NSIDE,NCOLOR,NF
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        REAL XL
        REAL YMIN,YMAX,DY
        REAL FMAX
        REAL XPOL(NCMAX),YPOL(NCMAX)
        REAL A(3),CHISQR
        CHARACTER*1 CCHAN
        CHARACTER*50 CDUMMY
        LOGICAL LCOLOR(MAX_ID_RED)
C
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        XL=REAL(NL)
        NCOLOR=1
        NSIDE=3
C
10      IF(NFIT.EQ.0)THEN
          WRITE(*,100)'Maximum X-shift around 0 '
          WRITE(CDUMMY,*)NL
          XL=READF(CDUMMY)
          NL=NINT(XL)
        END IF
C
        IF(LPLOT)THEN
          YMIN=1E30
          YMAX=-1E30
          DO I=1,N
            IF(ABS(X(I)).LE.XL)THEN
              IF(Y(I).LT.YMIN) YMIN=Y(I)
              IF(Y(I).GT.YMAX) YMAX=Y(I)
            END IF
          END DO
          DY=YMAX-YMIN
          YMIN=YMIN-DY/50.
          YMAX=YMAX+DY/10.
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL AUTOPLOT(N,X,Y,1,N,
     +       'X-offset','correlation value',
     +        'Correlation of the two data sets',
     +       .FALSE.,-XL,XL,YMIN,YMAX,0.00,
     +       0,.TRUE.,'BCNTS','BCNTS',
     +       101,5,
     +       0.,1.,0.0,0.5)
            CALL PGIDEN_RED
            CALL PGMTEXT('T',1.5,0.0,0.0,NAME1)
            CALL PGMTEXT('T',1.5,1.0,1.0,NAME2)
          END DO
        END IF
C
        IF(NFIT.EQ.0)THEN
          WRITE(*,100)'Change zoom (y/n) '
          CCHAN(1:1)=READC('n','yn')
          IF(CCHAN.EQ.'y') GOTO 10
        END IF
C
        FMAX=-1E30
        IMAX=0
        DO I=1,N
          IF(ABS(X(I)).LE.REAL(NL))THEN  !evitamos buscar el maximo fuera de NL
            IF(Y(I).GT.FMAX)THEN
              FMAX=Y(I)
              IMAX=I
            END IF
          END IF
        END DO
        IF(IMAX.EQ.0) STOP 'Unexpected IMAX=0 in FFTCORRZOOM'
C
20      IF(NFIT.EQ.0)THEN
          WRITE(*,100)'No. of point at each side of maximum to fit '//
     +     'polynomial '
          WRITE(CDUMMY,*)NSIDE
          NSIDE=READILIM(CDUMMY,1,(NCMAX-1)/2)
        ELSE
          NSIDE=NFIT
        END IF
C
        NF=2*NSIDE+1
        DO I=1,NF
          XPOL(I)=X(IMAX-NSIDE+I-1)
          YPOL(I)=Y(IMAX-NSIDE+I-1)
        END DO
C
        IF(LPLOT)THEN
          NCOLOR=NCOLOR+1
          IF(NCOLOR.GT.15) NCOLOR=2
          IF(NCOLOR.EQ.5) NCOLOR=6
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGSCI(NCOLOR)
            CALL PGBIN(NF,XPOL,YPOL,.TRUE.)
          END DO
        END IF
C
        CALL POLFIT(XPOL,YPOL,YPOL,NF,3,0,A,CHISQR)
        X0=-A(2)/(2.0*A(3))
        Y0=A(1)+A(2)*X0+A(3)*X0*X0
C
        IF(LPLOT)THEN
          WRITE(*,100)'Maximum is located at (x,y) = ('
          WRITE(CDUMMY,*)X0
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          WRITE(CDUMMY(L+1:),'(A1)')','
          WRITE(CDUMMY(L+2:),*)Y0
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          WRITE(*,100)CDUMMY(1:L)//') for '
          WRITE(*,101) NAME2(1:TRUELEN(NAME2))
        END IF
C
        IF(LPLOT)THEN
          DO I=1,NCMAX
            XPOL(I)=X(IMAX-NSIDE)+REAL(2*NSIDE)*REAL(I-1)/NCMAX
            YPOL(I)=A(1)+A(2)*XPOL(I)+A(3)*XPOL(I)*XPOL(I)
          END DO
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGLINE(NCMAX,XPOL,YPOL)
            CALL PGSLS(2)
            CALL PGMOVE(X0,YMIN)
            CALL PGDRAW(X0,YMAX)
            CALL PGSLS(1)
            CALL PGSCI(1)
          END DO
        END IF
C
        IF(NFIT.EQ.0)THEN
          WRITE(*,100)'Other fit (y/n) '
          CCHAN(1:1)=READC('n','yn')
          IF(CCHAN.EQ.'y') GOTO 20
        END IF
C
100     FORMAT(A,$)
101     FORMAT(A)
        END
