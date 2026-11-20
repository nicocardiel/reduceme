!------------------------------------------------------------------------------
! Version 20-November-2026                                    file: fitdata.f90
!------------------------------------------------------------------------------
! Copyright N. Cardiel & J. Gorgas, Departamento de Astrofisica
! Universidad Complutense de Madrid, 28040-Madrid, Spain
! E-mail: ncl@astrax.fis.ucm.es or fjg@astrax.fis.ucm.es
!------------------------------------------------------------------------------
! This program is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the Free
! Software Foundation; either version 2 of the License, or (at your option) any
! later version. See the file gnu-public-license.txt for details.
!------------------------------------------------------------------------------
!Comment
!
! Program: fitdata
! Classification: miscellany
! Description: Fits polynomials using tabulated data.
!
!Comment
!
! Programa para representar y ajustar polinomios a datos leidos de ficheros
! con columnas
!
        PROGRAM FITDATA
        USE RED_M
        USE FSPEC_M
        USE FUTILS_M
        USE REDLIB_M
        IMPLICIT NONE
!
        INCLUDE 'redlib.inc'
!
        INTEGER NMAXP
        PARAMETER (NMAXP=1000)                         !numero maximo de puntos
        INTEGER NMAXCOL
        PARAMETER (NMAXCOL=10)         !numero maximo de columnas en el fichero
!
        INTEGER K,I,IMIN
        INTEGER NP,NC,NCOL,NC1,NC2
        INTEGER NTERMS
        INTEGER NSYMB
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        REAL DATA(NMAXP,NMAXCOL)
        REAL XP(NMAXP),YP(NMAXP),YSIGMA(NMAXP)
        REAL XX(NMAXP),YY(NMAXP)
        REAL XMIN,XMAX,YMIN,YMAX,DX,DY
        REAL MEAN,SIGMA
        REAL A(20),X,POL,CHISQR
        REAL R,RMIN,XC,YC
        CHARACTER*1 CPLOT,CFIT,CLOCATE,CH
        CHARACTER*50 CDUMMY
        CHARACTER*75 CXLABEL,CYLABEL,CGLABEL
        CHARACTER*255 INFILE
        LOGICAL LCOLOR(MAX_ID_RED)
!------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
!
        THISPROGRAM='fitdata'
        CALL WELCOME('20-November-2026')
!
        CALL PIDEGTER(NTERM,IDN,LCOLOR)
!
        WRITE(*,100,ADVANCE='NO')'Input file name'
        INFILE=INFILEX(20,'@',0,0,0.,0.,3,.FALSE.)
        WRITE(*,100,ADVANCE='NO')'Number of columns '
        NCOL=READILIM('@',1,9999)
!
        NP=0
10      READ(20,*,END=12) (DATA(NP+1,NC),NC=1,NCOL)
        NP=NP+1
        WRITE(*,*) (DATA(NP,NC),NC=1,NCOL)
        GOTO 10
12      CLOSE(20)
        WRITE(*,110)'Number of rows: ',NP
!
        WRITE(*,100,ADVANCE='NO')'X column is #'
        NC1=READILIM('@',1,NCOL)
        WRITE(*,100,ADVANCE='NO')'Y column is #'
        NC2=READILIM('@',1,NCOL)
!
        DO K=1,NP
          XP(K)=DATA(K,NC1)
          YP(K)=DATA(K,NC2)
          YSIGMA(K)=0.
        END DO
!
        CALL FINDMM(NP,XP,XMIN,XMAX)
        CALL FINDMM(NP,YP,YMIN,YMAX)
!
        WRITE(*,100,ADVANCE='NO')'Xmin: '
        WRITE(*,*) XMIN
        WRITE(*,100,ADVANCE='NO')'Xmax: '
        WRITE(*,*) XMAX
        WRITE(*,100,ADVANCE='NO')'Ymin: '
        WRITE(*,*) YMIN
        WRITE(*,100,ADVANCE='NO')'Ymax: '
        WRITE(*,*) YMAX
!
        CALL STATISTICS(NP,XP,MEAN,SIGMA)
        WRITE(*,100,ADVANCE='NO')'X Mean..............: '
        WRITE(*,*)MEAN
        WRITE(*,100,ADVANCE='NO')'X standard deviation: '
        WRITE(*,*)SIGMA
        CALL STATISTICS(NP,YP,MEAN,SIGMA)
        WRITE(*,100,ADVANCE='NO')'Y Mean..............: '
        WRITE(*,*)MEAN
        WRITE(*,100,ADVANCE='NO')'Y standard deviation: '
        WRITE(*,*)SIGMA
!
        WRITE(*,100,ADVANCE='NO')'Plot data (y/n) '
        CPLOT(1:1)=READC('y','yn')
        IF(CPLOT.EQ.'n') STOP
!
        DX=XMAX-XMIN
        DY=YMAX-YMIN
        XMIN=XMIN-DX/50.
        XMAX=XMAX+DX/50.
        YMIN=YMIN-DY/50.
        YMAX=YMAX+DY/50.
!
        WRITE(CDUMMY,*)XMIN
        WRITE(*,100,ADVANCE='NO')'Xmin '
        XMIN=READF(CDUMMY)
        WRITE(CDUMMY,*)XMAX
        WRITE(*,100,ADVANCE='NO')'Xmax '
        XMAX=READF(CDUMMY)
        WRITE(CDUMMY,*)YMIN
        WRITE(*,100,ADVANCE='NO')'Ymin '
        YMIN=READF(CDUMMY)
        WRITE(CDUMMY,*)YMAX
        WRITE(*,100,ADVANCE='NO')'Ymax '
        YMAX=READF(CDUMMY)
!
        WRITE(*,100,ADVANCE='NO')'PGPLOT symbol number '
        NSYMB=READILIM('@',-1,31)
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGSCH(1.2)
          CALL PGSLW(3)
          CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,0)
          CALL PGIDEN_RED
          CALL PGPOINT(NP,XP,YP,NSYMB)
        END DO
        WRITE(*,100,ADVANCE='NO')'X -  label'
        CXLABEL(1:75)=READC('@','@')
        WRITE(*,100,ADVANCE='NO')'Y -  label'
        CYLABEL(1:75)=READC('@','@')
        WRITE(*,100,ADVANCE='NO')'Plot label'
        CGLABEL(1:75)=READC('@','@')
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGLABEL(CXLABEL,CYLABEL,CGLABEL)
        END DO
!
20      WRITE(*,100,ADVANCE='NO')'Locate data points with mouse (y/n) '
        CLOCATE(1:1)=READC('y','yn')
        IF(CLOCATE.EQ.'y')THEN
          CALL PGCURSE(XC,YC,CH)
          RMIN=(XC-XP(1))*(XC-XP(1))+(YC-YP(1))*(YC-YP(1))
          IMIN=1
          DO I=2,NP
            R=(XC-XP(I))*(XC-XP(I))+(YC-YP(I))*(YC-YP(I))
            IF(R.LT.RMIN)THEN
              RMIN=R
              IMIN=I
            END IF
          END DO
          CALL PGSCI(3)
          CALL PGPOINT(1,XP(IMIN),YP(IMIN),NSYMB)
          CALL PGSCI(1)
          WRITE(*,100,ADVANCE='NO')'Nearest point is #,X,Y: '
          WRITE(*,*)IMIN,XP(IMIN),YP(IMIN)
          WRITE(*,100,ADVANCE='NO')'Press <CR> to continue...'
          READ(*,*)
          CALL PGPOINT(1,XP(IMIN),YP(IMIN),NSYMB)
          GOTO 20
        END IF
!
        WRITE(*,100,ADVANCE='NO')'Fit polynomial (y/n) '
        CFIT(1:1)=READC('y','yn')
        IF(CFIT.EQ.'y')THEN
          WRITE(*,100,ADVANCE='NO')'Polynomial degree '
          NTERMS=READILIM('@',0,19)
          NTERMS=NTERMS+1
          CALL POLFIT(XP,YP,YP,NP,NTERMS,0,A,CHISQR)
          DO K=1,NTERMS
            IF(K.LT.11)THEN
              WRITE(*,'(A,I1,A)',ADVANCE='NO')'> a(',K-1,') : '
            ELSE
              WRITE(*,'(A,I2,A)',ADVANCE='NO')'> a(',K-1,'): '
            END IF
            WRITE(*,*) A(K)
          END DO
          DO I=1,NMAXP
            X=REAL(I-1)/REAL(NMAXP-1)*(XMAX-XMIN)+XMIN
            POL=A(NTERMS)
            DO K=NTERMS-1,1,-1
              POL=POL*X+A(K)
            END DO
            XX(I)=X
            YY(I)=POL
          END DO
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(2)
            CALL PGLINE(NMAXP,XX,YY)
          END DO
        END IF
!
        CALL PGEND
!
        STOP
100     FORMAT(A)
110     FORMAT(A,I6)
        END
!
!******************************************************************************
!
        SUBROUTINE STATISTICS(NP,X,MEAN,SIGMA)
        IMPLICIT NONE
        INTEGER NP
        REAL X(NP)
        REAL MEAN,SIGMA
!
        INTEGER K
        DOUBLE PRECISION DMEAN,DSIGMA
!------------------------------------------------------------------------------
        DMEAN=0.D0
        DO K=1,NP
          DMEAN=DMEAN+DBLE(X(K))
        END DO
        DMEAN=DMEAN/DBLE(NP)
        DSIGMA=0.D0
        IF(NP.GT.1)THEN
          DO K=1,NP
            DSIGMA=DSIGMA+(DBLE(X(K))-DMEAN)*(DBLE(X(K))-DMEAN)
          END DO
          DSIGMA=DSQRT(DSIGMA/DBLE(NP-1))
        END IF
        MEAN=REAL(DMEAN)
        SIGMA=REAL(DSIGMA)
        END
