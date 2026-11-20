!------------------------------------------------------------------------------
! Version 20-November-2026                                    file: histogram.f
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
! Program: histogram
! Classification: examination & statistics
! Description: Performs a histogram of data tabulated in a single column
! ASCII file.
!
!Comment
!
        PROGRAM HISTOGRAM
        IMPLICIT NONE
        INCLUDE 'futils.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'redlib.inc'
        INTEGER READILIM
        REAL READF
!
        INTEGER NMAXDATA
        PARAMETER (NMAXDATA=1000)
!
        INTEGER I,I0,NDAT
        INTEGER NBIN
        INTEGER NFILL(NMAXDATA)
        REAL A(NMAXDATA)
        REAL Y(NMAXDATA)
        REAL XMIN,XMAX,YMIN,YMAX,DX,DY
        REAL WIDTH
        REAL X1,X2,Y1,Y2
        REAL FMEAN0
        REAL FMEAN,FSIGMA
        REAL PERCENT16,PERCENT50,PERCENT84
        REAL FPERCENT
        CHARACTER*1 CAUTO
        CHARACTER*80 CDUMMY
        CHARACTER*255 INFILE
!------------------------------------------------------------------------------
        OUTFILEX=OUTFILEX
!
        THISPROGRAM='histogram'
        CALL WELCOME('20-November-2026')
!
        WRITE(*,100,ADVANCE='NO')'Input file name'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,3,.FALSE.)
        I=1
10      READ(20,*,END=20) A(I)
        I=I+1
        GOTO 10
20      CLOSE(20)
        NDAT=I-1
!------------------------------------------------------------------------------
        CALL FINDMM(NDAT,A,XMIN,XMAX)
        WRITE(*,100,ADVANCE='NO') 'Plot limits: [a]uto, [f]ixed (a/f) '
        CAUTO(1:1)=READC('a','af')
        IF(CAUTO.EQ.'a')THEN
          DX=XMAX-XMIN
          XMIN=XMIN-DX/20.
          XMAX=XMAX+DX/20.
        ELSE
          WRITE(CDUMMY,*) XMIN
          WRITE(*,100,ADVANCE='NO') 'Xmin '
          XMIN=READF(CDUMMY)
          WRITE(CDUMMY,*) XMAX
          WRITE(*,100,ADVANCE='NO') 'Xmax '
          XMAX=READF(CDUMMY)
        END IF
        WRITE(*,100,ADVANCE='NO') 'No. of bins '
        NBIN=READILIM('@',1,NMAXDATA)
        WIDTH=(XMAX-XMIN)/REAL(NBIN)
!------------------------------------------------------------------------------
        DO I=1,NBIN
          NFILL(I)=0
        END DO
!
        DO I=1,NDAT
          I0=INT((A(I)-XMIN)/WIDTH)+1
!cc          IF(I0.EQ.NBIN+1) I0=NBIN !incluimos en el ultimo bin ambos bordes
          IF((I0.GT.0).AND.(I0.LT.NBIN))THEN
            NFILL(I0)=NFILL(I0)+1
          END IF
        END DO
!
        DO I=1,NBIN
          Y(I)=REAL(NFILL(I))
        END DO
!------------------------------------------------------------------------------
        FMEAN=FMEAN0(NDAT,A,FSIGMA)
        PERCENT16=FPERCENT(NDAT,A,16.0)
        PERCENT50=FPERCENT(NDAT,A,50.0)
        PERCENT84=FPERCENT(NDAT,A,84.0)
!
        WRITE(*,100,ADVANCE='NO') 'Mean: '
        WRITE(*,*) FMEAN
        WRITE(*,100,ADVANCE='NO') 'rms.: '
        WRITE(*,*) FSIGMA
        WRITE(*,100,ADVANCE='NO') '16% percentile: '
        WRITE(*,*) PERCENT16
        WRITE(*,100,ADVANCE='NO') '50% percentile: '
        WRITE(*,*) PERCENT50
        WRITE(*,100,ADVANCE='NO') '84% percentile: '
        WRITE(*,*) PERCENT84
!------------------------------------------------------------------------------
        CALL FINDMM(NBIN,Y,YMIN,YMAX)
        DY=YMAX-YMIN
        YMIN=YMIN-DY/20.
        YMAX=YMAX+DY/20.
        WRITE(*,100,ADVANCE='NO') '>>> Ymin, Ymax: '
        WRITE(*,*) YMIN,YMAX
!
        CALL PGBEGIN(0,'?',1,1)
        CALL PGSLW(3)
        CALL PGSCF(2)
        CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,0)
        CALL PGLABEL('DATA value','No. of data',' ')
        CALL PGIDEN_RED
        DO I=1,NBIN
          X1=XMIN+REAL(I-1)*WIDTH
          X2=X1+WIDTH
          Y1=0.
          Y2=Y(I)
          CALL PGMOVE(X1,Y1)
          CALL PGDRAW(X1,Y2)
          CALL PGDRAW(X2,Y2)
          CALL PGDRAW(X2,Y1)
        END DO
        CALL PGEND
!------------------------------------------------------------------------------
        STOP
100     FORMAT(A)
        END
