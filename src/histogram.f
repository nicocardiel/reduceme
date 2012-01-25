C------------------------------------------------------------------------------
C Version 23-April-1999                                       file: histogram.f
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
C Program: histogram
C Classification: examination & statistics
C Description: Performs a histogram of data tabulated in a single column
C ASCII file.
C
Comment
C
        PROGRAM HISTOGRAM
        IMPLICIT NONE
        INCLUDE 'futils.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'redlib.inc'
        INTEGER READILIM
        REAL READF
C
        INTEGER NMAXDATA
        PARAMETER (NMAXDATA=1000)
C
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
        CHARACTER*80 INFILE,CDUMMY
C------------------------------------------------------------------------------
        OUTFILEX=OUTFILEX
C
        THISPROGRAM='histogram'
        CALL WELCOME('23-April-1999')
C
        WRITE(*,100)'Input file name'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,3,.FALSE.)
        I=1
10      READ(20,*,END=20) A(I)
        I=I+1
        GOTO 10
20      CLOSE(20)
        NDAT=I-1
C------------------------------------------------------------------------------
        CALL FINDMM(NDAT,A,XMIN,XMAX)
        WRITE(*,100) 'Plot limits: [a]uto, [f]ixed (a/f) '
        CAUTO(1:1)=READC('a','af')
        IF(CAUTO.EQ.'a')THEN
          DX=XMAX-XMIN
          XMIN=XMIN-DX/20.
          XMAX=XMAX+DX/20.
        ELSE
          WRITE(CDUMMY,*) XMIN
          WRITE(*,100) 'Xmin '
          XMIN=READF(CDUMMY)
          WRITE(CDUMMY,*) XMAX
          WRITE(*,100) 'Xmax '
          XMAX=READF(CDUMMY)
        END IF
        WRITE(*,100) 'No. of bins '
        NBIN=READILIM('@',1,NMAXDATA)
        WIDTH=(XMAX-XMIN)/REAL(NBIN)
C------------------------------------------------------------------------------
        DO I=1,NBIN
          NFILL(I)=0
        END DO
C
        DO I=1,NDAT
          I0=INT((A(I)-XMIN)/WIDTH)+1
ccc          IF(I0.EQ.NBIN+1) I0=NBIN !incluimos en el ultimo bin ambos bordes
          IF((I0.GT.0).AND.(I0.LT.NBIN))THEN
            NFILL(I0)=NFILL(I0)+1
          END IF
        END DO
C
        DO I=1,NBIN
          Y(I)=REAL(NFILL(I))
        END DO
C------------------------------------------------------------------------------
        FMEAN=FMEAN0(NDAT,A,FSIGMA)
        PERCENT16=FPERCENT(NDAT,A,16.0)
        PERCENT50=FPERCENT(NDAT,A,50.0)
        PERCENT84=FPERCENT(NDAT,A,84.0)
C
        WRITE(*,100) 'Mean: '
        WRITE(*,*) FMEAN
        WRITE(*,100) 'rms.: '
        WRITE(*,*) FSIGMA
        WRITE(*,100) '16% percentile: '
        WRITE(*,*) PERCENT16
        WRITE(*,100) '50% percentile: '
        WRITE(*,*) PERCENT50
        WRITE(*,100) '84% percentile: '
        WRITE(*,*) PERCENT84
C------------------------------------------------------------------------------
        CALL FINDMM(NBIN,Y,YMIN,YMAX)
        DY=YMAX-YMIN
        YMIN=YMIN-DY/20.
        YMAX=YMAX+DY/20.
        WRITE(*,100) '>>> Ymin, Ymax: '
        WRITE(*,*) YMIN,YMAX
C
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
C------------------------------------------------------------------------------
        STOP
100     FORMAT(A,$)
        END
