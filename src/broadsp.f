C------------------------------------------------------------------------------
C Version 28-November-1996                                      file: broadsp.f
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
C Program: broadsp
C Classification: arithmetic & manipulations
C Description: Broaden a single spectrum by convolving with a gaussian. If the
C input file is an image, the program extracts a single spectrum prior
C broadening.
C
Comment
C
C Ensancha espectros introduciendo una dispersion de velocidades.
C
        PROGRAM BROADSP
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER TRUELEN
        REAL READF
C
        REAL C
        PARAMETER (C=299792.458)
C
        INTEGER I,J
        INTEGER N1,N2
        INTEGER NADD,NCOLOR
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        INTEGER NSCAN_,NCHAN_
        REAL STWV_,DISP_
        REAL A(NCMAX,NSMAX),ERR(NCMAX,NSMAX)
        REAL S(NCMAX),SS(NCMAX),X(NCMAX)
        REAL ES(NCMAX),ESS(NCMAX)
        REAL SIGMA,SIGMAVEC(NCMAX)
        REAL XMIN,XMAX,YMIN,YMAX,DX,DY
        REAL FLDO
        CHARACTER*1 CSAVE,CERR,CTYPEERR
        CHARACTER*75 INFILE,OUTFILE,ERRFILE
        CHARACTER*255 POLFILE
        LOGICAL IFSCAN(NSMAX)
        LOGICAL LCOLOR(MAX_ID_RED)
C------------------------------------------------------------------------------
        THISPROGRAM='broadsp'
        CALL WELCOME('28-November-1996')
        CALL SHOWHLP('explanation')
C
        WRITE(*,100)'Work with error images (y/n) '
        CERR(1:1)=READC('n','yn')
C Salida grafica
        CALL PIDEGTER(NTERM,IDN,LCOLOR)
C------------------------------------------------------------------------------
C Fichero de trabajo
        WRITE(*,100)'Input file name'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        DO I=1,NSCAN
          READ(20) (A(J,I),J=1,NCHAN)
        END DO
        CLOSE(20)
        IF(CERR.EQ.'y')THEN
          WRITE(*,100)'Input error file name '
          CALL GUESSEF(INFILE,ERRFILE)
          ERRFILE=INFILEX(20,ERRFILE,NSCAN,NCHAN,STWV,DISP,21,.TRUE.)!....match
          DO I=1,NSCAN
            READ(20) (ERR(J,I),J=1,NCHAN)
          END DO
          CLOSE(20)
        END IF
C------------------------------------------------------------------------------
        DO J=1,NCHAN
          X(J)=REAL(J)
        END DO
C------------------------------------------------------------------------------
C Si solo hay un espectro, sera el espectro de trabajo
        IF(NSCAN.EQ.1)THEN
          DO J=1,NCHAN
            S(J)=A(J,1)
          END DO
          IF(CERR.EQ.'y')THEN
            DO J=1,NCHAN
              ES(J)=ERR(J,1)
            END DO
          END IF
          NADD=1
          GOTO 30
        END IF
C------------------------------------------------------------------------------
C Si es una imagen, seleccionamos espectro
        DO I=1,NSCAN
          IFSCAN(I)=.FALSE.
        END DO
        WRITE(*,101)'* Define scans to be added: '
10      WRITE(*,100)'Scans (0,0=EXIT) '
        CALL READ2I('0,0',N1,N2)
        IF((N1.EQ.0).AND.(N2.EQ.0)) GOTO 20
        IF((N1.LT.1).OR.(N2.GT.NSCAN).OR.(N1.GT.N2))THEN
          WRITE(*,101)'ERROR: numbers out of range. Try again.'
          GOTO 10
        END IF
        DO I=N1,N2
          IFSCAN(I)=.TRUE.
        END DO
        GOTO 10
20      NADD=0
        DO I=1,NSCAN
          IF(IFSCAN(I)) NADD=NADD+1
        END DO
        IF(NADD.EQ.0)THEN
          WRITE(*,101)'ERROR: number of scans added = 0!'
          GOTO 10
        END IF
C
        DO J=1,NCHAN
          S(J)=0.
          DO I=1,NSCAN
            IF(IFSCAN(I)) S(J)=S(J)+A(J,I)
          END DO
          S(J)=S(J)/REAL(NADD)
        END DO
        IF(CERR.EQ.'y')THEN
          DO J=1,NCHAN
            ES(J)=0.
            DO I=1,NSCAN
              IF(IFSCAN(I)) ES(J)=ES(J)+ERR(J,I)
            END DO
            ES(J)=ES(J)/REAL(NADD)
          END DO
        END IF
C------------------------------------------------------------------------------
30      XMIN=1.
        XMAX=REAL(NCHAN)
        CALL FINDMM(NCHAN,S,YMIN,YMAX)
        DX=XMAX-XMIN
        DY=YMAX-YMIN
        XMIN=XMIN-DX/50.
        XMAX=XMAX+DX/50.
        YMIN=YMIN-DY/50.
        YMAX=YMAX+DY/50.
C
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,0)
          CALL PGIDEN_RED
          CALL PGBIN(NCHAN,X,S,.TRUE.)
          IF(NADD.EQ.1)THEN
            CALL PGLABEL('channel','No. of counts',CHAR(32))
          ELSE
            CALL PGLABEL('channel','Averaged No. of counts',CHAR(32))
          END IF
          CALL PGMTEXT('T',1.5,0.,0.,'File: '//INFILE)
          CALL PGMTEXT('T',1.5,1.,1.,OBJECT(1:TRUELEN(OBJECT)))
        END DO
C
        IF(CERR.EQ.'y')THEN
          CALL SHOWHLP('error types')
          WRITE(*,101)'1- compute error spectrum like a normal spectrum'
          WRITE(*,101)'2- compute error of the broadened spectrum'
          WRITE(*,100)'Option (1/2) '
          CTYPEERR(1:1)=READC('1','12')
        END IF
C
        NCOLOR=1
50      WRITE(*,100)'Sigma (km/s) (-1=from external file,0=EXIT): '
        SIGMA=READF('@')
        IF(SIGMA.EQ.0.0)THEN
          CALL PGEND
          STOP
        ELSEIF(SIGMA.GT.0.0)THEN
          DO J=1,NCHAN
            SIGMAVEC(J)=SIGMA
          END DO
        ELSE
          WRITE(*,100) 'File name with spectrum with FWHM (Angs)'
          WRITE(*,100) ' for each channel'
          POLFILE=INFILEX(30,'@',NSCAN_,NCHAN_,STWV_,DISP_,1,.FALSE.)
          IF(NSCAN_.NE.1)THEN
            WRITE(*,101) 'ERROR: NSCAN.NE.1'
            CLOSE(30)
            STOP
          ELSEIF(NCHAN_.NE.NCHAN)THEN
            WRITE(*,101) 'ERROR: NCHAN_.NE.NCHAN'
            CLOSE(30)
            STOP
          END IF
          READ(30) (SIGMAVEC(J),J=1,NCHAN)
          CLOSE(30)
          DO J=1,NCHAN
            SIGMAVEC(J)=SIGMAVEC(J)/2.35 !FWHM (Angs.) -> sigma (Angs.)
            FLDO=REAL(J-1)*DISP+STWV
            SIGMAVEC(J)=SIGMAVEC(J)/FLDO*C !sigma (Angs.) -> sigma (km/s)
          END DO
        END IF
        CALL BROADEN(S,SS,NCHAN,STWV,DISP,SIGMAVEC,.FALSE.)
        IF(CERR.EQ.'y')THEN
          IF(CTYPEERR.EQ.'1')THEN
            CALL BROADEN(ES,ESS,NCHAN,STWV,DISP,SIGMAVEC,.FALSE.)
          ELSE
            CALL BROADEN(ES,ESS,NCHAN,STWV,DISP,SIGMAVEC,.TRUE.)
          END IF
        END IF
        NCOLOR=NCOLOR+1
        IF(NCOLOR.GT.14) NCOLOR=2
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM))THEN
            CALL PGSCI(NCOLOR)
          END IF
          CALL PGBIN(NCHAN,X,SS,.TRUE.)
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
        END DO
C
        WRITE(*,100)'Save broadened spectrum (y/n) '
        CSAVE(1:1)=READC('n','yn')
        IF(CSAVE.EQ.'y')THEN
          WRITE(*,100)'Output file name'
          OUTFILE=OUTFILEX(30,'@',1,NCHAN,STWV,DISP,1,.FALSE.)
          WRITE(30) (SS(J),J=1,NCHAN)
          CLOSE(30)
          IF(CERR.EQ.'y')THEN
            WRITE(*,100)'Output error file name '
            CALL GUESSEF(OUTFILE,ERRFILE)
            ERRFILE=OUTFILEX(30,ERRFILE,1,NCHAN,STWV,DISP,1,.TRUE.)
            WRITE(30) (ESS(J),J=1,NCHAN)
            CLOSE(30)
          END IF
        END IF
        GOTO 50
C------------------------------------------------------------------------------
        STOP
100     FORMAT(A,$)
101     FORMAT(A)
        END
