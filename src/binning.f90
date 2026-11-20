!------------------------------------------------------------------------------
! Version 20-November-2026                                    file: binning.f90
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
! Program: binning
! Classification: arithmetic & manipulations
! Description: Perform a constant binning in the spatial and wavelength 
! direction.
!
!Comment
!
! A partir de una imagen inicial genera otra imagen en la que se ha realizado
! un binning constante en la direccion espacial y en la direccion longitud de
! onda.
!
        PROGRAM BINNING
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER READILIM
!
        INTEGER I,J,IBIN,JBIN
        INTEGER XBIN,YBIN
        INTEGER NCHANBIN,NSCANBIN
        INTEGER NS1,NS2,NC1,NC2
        REAL A(NCMAX,NSMAX),S(NCMAX)
        REAL ERR(NCMAX,NSMAX)
        REAL STWVBIN,DISPBIN
        CHARACTER*1 CERR,COK
        CHARACTER*255 INFILE,ERRFILE,OUTFILE
!------------------------------------------------------------------------------
        THISPROGRAM='binning'
        CALL WELCOME('20-November-2026')
        CALL SHOWHLP('explanation')
!
        WRITE(*,100,ADVANCE='NO')'Work with error images (y/n) '
        CERR(1:1)=READC('n','yn')
!
        WRITE(*,100,ADVANCE='NO')'Input file name'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        DO I=1,NSCAN
          READ(20) (A(J,I),J=1,NCHAN)
        END DO
        CLOSE(20)
        IF(CERR.EQ.'y')THEN
          WRITE(*,100,ADVANCE='NO')'Input error file name '
          CALL GUESSEF(INFILE,ERRFILE)
          ERRFILE=INFILEX(21,ERRFILE,NSCAN,NCHAN,STWV,DISP,21,.TRUE.)!....match
          DO I=1,NSCAN
            READ(21) (ERR(J,I),J=1,NCHAN)
          END DO
          CLOSE(21)
        END IF
!------------------------------------------------------------------------------
10      WRITE(*,100,ADVANCE='NO')'X-binning '
        XBIN=READILIM('1',1,NCHAN)
        IF(MOD(NCHAN,XBIN).NE.0)THEN
          WRITE(*,101)'ERROR: invalid number: MOD(NCHAN,XBIN).NE.0'
          WRITE(*,101)'Try again.'
          GOTO 10
        END IF
12      WRITE(*,100,ADVANCE='NO')'Y-binning '
        YBIN=READILIM('1',1,NSCAN)
        IF(MOD(NSCAN,YBIN).NE.0)THEN
          WRITE(*,101)'ERROR: invalid number: MOD(NSCAN,YBIN).NE.0'
          WRITE(*,101)'Try again.'
          GOTO 12
        END IF
!
        NSCANBIN=NSCAN/YBIN
        NCHANBIN=NCHAN/XBIN
        STWVBIN=STWV+0.5*REAL(XBIN-1)*DISP
        DISPBIN=DISP*REAL(XBIN)
        WRITE(*,*)
        WRITE(*,110)'>>> NSCAN in output: ',NSCANBIN
        WRITE(*,110)'>>> NCHAN in output: ',NCHANBIN
        WRITE(*,100,ADVANCE='NO')'>>> STWV  in output:'
        WRITE(*,*) STWVBIN
        WRITE(*,100,ADVANCE='NO')'>>> DISP  in output:'
        WRITE(*,*) DISPBIN
        WRITE(*,100,ADVANCE='NO')'Are these values OK (y/n) '
        COK(1:1)=READC('y','yn')
        IF(COK.EQ.'n') GOTO 10
!------------------------------------------------------------------------------
        WRITE(*,100,ADVANCE='NO')'Output file name......'
        OUTFILE=OUTFILEX(30,'@',NSCANBIN,NCHANBIN,STWVBIN,DISPBIN,&
     &   1,.FALSE.)
!------------------------------------------------------------------------------
! Salvamos imagen con binning
        DO IBIN=1,NSCANBIN
          NS1=(IBIN-1)*YBIN+1
          NS2=NS1+YBIN-1
          DO JBIN=1,NCHANBIN
            NC1=(JBIN-1)*XBIN+1
            NC2=NC1+XBIN-1
            S(JBIN)=0.
            DO I=NS1,NS2
              DO J=NC1,NC2
                S(JBIN)=S(JBIN)+A(J,I)
              END DO
            END DO
          END DO
          WRITE(30) (S(JBIN),JBIN=1,NCHANBIN)
        END DO
        CLOSE(30)
!------------------------------------------------------------------------------
! Salvamos errores con binning
        IF(CERR.EQ.'y')THEN
          WRITE(*,100,ADVANCE='NO')'Output error file name '
          CALL GUESSEF(OUTFILE,ERRFILE)
          OUTFILE=OUTFILEX(31,ERRFILE,NSCANBIN,NCHANBIN,STWVBIN,&
     &     DISPBIN,1,.TRUE.)
          DO IBIN=1,NSCANBIN
            NS1=(IBIN-1)*YBIN+1
            NS2=NS1+YBIN-1
            DO JBIN=1,NCHANBIN
              NC1=(JBIN-1)*XBIN+1
              NC2=NC1+XBIN-1
              S(JBIN)=0.
              DO I=NS1,NS2
                DO J=NC1,NC2
                  S(JBIN)=S(JBIN)+ERR(J,I)*ERR(J,I)
                END DO
              END DO
            END DO
            DO JBIN=1,NCHANBIN
              S(JBIN)=SQRT(S(JBIN))
            END DO
            WRITE(31) (S(JBIN),JBIN=1,NCHANBIN)
          END DO
          CLOSE(31)
        END IF
!------------------------------------------------------------------------------
        STOP
100     FORMAT(A)
101     FORMAT(A)
110     FORMAT(A,I6)
        END
