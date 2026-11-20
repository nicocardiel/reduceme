!------------------------------------------------------------------------------
! Version 20-November-2026                                   file: calambda.f90
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
! Program: calambda
! Classification: wavelengths
! Description: determine the wavelength as a function of the channel number,
! using the wavelength calibration polynomial.
!
!Comment
!
! Calcula la longitud de onda de un determinado pixel a partir del polinomio
! de calibracion en l.d.o.
!
        PROGRAM CALAMBDA
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INTEGER READI
        INTEGER READILIM
!
        INTEGER K
        INTEGER NTERMS
        INTEGER NP
        REAL POL1,POL2
        REAL A(20)
        CHARACTER*255 POLFILE
!------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
        OUTFILEX=OUTFILEX
!
        THISPROGRAM='calambda'
        CALL WELCOME('20-November-2026')
        CALL SHOWHLP('explanation')
!
        WRITE(*,100,ADVANCE='NO')'Polynomial file name'
        POLFILE=INFILEX(20,'@',0,0,0.,0.,3,.FALSE.)
!
        NTERMS=0
10      READ(20,*,END=20)K,A(NTERMS+1)
        NTERMS=NTERMS+1
        GOTO 10
20      CLOSE(20)
!
        WRITE(*,100,ADVANCE='NO')'NCHAN'
        NCHAN=READI('@')
!
        WRITE(*,100,ADVANCE='NO')'Pol. degree: '
        WRITE(*,*) NTERMS-1
!
        NP=1
        POL1=A(NTERMS)
        DO K=NTERMS-1,1,-1
          POL1=POL1*REAL(NP)+A(K)
        END DO
        WRITE(*,'(A,I5,A)',ADVANCE='NO')'Pixel: ',NP,'  Lambda: '
        WRITE(*,*) POL1
        NP=NCHAN
        POL2=A(NTERMS)
        DO K=NTERMS-1,1,-1
          POL2=POL2*REAL(NP)+A(K)
        END DO
        WRITE(*,'(A,I5,A)',ADVANCE='NO')'Pixel: ',NP,'  Lambda: '
        WRITE(*,*) POL2
        WRITE(*,100,ADVANCE='NO')'Dispersion: '
        WRITE(*,*)(POL2-POL1)/REAL(NCHAN-1)
!
30      WRITE(*,100,ADVANCE='NO')'Pixel (0=EXIT) '
        NP=READILIM('0',0,NCHAN)
        IF(NP.EQ.0)STOP
        POL1=A(NTERMS)
        DO K=NTERMS-1,1,-1
          POL1=POL1*REAL(NP)+A(K)
        END DO
        WRITE(*,'(A,I5,A)',ADVANCE='NO')'Pixel: ',NP,'  Lambda: '
        WRITE(*,*) POL1
        GOTO 30
!
        STOP
100     FORMAT(A)
        END
