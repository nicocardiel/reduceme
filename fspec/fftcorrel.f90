!------------------------------------------------------------------------------
! Version 20-November-2026                                  File: fftcorrel.f90
!------------------------------------------------------------------------------
! Copyright N. Cardiel & J. Gorgas, Departamento de Astrofisica
! Universidad Complutense de Madrid, 28040-Madrid, Spain
! E-mail: ncl@astrax.fis.ucm.es or fjg@astrax.fis.ucm.es
!------------------------------------------------------------------------------
! This routine is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the Free
! Software Foundation; either version 2 of the License, or (at your option) any
! later version. See the file gnu-public-license.txt for details.
!------------------------------------------------------------------------------
!Comment
!
! SUBROUTINE FFTCORREL(N,DATA1,DATA2,XCORR,FCORR)
!
! Input N,DATA1,DATA2
! Output XCORR,FCORR
!
! Compute the correlation function FCORR of two real data vectors DATA1 and 
! DATA2, using FFT. We assume that both data sets have been properly filtered
! and dimensioned. We use the discrete correlation theorem, which says that
! the discrete correlation of two real functions f1 and f2 is one member of
! the discrete Fourier transform pair:
!                       Corr(f1,f2) <==> F1 F2*
! where F1 and F2 are the discrete Fourier transforms of f1 and f2,
! respectively, and asterisk denotes complex conjugation.
! Note that DATA1 and DATA2 are not modified by this routine.
!
! INTEGER N  -> number of points (must be a power of 2)
! REAL DATA1(N) -> first data set to be correlated
! REAL DATA2(N) -> second data set to be correlated
! REAL XCORR(N) -> abcissas of the output correlation function
! REAL FCORR(N) -> output correlation function
!
!Comment
!------------------------------------------------------------------------------
        SUBROUTINE FFTCORREL(N,DATA1,DATA2,XCORR,FCORR)
        IMPLICIT NONE
!
        INTEGER N
        REAL DATA1(N),DATA2(N)
        REAL XCORR(N),FCORR(N)
! local parameters
        INTEGER NMAX
        PARAMETER (NMAX=8192)
! local variables
        INTEGER I,N2
        REAL XR1(NMAX),XR2(NMAX)
        REAL XI1(NMAX),XI2(NMAX)
        REAL XR(NMAX),XI(NMAX)
!------------------------------------------------------------------------------
        IF(N.GT.NMAX)THEN
          WRITE(*,100,ADVANCE='NO')'N='
          WRITE(*,*)N
          WRITE(*,101)'FATAL ERROR: N is too large.'
          STOP
        END IF
!------------------------------------------------------------------------------
! calculamos la FFT de los dos conjuntos de datos (reales)
        DO I=1,N
          XR1(I)=DATA1(I)
          XR2(I)=DATA2(I)
          XI1(I)=0.
          XI2(I)=0.
        END DO
        CALL CFFT(N,XR1,XI1,1)
        CALL CFFT(N,XR2,XI2,1)
! calculamos el complejo conjugado de la segunda FFT
        DO I=1,N
          XI2(I)=-XI2(I)
        END DO
! multiplicamos FFT(DATA1) FFT(DATA2)*
        DO I=1,N
          XR(I)=XR1(I)*XR2(I)-XI1(I)*XI2(I)
          XI(I)=XR1(I)*XI2(I)+XI1(I)*XR2(I)
        END DO
! invertimos la FFT del producto anterior
        CALL CFFT(N,XR,XI,-1)
! reorganizamos la salida
        N2=N/2
        DO I=1,N2
          XCORR(I)=-REAL(N2-I+1)
          XCORR(N2+I)=REAL(I-1)
          FCORR(I)=XR(N2+I)
          FCORR(N2+I)=XR(I)
        END DO
!------------------------------------------------------------------------------
100     FORMAT(A)
101     FORMAT(A)
        END
