C------------------------------------------------------------------------------
C Version 7-July-1998                                         File: fftcorrel.f
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
C SUBROUTINE FFTCORREL(N,DATA1,DATA2,XCORR,FCORR)
C
C Input N,DATA1,DATA2
C Output XCORR,FCORR
C
C Compute the correlation function FCORR of two real data vectors DATA1 and 
C DATA2, using FFT. We assume that both data sets have been properly filtered
C and dimensioned. We use the discrete correlation theorem, which says that
C the discrete correlation of two real functions f1 and f2 is one member of
C the discrete Fourier transform pair:
C                       Corr(f1,f2) <==> F1 F2*
C where F1 and F2 are the discrete Fourier transforms of f1 and f2,
C respectively, and asterisk denotes complex conjugation.
C Note that DATA1 and DATA2 are not modified by this routine.
C
C INTEGER N  -> number of points (must be a power of 2)
C REAL DATA1(N) -> first data set to be correlated
C REAL DATA2(N) -> second data set to be correlated
C REAL XCORR(N) -> abcissas of the output correlation function
C REAL FCORR(N) -> output correlation function
C
Comment
C------------------------------------------------------------------------------
        SUBROUTINE FFTCORREL(N,DATA1,DATA2,XCORR,FCORR)
        IMPLICIT NONE
C
        INTEGER N
        REAL DATA1(N),DATA2(N)
        REAL XCORR(N),FCORR(N)
C local parameters
        INTEGER NMAX
        PARAMETER (NMAX=8192)
C local variables
        INTEGER I,N2
        REAL XR1(NMAX),XR2(NMAX)
        REAL XI1(NMAX),XI2(NMAX)
        REAL XR(NMAX),XI(NMAX)
C------------------------------------------------------------------------------
        IF(N.GT.NMAX)THEN
          WRITE(*,100)'N='
          WRITE(*,*)N
          WRITE(*,101)'FATAL ERROR: N is too large.'
          STOP
        END IF
C------------------------------------------------------------------------------
C calculamos la FFT de los dos conjuntos de datos (reales)
        DO I=1,N
          XR1(I)=DATA1(I)
          XR2(I)=DATA2(I)
          XI1(I)=0.
          XI2(I)=0.
        END DO
        CALL CFFT(N,XR1,XI1,1)
        CALL CFFT(N,XR2,XI2,1)
C calculamos el complejo conjugado de la segunda FFT
        DO I=1,N
          XI2(I)=-XI2(I)
        END DO
C multiplicamos FFT(DATA1) FFT(DATA2)*
        DO I=1,N
          XR(I)=XR1(I)*XR2(I)-XI1(I)*XI2(I)
          XI(I)=XR1(I)*XI2(I)+XI1(I)*XR2(I)
        END DO
C invertimos la FFT del producto anterior
        CALL CFFT(N,XR,XI,-1)
C reorganizamos la salida
        N2=N/2
        DO I=1,N2
          XCORR(I)=-REAL(N2-I+1)
          XCORR(N2+I)=REAL(I-1)
          FCORR(I)=XR(N2+I)
          FCORR(N2+I)=XR(I)
        END DO
C------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
        END
