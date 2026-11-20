!------------------------------------------------------------------------------
! Version 20-November-2026                                       File: cfft.f90
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
! SUBROUTINE CFFT(N,XR,XI,IMODE)
!
! Input N,XR,XI,IMODE
! Output XR,YR
!
! If IMODE=1, this subroutine computes the FFT of the input complex vector 
! XR + i XI, where XR and XI are both real variables. As output XR and XI are 
! the real and imaginary part of the transform. If IMODE=-1 the subroutine
! evaluates the inverse FFT. See E. O. Brigham, The Fast Fourier Transform, 
! pag.160.
!
! INTEGER N  -> number of points (must be a power of 2)
! REAL XR(N) -> real part of the input data
! REAL XI(N) -> imaginary part of the input data
! INTEGER IMODE -> +1 (direct FFT) or -1 (inverse FFT)
!
!Comment
!------------------------------------------------------------------------------
        SUBROUTINE CFFT(N,XR,XI,IMODE)
        IMPLICIT NONE
!
        INTEGER N
        REAL XR(N),XI(N)
        INTEGER IMODE
! local variables
        INTEGER NU,N0
        INTEGER N2,NU1,K,L,I,P,K1,K1N2
        REAL XR0,XI0,FN
        DOUBLE PRECISION DARG,CC,SS,TREAL,TIMAG
!------------------------------------------------------------------------------
! NU: power of 2 such as 2**NU=N
        NU=0
        N0=N
        DO WHILE(N0.GT.2)
          IF(MOD(N0,2).NE.0)THEN
            WRITE(*,*)
            WRITE(*,100,ADVANCE='NO')'>>> N='
            WRITE(*,*) N
            WRITE(*,101)'FATAL ERROR: N is not a power of 2 in CFFT'
            STOP
          END IF
          N0=N0/2
          NU=NU+1
        END DO
        NU=NU+1
!------------------------------------------------------------------------------
! To compute the inverse FFT remember that, given a function f and its Fourier 
! transform F: F = FFT(f) and  f = 1/N [FFT(f*)]*, where the asterisks mean 
! complex conjugate.
        IF(IMODE.EQ.-1)THEN
          DO I=1,N
            XI(I)=-XI(I)
          END DO
        END IF
!------------------------------------------------------------------------------
! initialization
        N2=N/2
        NU1=NU-1
        K=0
!------------------------------------------------------------------------------
! computing the FFT
        DO L=1,NU
          DO WHILE(K.LT.N)
            DO I=1,N2
              P=IBITRFFT(K/2**NU1,NU)
              DARG=6.28318530717959D0*DBLE(P)/DBLE(N)
              CC=DCOS(DARG)
              SS=DSIN(DARG)
              K1=K+1
              K1N2=K1+N2
              TREAL=DBLE(XR(K1N2))*CC+DBLE(XI(K1N2))*SS
              TIMAG=DBLE(XI(K1N2))*CC-DBLE(XR(K1N2))*SS
              XR(K1N2)=XR(K1)-REAL(TREAL)
              XI(K1N2)=XI(K1)-REAL(TIMAG)
              XR(K1)=XR(K1)+REAL(TREAL)
              XI(K1)=XI(K1)+REAL(TIMAG)
              K=K+1
            END DO
            K=K+N2
          END DO
          K=0
          NU1=NU1-1
          N2=N2/2
        END DO
!------------------------------------------------------------------------------
! unscrambling the FFT
        DO K=1,N
          I=IBITRFFT(K-1,NU)+1
          IF(I.GT.K)THEN
            XR0=XR(K)
            XI0=XI(K)
            XR(K)=XR(I)
            XI(K)=XI(I)
            XR(I)=XR0
            XI(I)=XI0
          END IF
        END DO
!------------------------------------------------------------------------------
! inverse FFT
        IF(IMODE.EQ.-1)THEN
          DO I=1,N
            XI(I)=-XI(I)
          END DO
          FN=1./REAL(N)
          DO I=1,N
            XR(I)=XR(I)*FN
            XI(I)=XI(I)*FN
          END DO
        END IF
!------------------------------------------------------------------------------
100     FORMAT(A)
101     FORMAT(A)
        END
!
!******************************************************************************
! Bit-reversing function
        INTEGER FUNCTION IBITRFFT(J,NU)
        IMPLICIT NONE
!
        INTEGER J,NU
! local variables
        INTEGER I,J1,J2
!------------------------------------------------------------------------------
        J1=J
        IBITRFFT=0
        DO I=1,NU
          J2=J1/2
          IBITRFFT=IBITRFFT*2+(J1-2*J2)
          J1=J2
        END DO
!
        END
