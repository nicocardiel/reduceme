!------------------------------------------------------------------------------
! Version 20-November-2026                                  file: fft2power.f90
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
! SUBROUTINE FFT2POWER(N0,N)
!
! Input: N0
! Output: N
!
! Given an integer N0, this subroutine asks for a power of 2, such as
! N=2**K, being K integer and N.GE.N0. This value of N is employed by other
! subroutines to perform zero padding prior computing FFT.
!
!Comment
!------------------------------------------------------------------------------
! Calcula la el valor de N=2**K que sobrepasa a N0
        SUBROUTINE FFT2POWER(N0,N)
        IMPLICIT NONE
        INTEGER N0,N
!
        INCLUDE 'futils.inc'
        INTEGER READILIM
!
        INTEGER NMAX         !si se cambia, hacerlo tambien en otras subrutinas
        PARAMETER (NMAX=8192)
!
        INTEGER K,NZPAD
        CHARACTER*1 COK
        CHARACTER*50 CDUMMY
!------------------------------------------------------------------------------
5       K=0
        N=1
        DO WHILE(N.LT.N0)
          K=K+1
          N=N*2
        END DO
!
        WRITE(*,*)
        WRITE(*,101)'>>> Zero padding:'
10      NZPAD=N-N0
        IF(NZPAD.LT.0)THEN
          WRITE(*,101)'ERROR: invalid power of 2.'
          GOTO 5
        END IF
        WRITE(*,110)'Number of data points........................: ',&
     &   N0
        WRITE(*,110)'Nearest power of 2 above no. of data points..: ',&
     &   K
        WRITE(*,110)'Initial data dimension power of 2............: ',&
     &   N
        WRITE(*,110)'Zero padding will expand over (channels).....: ',&
     &   NZPAD
        WRITE(*,100,ADVANCE='NO')'Are these numbers OK (y/n) '
        IF(NZPAD.GT.0)THEN
          COK(1:1)=READC('y','yn')
        ELSE
          COK(1:1)=READC('n','yn')
        END IF
        IF(COK.EQ.'n')THEN
20        WRITE(*,100,ADVANCE='NO')'New power of 2 '
          WRITE(CDUMMY,*)K+1
          K=READILIM(CDUMMY,0,9999)
          N=2**K
          IF(N.GT.NMAX)THEN
            WRITE(*,101)'ERROR: number out of limits. Try again.'
            K=12
            GOTO 20
          END IF
          GOTO 10
        END IF
!
100     FORMAT(A)
101     FORMAT(A)
110     FORMAT(A,I6)
!
        END
