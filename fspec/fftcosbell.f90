!------------------------------------------------------------------------------
! Version 20-November-2026                                   file: fftcosbell.f
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
! SUBROUTINE FFTCOSBELL(N,COSBELL,FL)
!
! Input: N, FL
! Output: COSBELL
!
! Compute a cosine bell expanding N pixels, being FL the fraction of pixels (at
! the beginning and at the end of the cosine bell) employed to perform the
! transition from zero to one. See Brault & White, A&A, 13, 169.
!
! INTEGER N -> dimension of COSBELL (pixels)
! REAL COSBELL(N) -> cosine bell
! REAL FL -> fraction of pixels at the borders of the cosine bell
!
!Comment
!------------------------------------------------------------------------------
        SUBROUTINE FFTCOSBELL(N,COSBELL,FL)
        IMPLICIT NONE
        REAL READF
!
        INTEGER N
        REAL COSBELL(N)
        REAL FL
!
        REAL PI
        PARAMETER (PI=3.141592654)
!
        INTEGER NL,J
        CHARACTER*50 CDUMMY
!------------------------------------------------------------------------------
! parametro de la campana de coseno
        WRITE(*,*)
        WRITE(*,101)'>>> Cosine bell:'
10      WRITE(*,100,ADVANCE='NO')&
     & 'Length fraction over which the data will '//&
     &   'be masked '
        WRITE(CDUMMY,*)FL
        FL=READF(CDUMMY)
        IF((FL.LT.0.0).OR.(FL.GT.0.5))THEN
          WRITE(*,101)'ERROR: invalid entry. Try again.'
          GOTO 10
        END IF
!
        NL=FL*REAL(N)
        WRITE(*,100,ADVANCE='NO')&
     &     'Cosine bell will be applied over the first and '//&
     &   'last (channels): '
        WRITE(*,*) NL
!------------------------------------------------------------------------------
! calculamos la campana de coseno
        DO J=1,NL
          COSBELL(J)=0.5*(1.-COS(PI*REAL(J)/REAL(NL)))
        END DO
!
        DO J=NL+1,N-NL
          COSBELL(J)=1.0
        END DO
!
        DO J=N-NL+1,N
          COSBELL(J)=0.5*(1.-COS(PI*REAL(N-J)/REAL(NL)))
        END DO
!------------------------------------------------------------------------------
100      FORMAT(A)
101      FORMAT(A)
        END
