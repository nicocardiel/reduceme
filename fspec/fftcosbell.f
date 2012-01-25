C------------------------------------------------------------------------------
C Version 24-June-1998                                       file: fftcosbell.f
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
C SUBROUTINE FFTCOSBELL(N,COSBELL,FL)
C
C Input: N, FL
C Output: COSBELL
C
C Compute a cosine bell expanding N pixels, being FL the fraction of pixels (at
C the beginning and at the end of the cosine bell) employed to perform the
C transition from zero to one. See Brault & White, A&A, 13, 169.
C
C INTEGER N -> dimension of COSBELL (pixels)
C REAL COSBELL(N) -> cosine bell
C REAL FL -> fraction of pixels at the borders of the cosine bell
C
Comment
C------------------------------------------------------------------------------
        SUBROUTINE FFTCOSBELL(N,COSBELL,FL)
        IMPLICIT NONE
        REAL READF
C
        INTEGER N
        REAL COSBELL(N)
        REAL FL
C
        REAL PI
        PARAMETER (PI=3.141592654)
C
        INTEGER NL,J
        CHARACTER*50 CDUMMY
C------------------------------------------------------------------------------
C parametro de la campana de coseno
        WRITE(*,*)
        WRITE(*,101)'>>> Cosine bell:'
10      WRITE(*,100)'Length fraction over which the data will '//
     +   'be masked '
        WRITE(CDUMMY,*)FL
        FL=READF(CDUMMY)
        IF((FL.LT.0.0).OR.(FL.GT.0.5))THEN
          WRITE(*,101)'ERROR: invalid entry. Try again.'
          GOTO 10
        END IF
C
        NL=FL*REAL(N)
        WRITE(*,100)'Cosine bell will be applied over the first and '//
     +   'last (channels): '
        WRITE(*,*) NL
C------------------------------------------------------------------------------
C calculamos la campana de coseno
        DO J=1,NL
          COSBELL(J)=0.5*(1.-COS(PI*REAL(J)/REAL(NL)))
        END DO
C
        DO J=NL+1,N-NL
          COSBELL(J)=1.0
        END DO
C
        DO J=N-NL+1,N
          COSBELL(J)=0.5*(1.-COS(PI*REAL(N-J)/REAL(NL)))
        END DO
C------------------------------------------------------------------------------
100      FORMAT(A,$)
101      FORMAT(A)
        END
