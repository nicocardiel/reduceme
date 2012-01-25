C------------------------------------------------------------------------------
C Version 3-July-1998                                        file: fftkfilter.f
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
C SUBROUTINE FFTKFILTER(N,KFILTER,K1,K2,K3,K4)
C
C Input: N,K1,K2,K3,K4
C Output: KFILTER
C
C Given K1, K2, K3 and K4, this subroutine computes the filter KFILTER, such
C as:
C
C KFILTER(I)=0., I=1,...,K1
C KFILTER(I) is a straight line from 0. to 1., I=K1,...,K2
C KFILTER(I)=1., I=K2,...,K3
C KFILTER(I) is a straight line from 1. to 0., I=K3,...,K4
C KFILTER(I)=0., I=K4,...,N/2+1
C
C KFILTER(I)=0., I=N,...,N-K1+2
C KFILTER(I) is a straight line from 0. to 1., I=N-K1+2,...,N-K2+2
C KFILTER(I)=1., I=N-K2+2,...,N-K3+2
C KFILTER(I) is a straight line from 1. to 0., I=N-K3+2,...,N-K4+2
C KFILTER(I)=0., I=N-K4+2,...,N/2+1
C
C N must be a power of 2, and 1.le.K1.le.K2.le.K3.le.K4.le.(N/2-1).
C This filter is employed by other routines to perform a frequency filter
C in the frequency domain of the FFT.
C
Comment
C------------------------------------------------------------------------------
        SUBROUTINE FFTKFILTER(N,KFILTER,K1,K2,K3,K4)
        IMPLICIT NONE
        INCLUDE 'futils.inc'
        INTEGER READI
C
        INTEGER N
        REAL KFILTER(N)
        INTEGER K1,K2,K3,K4
C
        INTEGER J,L
        INTEGER N0,K
        CHARACTER*1 CK
        CHARACTER*120 CDUMMY
C------------------------------------------------------------------------------
C chequeamos que N es una potencia de 2
        K=0
        N0=1
        DO WHILE(N0.LT.N)
          K=K+1
          N0=2*N0
        END DO
        IF(N0.NE.N)THEN
          WRITE(*,101)'FATAL ERROR in subroutine FFTKFILTER: '
          WRITE(*,100)'N is not a power of 2!'
          STOP
        END IF
C------------------------------------------------------------------------------
C Definimos los parametros del filtro
        WRITE(*,*)
        WRITE(*,101)'>>> Filtering:'
        WRITE(CDUMMY,'(4(A4,I8))')' K1=',K1,',K2=',K2,',K3=',K3,
     +   ',K4=',K4
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,101)CDUMMY(1:L)
C
        WRITE(*,100)'Ki must be in the range 1,'
        WRITE(CDUMMY,*)N/2+1
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,101)CDUMMY(1:L)
C
10      IF((K1.LT.1).OR.(K4.GT.N/2+1).OR.(K1.GT.K2).OR.(K2.GT.K3).OR.
     +   (K3.GT.K4))THEN
          WRITE(*,101)'ERROR: initial K1, K2, K3 or K4 are not correct.'
          CK='n'
        ELSE
          WRITE(*,100)'Are these numbers OK (y/n) '
          CK(1:1)=READC('y','yn')
        END IF
C
        IF(CK.EQ.'n')THEN
          WRITE(*,100)'K1 '
          WRITE(CDUMMY,*) K1
          K1=READI(CDUMMY)
          WRITE(*,100)'K2 '
          WRITE(CDUMMY,*) K2
          K2=READI(CDUMMY)
          WRITE(*,100)'K3 '
          WRITE(CDUMMY,*) K3
          K3=READI(CDUMMY)
          WRITE(*,100)'K4 '
          WRITE(CDUMMY,*) K4
          K4=READI(CDUMMY)
          GOTO 10
        END IF
C------------------------------------------------------------------------------
C Definimos el filtro, de tal forma que hasta K1 es cero, de K1 a K2
C aumenta linealmente de cero a 1, entre K2 y K3 es uno, de K3 a K4 cae
C linealmente de 1 a cero y finalmente, de K4 en adelante (hasta N/2+1) es 
C cero. Si K1=K2, KFILTER(K1)=1.0. Analogamente, si K3=K4, KFILTER(K4)=1.0.
        IF(K1.GT.1)THEN
          DO J=1,K1-1
            KFILTER(J)=0.
          END DO
        END IF
        IF(K2.GT.K1)THEN
          DO J=K1,K2
            KFILTER(J)=REAL(J-K1)/REAL(K2-K1)
          END DO
        END IF
        DO J=K2,K3
          KFILTER(J)=1.0
        END DO
        IF(K4.GT.K3)THEN
          DO J=K3,K4
            KFILTER(J)=REAL(K4-J)/REAL(K4-K3)
          END DO
        END IF
        IF(K4.LT.N/2+1)THEN
          DO J=K4+1,N/2+1
            KFILTER(J)=0.
          END DO
        END IF
C definimos ahora el filtro para las frecuencias negativas, teniendo en
C cuenta que en el espacio transformado, las frecuencias son simetricas
C respecto al punto N/2+1
        DO J=N/2+2,N
          KFILTER(J)=KFILTER(N-J+2)
        END DO
C------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
        END
