! **********************************************************************
!                                                       SUBROUTINE SEXAG
!                                                       ****************
      SUBROUTINE SEXAG(A,B,C,D)
!
! Paso de DD.dddd a DD.MMSS
!
      IMPLICIT NONE
!---> argumentos ficticios: ENTRADA
      DOUBLE PRECISION A
!---> argumentos ficticios: SALIDA
      INTEGER B,C
      REAL D
!---> variables locales
      DOUBLE PRECISION AA,BR
!
      AA=DABS(A)
      B=INT(AA)
      BR=(AA-DBLE(B))*60.D0
      C=INT(BR)
      D=REAL((BR-DINT(BR))*60.D0)
      IF(D.GE.60.)THEN
        C=C+1
        D=0.0
        IF(C.GE.60)THEN
          B=B+1
          C=0
        END IF
      END IF
! hay que tener cuidado con la siguiente asignacion de signos
! (necesaria por si B=0 o/y C=0)
      IF(A.LT.0.)THEN
        B=-B
        C=-C
        D=-D
      END IF
      END
