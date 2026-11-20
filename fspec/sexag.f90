C **********************************************************************
C                                                       SUBROUTINE SEXAG
C                                                       ****************
      SUBROUTINE SEXAG(A,B,C,D)
C
C Paso de DD.dddd a DD.MMSS
C
      IMPLICIT NONE
C---> argumentos ficticios: ENTRADA
      DOUBLE PRECISION A
C---> argumentos ficticios: SALIDA
      INTEGER B,C
      REAL D
C---> variables locales
      DOUBLE PRECISION AA,BR
C
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
C hay que tener cuidado con la siguiente asignacion de signos
C (necesaria por si B=0 o/y C=0)
      IF(A.LT.0.)THEN
        B=-B
        C=-C
        D=-D
      END IF
      END
