C
C******************************************************************************
C Ajusta un polinomio a N puntos X(N),Y(N), eliminando puntos que se alejen
C mas de TSGIMA veces. La rutina devuelve en LFIT2 que puntos han sido
C finalmente utilizados en el ajuste.
        SUBROUTINE POLFITSIG(N,X,Y,TSIGMA,NDEG,A,LFIT2)
        IMPLICIT NONE
        INTEGER N
        REAL X(N),Y(N)
        INTEGER NDEG
        REAL TSIGMA,A(NDEG+1)
        LOGICAL LFIT2(N)
C
        INTEGER NMAX
        PARAMETER (NMAX=5000)
C
        REAL FPOLY
C
        INTEGER I
        INTEGER NFIT,K
        REAL XFIT(NMAX),YFIT(NMAX)
        REAL CHISQR,RMS,YDUM
        LOGICAL LFIT1(NMAX)
        LOGICAL LEXIT
C------------------------------------------------------------------------------
        IF(N.GT.NMAX)THEN
          WRITE(*,100) 'N, NMAX: '
          WRITE(*,*) N,NMAX
          WRITE(*,100) 'FATAL ERROR in subroutine POLFITSIG:'
          WRITE(*,101) 'N.GT.NMAX'
          STOP
        END IF
C primero ajustamos con todos los puntos
        DO I=1,N
          XFIT(I)=X(I)
          YFIT(I)=Y(I)
          LFIT1(I)=.TRUE.
        END DO
        NFIT=N
C ajustamos el polinomio
10        CALL POLFIT(XFIT,YFIT,YFIT,NFIT,NDEG+1,0,A,CHISQR)
        IF(N.EQ.1)THEN
          LFIT2(1)=.TRUE.
          RETURN
        END IF
C calculamos varianza residual (con los puntos ajustados) y los puntos que se 
C desvian (de toda la muestra)
        RMS=0.
        DO I=1,NFIT
          YDUM=FPOLY(NDEG,A,XFIT(I))
          RMS=RMS+(YFIT(I)-YDUM)*(YFIT(I)-YDUM)
        END DO
        RMS=SQRT(RMS/REAL(NFIT-1))
        DO I=1,N
          YDUM=FPOLY(NDEG,A,X(I))
          LFIT2(I)=(ABS(Y(I)-YDUM).LT.TSIGMA*RMS) 
        END DO
C comprobamos si ha cambiado algun punto
        LEXIT=.TRUE.
        DO I=1,N
          IF(LFIT1(I).NEQV.LFIT2(I)) LEXIT=.FALSE.
        END DO
C
        IF(.NOT.LEXIT)THEN
          DO I=1,N
            LFIT1(I)=LFIT2(I)
          END DO
          K=0
          DO I=1,N
            IF(LFIT1(I))THEN
              K=K+1
              XFIT(K)=X(I)
              YFIT(K)=Y(I)
            END IF
          END DO
          NFIT=K
          GOTO 10
        END IF
C
100     FORMAT(A,$)
101     FORMAT(A)
        END
C
C******************************************************************************
C
C Calcula el polinomio en el valor de X
        REAL FUNCTION FPOLY(NDEG,COEFF,X)
        IMPLICIT NONE
        INTEGER NDEG
        REAL COEFF(NDEG+1)
        REAL X
C
        INTEGER K
C------------------------------------------------------------------------------
        FPOLY=COEFF(NDEG+1)
        IF(NDEG.GT.0)THEN
          DO K=NDEG,1,-1
            FPOLY=FPOLY*X+COEFF(K)
          END DO
        END IF
C
        END
