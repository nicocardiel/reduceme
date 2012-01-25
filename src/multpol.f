C------------------------------------------------------------------------------
C Version 02-September-2008                                     file: multpol.f
C @ ncl & fjg
C------------------------------------------------------------------------------
Comment
C
C Program: multpol
C Classification: wavelengths
C Description: Determines a new wavelength calibration polynomial from an 
C initial polynomial and a second polynomial (second-order correction).
C
Comment
C------------------------------------------------------------------------------
C
C Programa para calcular un nuevo polinomio de calibracion en longitud de onda
C a partir de un polinomio de calibracion inicial y de un segundo polinomio que
C contenga la correccion a realizar. Notar que el primer polinomio es de la
C forma:
C pol = sum_{i=0}^{N} a_i * (x**)^i, donde x es numero de canal y pol viene
C                                    dado en longitud de onda
C mientras que el segundo polinomio es de la forma:
C corr = sum_{j=0}^{M} b_j * (lambda)^j, donde lambda es longitud de onda
C------------------------------------------------------------------------------
C g77 -Wall -o multpol multpol.f 
        PROGRAM MULTPOL
        IMPLICIT NONE
C
        INTEGER NMAX,MMAX
        PARAMETER (NMAX=10,MMAX=10)
C
        INTEGER I,J,K
        INTEGER JJ
        INTEGER N,M
        INTEGER IDUM,JDUM
        DOUBLE PRECISION A(0:NMAX),B(0:MMAX)
        DOUBLE PRECISION AA(0:NMAX*MMAX),AAA(0:NMAX*MMAX)
        DOUBLE PRECISION T(0:NMAX*MMAX),C(0:NMAX*MMAX)
        CHARACTER*255 POLFILE
        CHARACTER*255 OUTFILE
        LOGICAL LOGFILE
C------------------------------------------------------------------------------
        LOGFILE=.FALSE.
        DO WHILE(.NOT.LOGFILE)
          WRITE(*,101) 'Enter file name with initial wavelength '//
     +     'calibration polynomial:'
          READ(*,100) POLFILE
          INQUIRE(FILE=POLFILE,EXIST=LOGFILE)
          IF(.NOT.LOGFILE)THEN
            WRITE(*,101) 'ERROR: this file does not exist. Try again.'
          END IF
        END DO
        OPEN(10,FILE=POLFILE,STATUS='OLD',FORM='FORMATTED')
        I=0
10      READ(10,*,END=12) IDUM,A(I)
        I=I+1
        IF(I.GT.NMAX) STOP 'FATAL ERROR: redim NMAX'
        GOTO 10
12      CLOSE(10)
        N=I-1
        WRITE(*,100) '>>> Polynomial degree: '
        WRITE(*,*) N
        IF(N.LT.1) STOP 'FATAL ERROR: N.LT.1'
        DO I=0,N
          WRITE(*,'(A2,I2.2,A2,$)') 'a(',I,')='
          WRITE(*,*) A(I)
        END DO
        WRITE(*,*)
C..............................................................................
        LOGFILE=.FALSE.
        DO WHILE(.NOT.LOGFILE)
          WRITE(*,101) 'Enter file name with wavelength '//
     +     'correction polynomial:'
          READ(*,100) POLFILE
          INQUIRE(FILE=POLFILE,EXIST=LOGFILE)
          IF(.NOT.LOGFILE)THEN
            WRITE(*,101) 'ERROR: this file does not exist. Try again.'
          END IF
        END DO
        OPEN(10,FILE=POLFILE,STATUS='OLD',FORM='FORMATTED')
        J=0
20      READ(10,*,END=22) JDUM,B(J)
        J=J+1
        IF(J.GT.MMAX) STOP 'FATAL ERROR: redim MMAX'
        GOTO 20
22      CLOSE(10)
        M=J-1
        WRITE(*,100) '>>> Polynomial degree: '
        WRITE(*,*) M
        DO J=0,M
          WRITE(*,'(A2,I2.2,A2,$)') 'b(',J,')='
          WRITE(*,*) B(J)
        END DO
        WRITE(*,*)
        WRITE(*,100) 'Press <CR> to continue...'
        READ(*,*)
C------------------------------------------------------------------------------
        DO K=0,M*N
          T(K)=0.D0
        END DO
C
        DO J=0,M
          IF(J.EQ.0)THEN
            T(J)=T(J)+B(J)
          ELSE
            IF(J.EQ.1)THEN
              DO I=0,N
                AA(I)=A(I)
              END DO
            ELSE
              DO I=0,N
                AA(I)=A(I)
              END DO
              DO JJ=1,J-1
                CALL MULT2POL(N,A,N*JJ,AA,AAA)
                DO I=0,N*(JJ+1)
                  AA(I)=AAA(I)
                END DO
              END DO
            END IF
            DO K=0,N*J
              T(K)=T(K)+B(J)*AA(K)
            END DO
          END IF
        END DO
C
        IF(M.EQ.0)THEN
          C(0)=A(0)+T(0)
          DO K=1,N
            C(K)=A(K)
          END DO
        ELSE
          DO K=0,N*M
            IF(K.LE.N)THEN
              C(K)=A(K)+T(K)
            ELSE
              C(K)=T(K)
            END IF
          END DO
        END IF
C
        DO K=0,MAX(N,N*M)
          WRITE(*,*) K,C(K)
        END DO
C
        WRITE(*,100) 'Output file name? '
        READ(*,101) OUTFILE
        OPEN(30,FILE=OUTFILE,STATUS='UNKNOWN',FORM='FORMATTED')
        DO K=0,MAX(N,N*M)
          WRITE(30,*) K,C(K)
        END DO
        CLOSE(30)
C------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
        END
C
C******************************************************************************
C Multiplica dos polinomios de la forma
C p1 = sum_{i=0}^N a(i) * (x**i)
C p2 = sum_{j=0}^M b(j) * (x**j)
C El resultado será otro polinomio de grado N+M de la forma
C p  = sum_{k=0}^{N+M} c(k) (x**k)
C
        SUBROUTINE MULT2POL(N,A,M,B,C)
        IMPLICIT NONE
        INTEGER N
        DOUBLE PRECISION A(0:N)
        INTEGER M
        DOUBLE PRECISION B(0:M)
        DOUBLE PRECISION C(0:N+M)
C
        INTEGER I,J,K
C------------------------------------------------------------------------------
C inicializamos a cero el polinomio de salida
        DO K=0,N+M
          C(K)=0.D0
        END DO
C calculamos el producto
        DO I=0,N
          DO J=0,M
            C(I+J)=C(I+J)+A(I)*B(J)
          END DO
        END DO
C
        END
