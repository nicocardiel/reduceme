C Version 31-January-2001
C******************************************************************************
C Compute k_lambda as a function of lambda (and also as a function of Rv for
C those curves which depend on it).
C
        SUBROUTINE SCOMPEXT(COPC,LAMBDA,RV,K_LAMBDA,LOK)
        IMPLICIT NONE
        REAL LAMBDA,RV,K_LAMBDA
        CHARACTER*1 COPC
        LOGICAL LOK
C------------------------------------------------------------------------------
C COPC:
C (1) Galaxy: Savage & Mathis (1979)
C (2) Galaxy: Seaton (1979)
C (3) Galaxy: Cardelli, Clayton and Mathis (1989) + O`Donnell (1994)
C (4) Galaxy: Fitzpatrick (1999)
C (a) LMC: Howarth (1983)
C (b) LMC (30 Doradus): Fitzpatrick (1985)
C (m) SMC: Prevot et al. (1984) + Bouchet et al. (1985)
C (p) Starburst: Calzetti (1997)
C (q) Starburst: Calzetti et al. (2000)
C
C LAMBDA: wavelength, in Angstroms
C RV: = A(V)/E(B-V), which is typically 3.1
C K_LAMBDA: output of the subroutine
C LOK: .TRUE. if the wavelength range is OK, .FALSE. if it falls outside
C the valid range of the curve
C
C The output of this subroutine can be employed to correct from
C galactic/internal extinction, using:
C Flux_corrected = 10**(0.4 * E(B-V) * K_lambda) * Flux_observed
C------------------------------------------------------------------------------
        INTEGER NPMAX
        PARAMETER (NPMAX=100)
C
        INTEGER IFLAG,N1,N2
        INTEGER NPTOS
        INTEGER I0,IMODE
        REAL XX,YY,FA,FB
        REAL X(NPMAX),Y(NPMAX)
        REAL LININTERP__
        REAL FITZPATRICK99
        REAL SPLS(NPMAX),SPLA(NPMAX),SPLB(NPMAX),SPLC(NPMAX)
C------------------------------------------------------------------------------
        LOK=.FALSE. !salvo que se demuestre lo contrario
        XX=10000./LAMBDA !1/LAMBDA, with LAMBDA in microns
        IF(COPC.EQ.'1')THEN
C..............................................................................
C Galaxy: Savage and Mathis (1979) -> Linear interpolation in Table 2.
          IF((LAMBDA.GE.1000.0).AND.(LAMBDA.LE.34000.))THEN
            X(01)= 1000.000
            Y(01)=14.40
            X(02)= 1050.000
            Y(02)=12.90
            X(03)= 1110.000
            Y(03)=11.55
            X(04)= 1180.000
            Y(04)=10.55
            X(05)= 1250.000
            Y(05)= 9.65
            X(06)= 1390.000
            Y(06)= 8.49
            X(07)= 1490.000
            Y(07)= 8.15
            X(08)= 1600.000
            Y(08)= 8.12
            X(09)= 1700.000
            Y(09)= 7.87
            X(10)= 1800.000
            Y(10)= 7.75
            X(11)= 1900.000
            Y(11)= 8.00
            X(12)= 2000.000
            Y(12)= 8.62
            X(13)= 2100.000
            Y(13)= 9.33
            X(14)= 2190.000
            Y(14)= 9.67
            X(15)= 2300.000
            Y(15)= 8.87
            X(16)= 2400.000
            Y(16)= 8.00
            X(17)= 2500.000
            Y(17)= 7.29
            X(18)= 2740.000
            Y(18)= 6.20
            X(19)= 3440.000
            Y(19)= 4.90
            X(20)= 4000.000
            Y(20)= 4.40
            X(21)= 4400.000
            Y(21)= 4.10
            X(22)= 5500.000
            Y(22)= 3.10
            X(23)= 7000.000
            Y(23)= 2.32
            X(24)= 9000.000
            Y(24)= 1.50
            X(25)=12500.000
            Y(25)= 0.87
            X(26)=22000.000
            Y(26)= 0.38
            X(27)=34000.000
            Y(27)= 0.16
            NPTOS=27
            K_LAMBDA=LININTERP__(NPTOS,X,Y,LAMBDA,IFLAG,N1,N2)
            LOK=.TRUE.
          END IF
C..............................................................................
        ELSEIF(COPC.EQ.'2')THEN
C Galaxy: Seaton (1979)
C Since Seaton(1979) data are normalized to RV=3.20, it is important to
C re-normalize everything to the adopted RV value -> Expressions given in 
C Table 2 and linear interpolation in Table 3.
          IF(XX.LT.1.0)THEN
          ELSEIF(XX.LT.2.7)THEN
            X(01)=1.0
            Y(01)=1.36
            X(02)=1.1
            Y(02)=1.44
            X(03)=1.2
            Y(03)=1.84
            X(04)=1.3
            Y(04)=2.04
            X(05)=1.4
            Y(05)=2.24
            X(06)=1.5
            Y(06)=2.44
            X(07)=1.6
            Y(07)=2.66
            X(08)=1.7
            Y(08)=2.88
            X(09)=1.8
            Y(09)=3.14
            X(10)=1.9
            Y(10)=3.36
            X(11)=2.0
            Y(11)=3.56
            X(12)=2.1
            Y(12)=3.77
            X(13)=2.2
            Y(13)=3.96
            X(14)=2.3
            Y(14)=4.15
            X(15)=2.4
            Y(15)=4.26
            X(16)=2.5
            Y(16)=4.40
            X(17)=2.6
            Y(17)=4.52
            X(18)=2.7
            Y(18)=4.64
            NPTOS=18
            K_LAMBDA=LININTERP__(NPTOS,X,Y,XX,IFLAG,N1,N2)
            K_LAMBDA=K_LAMBDA+(RV-3.2)
            LOK=.TRUE.
          ELSEIF(XX.LT.3.65)THEN
            K_LAMBDA=1.56+1.048*XX+1.01/((XX-4.60)*(XX-4.60)+0.280)
            K_LAMBDA=K_LAMBDA+(RV-3.2)
            LOK=.TRUE.
          ELSEIF(XX.LT.7.14)THEN
            K_LAMBDA=2.29+0.848*XX+1.01/((XX-4.60)*(XX-4.60)+0.280)
            K_LAMBDA=K_LAMBDA+(RV-3.2)
            LOK=.TRUE.
          ELSEIF(XX.LE.10.0)THEN
            K_LAMBDA=16.17-3.20*XX+0.2975*XX*XX
            K_LAMBDA=K_LAMBDA+(RV-3.2)
            LOK=.TRUE.
          END IF
C..............................................................................
        ELSEIF(COPC.EQ.'3')THEN
C Galaxy: Cardelli, Clayton and Mathis (1989) -> Eqs. (2a), (2b), (4a), (4b),
C (5a), and (5b); with the update of O`Donnell (1994) for the Optical/NIR ->
C section 3
          YY=XX-1.82
          IF(XX.LT.0.3)THEN
          ELSEIF(XX.LT.1.1)THEN
            K_LAMBDA=0.574*(XX**1.61)*RV-0.527*(XX**1.61)
            LOK=.TRUE.
          ELSEIF(XX.LT.3.3)THEN
ccc         K_LAMBDA=(1.0+0.17699*YY-0.50447*YY*YY-0.02427*YY*YY*YY+
ccc     +       0.72085*YY*YY*YY*YY+0.01979*YY*YY*YY*YY*YY-
ccc     +       0.77530*YY*YY*YY*YY*YY*YY+0.32999*YY*YY*YY*YY*YY*YY*YY)*RV+
ccc     +       1.41338*YY+2.28305*YY*YY+1.07233*YY*YY*YY-
ccc     +       5.38434*YY*YY*YY*YY-0.62251*YY*YY*YY*YY*YY+
ccc     +       5.30260*YY*YY*YY*YY*YY*YY-2.09002*YY*YY*YY*YY*YY*YY*YY
            K_LAMBDA=(1.0+0.104*YY-0.609*YY*YY+0.701*YY*YY*YY+
     +       1.137*YY*YY*YY*YY-1.718*YY*YY*YY*YY*YY
     +       -0.827*YY*YY*YY*YY*YY*YY+1.647*YY*YY*YY*YY*YY*YY*YY
     +       -0.505*YY*YY*YY*YY*YY*YY*YY*YY)*RV+
     +       1.952*YY+2.908*YY*YY-3.989*YY*YY*YY-7.985*YY*YY*YY*YY
     +       +11.102*YY*YY*YY*YY*YY+5.491*YY*YY*YY*YY*YY*YY
     +       -10.805*YY*YY*YY*YY*YY*YY*YY+3.347*YY*YY*YY*YY*YY*YY*YY*YY
            LOK=.TRUE.
          ELSEIF(XX.LT.8.0)THEN
            IF(XX.LT.5.9)THEN
              FA=0.0
              FB=0.0
            ELSE
              FA=-0.04473*(XX-5.9)*(XX-5.9)-
     +         0.009779*(XX-5.9)*(XX-5.9)*(XX-5.9)
              FB=0.2130*(XX-5.9)*(XX-5.9)+
     +         0.1207*(XX-5.9)*(XX-5.9)*(XX-5.9)
            END IF
            K_LAMBDA=(1.752-0.316*XX-0.104/((XX-4.67)*(XX-4.67)+0.341)+
     +       FA)*RV-3.090+1.825*XX+1.206/((XX-4.62)*(XX-4.62)+0.263)+FB
            LOK=.TRUE.
          ELSEIF(XX.LT.10.0)THEN
            K_LAMBDA=(-1.073-0.628*(XX-8.0)+0.137*(XX-8.0)*(XX-8.0)-
     +       0.070*(XX-8.0)*(XX-8.0)*(XX-8.0))*RV+13.670+4.257*(XX-8.0)-
     +       0.420*(XX-8.0)*(XX-8.0)+0.374*(XX-8.0)*(XX-8.0)*(XX-8.0)
            LOK=.TRUE.
          END IF
C..............................................................................
        ELSEIF(COPC.EQ.'4')THEN
C Galaxy: Fitzpatrick (1999)
          IF(XX.LT.1./6.)THEN !unreliable beyond the M band at 6 microns
          ELSEIF(XX.LT.1./.2700)THEN !splines
            X(1)=0.00000
            Y(1)=0.00000*RV/3.1
            X(2)=1./2.6500
            Y(2)=0.26469*RV/3.1
            X(3)=1./1.2200
            Y(3)=0.82925*RV/3.1
            X(4)=1./.6000
            Y(4)=-0.422809+1.00270*RV+2.13572E-04*RV*RV
            X(5)=1./.5470
            Y(5)=-5.13540E-02+1.00216*RV-7.35778E-05*RV*RV
            X(6)=1./.4670
            Y(6)=7.00127E-01+1.00184*RV-3.32598e-05*RV*RV
            X(7)=1./.4110
            Y(7)=1.19456+1.01707*RV-5.46959E-03*RV*RV+
     +       7.97809E-04*RV*RV*RV-4.45636E-05*RV*RV*RV*RV
            X(8)=1./.2700
            Y(8)=FITZPATRICK99(X(8),RV)
            X(9)=1./.2600
            Y(9)=FITZPATRICK99(X(9),RV)
            NPTOS=9
            IMODE=1
            CALL CUBSPL__(X,Y,NPTOS,IMODE,SPLS,SPLA,SPLB,SPLC)
            I0=5
            CALL CUBSPLX__(X,Y,SPLA,SPLB,SPLC,NPTOS,I0,XX,K_LAMBDA)
            LOK=.TRUE.
          ELSEIF(XX.LT.9.0)THEN !fitting function
            K_LAMBDA=FITZPATRICK99(XX,RV)
            LOK=.TRUE.
          END IF
C..............................................................................
        ELSEIF(COPC.EQ.'a')THEN
C LMC: Howarth (1983)
          IF(XX.LT.1.83)THEN
            IF(RV.EQ.3.1)THEN
              K_LAMBDA=((1.86-0.48*XX)*XX-0.1)*XX
              LOK=.TRUE.
            END IF
          ELSEIF(XX.LT.2.75)THEN
            K_LAMBDA=RV+2.04*(XX-1.83)+0.094*(XX-1.83)*(XX-1.83)
            LOK=.TRUE.
          ELSEIF(XX.LT.9.0)THEN
            K_LAMBDA=RV-0.236+0.462*XX+0.105*XX*XX+
     +       0.454/((XX-4.557)*(XX-4.557)+0.293)
            LOK=.TRUE.
          END IF
C..............................................................................
        ELSEIF(COPC.EQ.'b')THEN
C LMC (30 Doradus): Fitzpatrick (1985)
          IF(XX.LT.3.3)THEN
          ELSEIF(XX.LE.5.9)THEN
            K_LAMBDA=0.493/((XX-4.60)*(XX-4.60)+0.501*0.501)
     +       -2.447+1.415*XX+RV
            LOK=.TRUE.
          ELSEIF(XX.LT.9.0)THEN
            K_LAMBDA=0.493/((XX-4.60)*(XX-4.60)+0.501*0.501)
     +       -2.447+1.415*XX+0.425*(XX-5.9)*(XX-5.9)
     +       -0.083*(XX-5.9)*(XX-5.9)*(XX-5.9)+RV
            LOK=.TRUE.
          END IF
C..............................................................................
        ELSEIF(COPC.EQ.'m')THEN
C SMC: Prevot et al (1984) + Bouchet et al. (1985)
          IF(XX.LT.1./2.2)THEN !limit K band
          ELSEIF(XX.LE.7.84)THEN
            X(30)=7.84
            Y(30)=13.54
            X(29)=7.52
            Y(29)=12.52
            X(28)=7.23
            Y(28)=11.51
            X(27)=6.98
            Y(27)=10.80
            X(26)=6.72
            Y(26)=9.84
            X(25)=6.48
            Y(25)=9.28
            X(24)=6.27
            Y(24)=9.06
            X(23)=6.07
            Y(23)=8.49
            X(22)=5.88
            Y(22)=8.01
            X(21)=5.70
            Y(21)=7.71
            X(20)=5.52
            Y(20)=7.17
            X(19)=5.38
            Y(19)=6.90
            X(18)=5.24
            Y(18)=6.76
            X(17)=5.00
            Y(17)=6.38
            X(16)=4.73
            Y(16)=5.85
            X(15)=4.50
            Y(15)=5.30
            X(14)=4.28
            Y(14)=4.53
            X(13)=4.09
            Y(13)=4.24
            X(12)=3.92
            Y(12)=3.91
            X(11)=3.75
            Y(11)=3.49
            X(10)=3.60
            Y(10)=3.15
            X(09)=3.46
            Y(09)=3.00
            X(08)=3.34
            Y(08)=2.65
            X(07)=3.22
            Y(07)=2.29
            X(06)=2.70
            Y(06)=1.67
            X(05)=2.35
            Y(05)=1.00
            X(04)=1.89
            Y(04)=0.00
            X(03)=1./1.25 !J band
            Y(03)=-2.02
            X(02)=1./1.65 !H band
            Y(02)=-2.36
            X(01)=1./2.20 !K band
            Y(01)=-2.47
            NPTOS=30
            K_LAMBDA=LININTERP__(NPTOS,X,Y,XX,IFLAG,N1,N2)+RV
            LOK=.TRUE.
          END IF
C..............................................................................
        ELSEIF(COPC.EQ.'p')THEN
C Starbursts: Calzetti (1997, proceeding)
          IF(XX.LT.1.0/1.0)THEN
          ELSEIF(XX.LE.1./0.63)THEN
            K_LAMBDA=((1.86-0.48*XX)*XX-0.1)*XX+1.73
            LOK=.TRUE.
          ELSEIF(XX.LE.1./0.12)THEN
            K_LAMBDA=2.656*(-2.156+1.509*XX-0.198*XX*XX+0.011*XX*XX*XX)
     +       +4.88
            LOK=.TRUE.
          END IF
C..............................................................................
        ELSEIF(COPC.EQ.'q')THEN
C Starbursts: Calzetti et al. (2000)
          IF(XX.LT.1.0/2.20)THEN
          ELSEIF(XX.LE.1./0.63)THEN
            K_LAMBDA=2.659*(-1.857+1.040*XX)+4.05
            LOK=.TRUE.
          ELSEIF(XX.LE.1./0.12)THEN
            K_LAMBDA=2.659*(-2.156+1.509*XX-0.198*XX*XX+0.011*XX*XX*XX)
     +       +4.05
            LOK=.TRUE.
          END IF
C..............................................................................
        END IF
C------------------------------------------------------------------------------
        END
C
C******************************************************************************
C
        REAL FUNCTION FITZPATRICK99(XX,RV)
        IMPLICIT NONE
        REAL XX,RV
C
        REAL C1,C2,C3,C4,DX,FX
C------------------------------------------------------------------------------
        C2=-0.824+4.717/RV
        C1=2.030-3.007*C2
        C3=3.23
        C4=0.41
        DX=XX*XX/
     +   ((XX*XX-4.596*4.596)*(XX*XX-4.596*4.596)+XX*XX*0.99*0.99)
        IF(XX.GT.5.9)THEN
          FX=0.5392*(XX-5.9)*(XX-5.9)+
     +     0.05644*(XX-5.9)*(XX-5.9)*(XX-5.9)
        ELSE
          FX=0.0
        END IF
        FITZPATRICK99=C1+C2*XX+C3*DX+C4*FX+RV
C
        END
C
C******************************************************************************
C REAL FUNCTION LININTERP__(N,X,Y,X0,IFLAG,N1,N2)
C
C Input: N,X,Y,X0
C Output: LININTERP__(function), IFLAG, N1, N2
C
C Performs a linear interpolation in the table X(N),Y(N) at x=X0. Note that the
C X matrix must be sorted in ascending order, although the 
C
C INTEGER N -> input number of data in X and Y
C REAL    X(N) -> data matrix 
C REAL    Y(N) -> data matrix 
C REAL    X0 -> x-point at which the linear interpolation is evaluated
C INTEGER IFLAG -> = 0 : interpolation
C                  = -1 : extrapolation towards lower X values
C                  = +1 : extrapolation towards higher X values
C                  = +9 : error (division by zero)
C INTEGER N1,N2 -> index of the X array within which X0 is found
C
C------------------------------------------------------------------------------
        REAL FUNCTION LININTERP__(N,X,Y,X0,IFLAG,N1,N2)
        IMPLICIT NONE
C
        INTEGER N
        REAL X(N),Y(N),X0
        INTEGER IFLAG
        INTEGER N1,N2
C local variables
        INTEGER I
C------------------------------------------------------------------------------
        IF(X0.LT.X(1))THEN
          IFLAG=-1
          N1=1
          N2=2
        ELSEIF(X0.GT.X(N))THEN
          IFLAG=+1
          N1=N-1
          N2=N
        ELSE
          IFLAG=0
          DO I=1,N
            IF(X(I).LE.X0) N1=I
          END DO
          DO I=N,1,-1
            IF(X(I).GE.X0) N2=I
          END DO
        END IF
C
        IF(N1.EQ.N2)THEN
          LININTERP__=Y(N1)
        ELSE
          IF(X(N1).NE.X(N2))THEN
            LININTERP__=Y(N1)+((X0-X(N1))/(X(N2)-X(N1)))*(Y(N2)-Y(N1))
          ELSE
            IFLAG=9
            LININTERP__=0.
          END IF
        END IF
C        
        END
C
C******************************************************************************
C SUBROUTINE CUBSPL__(X,Y,N,IMODE,S,A,B,C)
C
C Input: X,Y,N,IMODE,S
C Output: A,B,C,S
C
C This subroutine computes the coefficients of a cubic spline. See C.F. Gerald
C and P. O. Wheatley, in Applied Numerical Analysis, 4th edition, pag. 207.
C The subroutine returns the spline coefficients, where the spline defined
C in the interval between X(I),Y(I) and X(I+1),Y(I+1) is given by:
C
C      Y = A(I)*(X-X(I))**3 + B(I)*(X-X(I))**2 + C(I)*(X-X(I)) + D(I)
C
C REAL X(N) -> X-values to be fitted
C REAL Y(N) -> Y-values to be fitted
C INTEGER N -> number of data points
C INTEGER IMODE -> End conditions mode: if S(I) represent the second derivative
C                  at the point X(I),Y(I), the following four possibilites 
C                  are available:
C                  1) IMODE=1: S(0)=0, S(N)=0. This is called natural cubic
C                     spline. It is equivalent to assuming that the end cubics
C                     aproach linearity at their extremities.
C                  2) IMODE=2: S(0)=S(1), S(N)=S(N-1). This is equivalent to
C                     assuming that the cubics approach parabolas at their
C                     extremities.
C                  3) IMODE=3: S(0) is a linear extrapolation from S(1) and
C                     S(2), and S(N) is a linear extrapolation from S(N-2) 
C                     and S(N-1).
C                  4) IMODE=4: Force the slopes at each end to assume certain
C                     values.
C REAL S(N) -> if IMODE=4, in input S(1) and S(N) contain the first derivatives
C              at X(1) and X(N). In output, this matrix contains the second
C              derivatives
C REAL A(N) -> spline coefficients
C REAL B(N) -> spline coefficients
C REAL C(N) -> spline coefficients
C
Comment
C------------------------------------------------------------------------------
        SUBROUTINE CUBSPL__(X,Y,N,IMODE,S,A,B,C)
        IMPLICIT NONE
C
        INTEGER N
        REAL X(N),Y(N)
        INTEGER IMODE
        REAL S(N)
        REAL A(N),B(N),C(N)
C local parameters
        INTEGER NMAX
        PARAMETER (NMAX=100)
C local variables
        INTEGER I,I1,I2
        REAL DX1,DX2,DY1,DY2
        REAL SS(0:NMAX,4)
        REAL H
C------------------------------------------------------------------------------
C verificamos que al menos tenemos dos puntos
        IF(N.LT.2)THEN
          WRITE(*,100)'FATAL ERROR in subroutine CUBSPL__: '
          WRITE(*,100)' no. of data points too small: '
          WRITE(*,*)N
          STOP
        END IF
C verificamos que el numero de puntos no sea demasiado grande
        IF(N.GT.NMAX)THEN        !en realidad vale con N-1, pero lo dejamos asi
          WRITE(*,100)'FATAL ERROR in subroutine CUBSPL__: '
          WRITE(*,100)' no. of data points too large: '
          WRITE(*,*)N
          STOP
        END IF
C verificamos que IMODE toma un valor posible
        IF((IMODE.LT.1).OR.(IMODE.GT.4))THEN
          WRITE(*,100)'FATAL ERROR in subroutine CUBSPL__: '
          WRITE(*,100)' invalid IMODE value:'
          WRITE(*,*)IMODE
          STOP
        END IF
C------------------------------------------------------------------------------
C calculamos la matriz reducida, que contiene N-2 filas x 4 columnas
        DX1=X(2)-X(1)
        DY1=6.*(Y(2)-Y(1))/DX1
        DO I=1,N-2
          DX2=X(I+2)-X(I+1)
          DY2=6.*(Y(I+2)-Y(I+1))/DX2
          SS(I,1)=DX1
          SS(I,2)=2.*(DX1+DX2)
          SS(I,3)=DX2
          SS(I,4)=DY2-DY1
          DX1=DX2
          DY1=DY2
        END DO
C------------------------------------------------------------------------------
C modificamos la primera y ultima fila segun el valor de IMODE
        IF(IMODE.EQ.1)THEN              !no hay que modificar nada en este caso
c..............................................................................
        ELSEIF(IMODE.EQ.2)THEN
          SS(1,2)=SS(1,2)+(X(2)-X(1))
          SS(N-2,2)=SS(N-2,2)+(X(N)-X(N-1))
c..............................................................................
        ELSEIF(IMODE.EQ.3)THEN
          DX1=X(2)-X(1)
          DX2=X(3)-X(2)
          SS(1,2)=(DX1+DX2)*(DX1+2.*DX2)/DX2
          SS(1,3)=(DX2*DX2-DX1*DX1)/DX2
          DX1=X(N-1)-X(N-2)
          DX2=X(N)-X(N-1)
          SS(N-2,1)=(DX1*DX1-DX2*DX2)/DX1
          SS(N-2,2)=(DX2+DX1)*(DX2+2.*DX1)/DX1
c..............................................................................
        ELSEIF(IMODE.EQ.4)THEN
C notar que en este caso tambien cambia el taman~o de la matriz
          DX1=X(2)-X(1)
          DY1=(Y(2)-Y(1))/DX1
          SS(0,1)=1.
          SS(0,2)=2.*DX1
          SS(0,3)=DX1
          SS(0,4)=6.*(DY1-S(1))
          DX1=X(N)-X(N-1)
          DY1=(Y(N)-Y(N-1))/DX1
          SS(N-1,1)=DX1
          SS(N-1,2)=2.*DX1
          SS(N-1,3)=0.
          SS(N-1,4)=6.*(S(N)-DY1)
c..............................................................................
        END IF
C------------------------------------------------------------------------------
C dimensiones de la matriz a resolver
        IF(IMODE.EQ.4)THEN
          I1=1
          I2=N-1
        ELSE
          I1=2
          I2=N-2
        END IF
C resolvemos el sistema tridiagonal, eliminando en primer lugar los elementos
C que se encuentran por debajo de la diagonal
        DO I=I1,I2
          SS(I,1)=SS(I,1)/SS(I-1,2)
          SS(I,2)=SS(I,2)-SS(I,1)*SS(I-1,3)
          SS(I,4)=SS(I,4)-SS(I,1)*SS(I-1,4)
        END DO
C y ahora hacemos la sustitucion desde atras
        SS(I2,4)=SS(I2,4)/SS(I2,2)
        DO I=I2-1,I1-1,-1
          SS(I,4)=(SS(I,4)-SS(I,3)*SS(I+1,4))/SS(I,2)
        END DO
C------------------------------------------------------------------------------
C los valores de las derivadas segundas se almacenan en la matriz S
        DO I=I1-1,I2
          S(I+1)=SS(I,4)
        END DO
C
        IF(IMODE.EQ.1)THEN
          S(1)=0.
          S(N)=0.
        ELSEIF(IMODE.EQ.2)THEN
          S(1)=S(2)
          S(N)=S(N-1)
        ELSEIF(IMODE.EQ.3)THEN
          DX1=X(2)-X(1)
          DX2=X(3)-X(2)
          S(1)=((DX1+DX2)*S(2)-DX1*S(3))/DX2
          DX1=X(N-1)-X(N-2)
          DX2=X(N)-X(N-1)
          S(N)=((DX1+DX2)*S(N-1)-DX2*S(N-2))/DX1
        ELSEIF(IMODE.EQ.4)THEN                         !no hay nada que cambiar
        END IF
C------------------------------------------------------------------------------
C finalmente calculamos los coeficientes
        DO I=1,N-1
          H=X(I+1)-X(I)
          A(I)=(S(I+1)-S(I))/(6.*H)
          B(I)=S(I)/2.
          C(I)=(Y(I+1)-Y(I))/H-H*(2.*S(I)+S(I+1))/6.
        END DO
C------------------------------------------------------------------------------
100     FORMAT(A,$)
        END
C
C******************************************************************************
C SUBROUTINE CUBSPLX__(X,Y,A,B,C,N,I0,X0,Y0)
C
C Input: X,Y,A,B,C,N,I0,X0
C Output: Y0
C
C The subroutine returns the cubic spline evaluated at X0, using the
C coefficients determined in a previous call to CUBSPL__. The spline defined in 
C the interval between X(I),Y(I) and X(I+1),Y(I+1) is given by:
C
C      Y = A(I)*(X-X(I))**3 + B(I)*(X-X(I))**2 + C(I)*(X-X(I)) + D(I)
C
C If X0.LT.X(1), I=1 is employed (first computed spline)
C If X0.GT.X(N), I=N-1 is employed (last computed spline)
C
C REAL X(N) -> X-values fitted with CUBSPL__
C REAL Y(N) -> Y-values fitted with CUBSPL__
C REAL A(N) -> spline coefficients
C REAL B(N) -> spline coefficients
C REAL C(N) -> spline coefficients
C INTEGER N -> number of data points
C INTEGER I0 -> initial location to start the search of the place of X0 in
C               the X array
C REAL X0 -> X-value where the spline function will be evaluated
C REAL Y0 -> spline value at X0
C
Comment
C------------------------------------------------------------------------------
        SUBROUTINE CUBSPLX__(X,Y,A,B,C,N,I0,X0,Y0)
        IMPLICIT NONE
C
        INTEGER N
        REAL X(N),Y(N)
        REAL A(N),B(N),C(N)
        INTEGER I0
        REAL X0,Y0
C local variables
        REAL DX
C------------------------------------------------------------------------------
C buscamos el lugar en la tabla en la que se encuentra X0, para lo cual
C empleamos la subrutina BINSEARCH__, la cual permite emplear un valor de prueba
C para iniciar la busqueda, lo cual acelera el proceso cuando se realizan
C llamadadas sucesivas a esta funcion, con valores de X0 consecutivos.
        CALL BINSEARCH__(X,N,X0,I0)
        IF(I0.EQ.0) I0=1
        IF(I0.EQ.N) I0=N-1
C------------------------------------------------------------------------------
C evaluate the spline
        DX=X0-X(I0)
        Y0=Y(I0)+DX*(C(I0)+DX*(B(I0)+DX*A(I0)))
C------------------------------------------------------------------------------
        END
C
C******************************************************************************
C SUBROUTINE BINSEARCH__(X,N,X0,N0)
C
C Input: X,N,X0,N0
C Output: N0
C
C Given the array X(N), and the test value X0, this subroutine returns an
C integer N0, such that X0 is between X(N0) and X(N0+1). As input N0 is
C employed to start the searching. If X0.LT.X(1) then N0=0 on output, whereas 
C if X0.GT.X(N) then N0=N.
C
C REAL    X(N) -> ordered input array (not necesarilly equally-spaced)
C INTEGER N -> no. of points in input array
C REAL    X0 -> argument to be searched for
C INTEGER N0 -> location of X0 in the input array
C
Comment
C------------------------------------------------------------------------------
        SUBROUTINE BINSEARCH__(X,N,X0,N0)
        IMPLICIT NONE
C
        INTEGER N,N0
        REAL X(N),X0
C local variables
        INTEGER L,U,I
        INTEGER STEP
        LOGICAL LOOP
C------------------------------------------------------------------------------
        IF(N0.LT.1)THEN
ccc       WRITE(*,100)'* WARNING: in subroutine BINSEARCH__: '
ccc       WRITE(*,101)'N0.LT.1'
          N0=1
        END IF
        IF(N0.GT.N)THEN
ccc       WRITE(*,100)'* WARNING: in subroutine BINSEARCH__: '
ccc       WRITE(*,101)'N0.GT.N'
          N0=N
        END IF
C------------------------------------------------------------------------------
C Buscamos el intervalo inicial duplicando el paso de busqueda
        STEP=1
        L=N0
        LOOP=.TRUE.
c..............................................................................
        IF((X(1).LT.X(N)).EQV.(X0.GE.X(L)))THEN
          DO WHILE(LOOP)
            U=L+STEP
            IF(U.GT.N)THEN
              U=N+1
              LOOP=.FALSE.
            ELSE
              IF((X(1).LT.X(N)).EQV.(X0.GE.X(U)))THEN
                L=U
                STEP=2*STEP
              ELSE
                LOOP=.FALSE.
              END IF
            END IF
          END DO
c..............................................................................
        ELSE
          U=L
          DO WHILE(LOOP)
            L=U-STEP
            IF(L.LT.1)THEN
              L=0
              LOOP=.FALSE.
            ELSE
              IF((X(1).LT.X(N)).EQV.(X0.LT.X(L)))THEN
                U=L
                STEP=2*STEP
              ELSE
                LOOP=.FALSE.
              END IF
            END IF
          END DO
c..............................................................................
        END IF
C------------------------------------------------------------------------------
C Ahora buscamos el valor de N0 dividiendo el paso de busqueda
        DO WHILE(U-L.GT.1)
          I=(U+L)/2
          IF(X0.EQ.X(I))THEN
            N0=I
            RETURN
          END IF
          IF((X(1).LT.X(N)).EQV.(X0.GT.X(I)))THEN
            L=I
          ELSE
            U=I
          END IF
        END DO
C
        IF(U.LT.L)THEN
          WRITE(*,101) 'FATAL ERROR: in subroutine BINSEARCH__.'
          STOP
        END IF
C
        N0=L
C
ccc100     FORMAT(A,$)
101     FORMAT(A)
        END
