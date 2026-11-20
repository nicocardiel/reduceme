C Version 03-April-2003
C******************************************************************************
C Calcula la extinción atmosférica para una longitud de onda X (en micras),
C con F1, F2 y F3 los factores correspondientes a la dispersión Rayleigh, los
C aerosoles y el Ozono, respectivamente. H es la altura del observatorio
C sobre el nivel del mar (en km). Para las fórmulas se ha utilizado:
C - Hayes & Latham, 1975, ApJ, 197, 593
C - King, 1985, RGO/La Palma technical note no 31
C - Hopp & Fernández, 2002, Calar Alto Newsletter no 4
C
C Para Calar Alto, y ajustando a los datos de la Figura 4 de Hopp & Fernández,
C hay que usar H=2.168 y
C F1=1.50, F2=1.25, F3=0.05
C F1=1.25, F2=1.00, F3=0.00 (en invierno)
C
        REAL FUNCTION FATMEXT(X,F1,F2,F3,H)
        IMPLICIT NONE
        REAL X
        REAL F1,F2,F3
        REAL H
C
        INTEGER NDAT_OZONO
        PARAMETER (NDAT_OZONO=19)
C
        REAL LININTERP_EXT
C
        INTEGER I
        INTEGER IFLAG,N1,N2
        REAL XL3(NDAT_OZONO),AL3(NDAT_OZONO)
        REAL A1,A2,A3
        REAL N
C------------------------------------------------------------------------------
C Curva de extinción por Ozono extraída de la digitalización de la Figura 1
C de Hopp y Fernández (2002)
        DATA(XL3(I),AL3(I),I=1,NDAT_OZONO)/
     +   0.400858492, 0.0000000000,
     +   0.445009947, 0.0034196426,
     +   0.467436016, 0.0078700306,
     +   0.490562797, 0.0140992925,
     +   0.512988687, 0.0221072268,
     +   0.534012854, 0.0318937153,
     +   0.557840288, 0.0416807123,
     +   0.579565465, 0.0470205247,
     +   0.597786546, 0.0496916547,
     +   0.618110359, 0.0479162075,
     +   0.635631025, 0.0425827615,
     +   0.651049197, 0.0381383374,
     +   0.665766656, 0.0319150127,
     +   0.680484116, 0.0248022731,
     +   0.695902407, 0.0185790136,
     +   0.719730496, 0.0123570720,
     +   0.741456151, 0.0061347620,
     +   0.765284181, 0.0025809477,
     +   0.808735013, 0.0000000000/
C------------------------------------------------------------------------------
        N=0.23465+1.076E2/(146.0-1./(X*X))+0.93161/(41.-1./(X*X))
        A1=9.4977E-3*N*N*EXP(-H/7.996)/(X*X*X*X)
        A1=A1*0.632817745 !ajuste a Figure 1 de Hopp & Fernandez (2002)
        A2=(X**(-0.8))*EXP(-H/1.5)
        A2=A2*0.198781848 !ajuste a Figure 1 de Hopp & Fernandez (2002)
        IF((X.LT.XL3(1)).OR.(X.GT.XL3(19)))THEN
          A3=0.0
        ELSE
          A3=LININTERP_EXT(NDAT_OZONO,XL3,AL3,X,IFLAG,N1,N2)/0.05
          IF(IFLAG.NE.0) STOP 'IFLAG.NE.0'
        END IF
C
        FATMEXT=F1*A1+F2*A2+F3*A3
        END
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
        REAL FUNCTION LININTERP_EXT(N,X,Y,X0,IFLAG,N1,N2)
        IMPLICIT NONE
C       
        INTEGER N
        REAL X(N),Y(N),X0
        INTEGER IFLAG
        INTEGER N1,N2
C------------------------------------------------------------------------------
        IF(X0.EQ.X(N))THEN !extremo superior (da division por cero abajo)
          IFLAG=0
          LININTERP_EXT=Y(N)
          RETURN
        ELSEIF(X0.LT.X(1))THEN !extrapolacion a la izquierda
          IFLAG=-1
          N1=1
          N2=2
        ELSEIF(X0.GT.X(N))THEN !extrapolacion a la derecha
          IFLAG=1
          N1=N-1
          N2=N
        ELSE !caso general
          IFLAG=0
          N1=N/2 !como inicio de busqueda tomamos el centro de la tabla
          CALL BINSEARCH_EXT(X,N,X0,N1)
          N2=N1+1
        END IF
C       
        IF(X(N1).NE.X(N2))THEN
          LININTERP_EXT=Y(N1)+((X0-X(N1))/(X(N2)-X(N1)))*(Y(N2)-Y(N1))
        ELSE
          IFLAG=9
          LININTERP_EXT=0.
        END IF
C       
        END
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
        SUBROUTINE BINSEARCH_EXT(X,N,X0,N0)
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
ccc       WRITE(*,100)'* WARNING: in subroutine BINSEARCH_EXT: '
ccc       WRITE(*,101)'N0.LT.1'
          N0=1
        END IF
        IF(N0.GT.N)THEN
ccc       WRITE(*,100)'* WARNING: in subroutine BINSEARCH_EXT: '
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
          WRITE(*,101)'FATAL ERROR: in subroutine BINSEARCH_EXT.'
          STOP
        END IF
C
        N0=L
C
ccc100     FORMAT(A,$)
101     FORMAT(A)
        END
