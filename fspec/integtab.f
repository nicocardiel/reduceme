C------------------------------------------------------------------------------
C Version 08-March-2005                                        File: integtab.f
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
C REAL FUNCTION INTEGTAB(N,X,Y,X1,X2,IFLAG1,IFLAG2)
C
C Input: N,X,Y,X1,X2
C Output: INTEGTAB(function), IFLAG1,IFLAG2
C
C Performs the integration of a function given in a tabular form, between the
C limits X1 and X2. Note that the X matrix must be sorted in ascending order.
C
C INTEGER N -> input number of data in X and Y
C REAL    X(N) -> data matrix 
C REAL    Y(N) -> data matrix 
C REAL    X1 -> first limit of the integral
C REAL    X2 -> second limit of the integral
C INTEGER IFLAG1 -> = 0 : interpolation
C                   = -1 : extrapolation towards lower X values
C                   = +1 : extrapolation towards higher X values
C                   = +9 : error (division by zero)
C INTEGER IFLAG2 -> = 0 : interpolation
C                   = -1 : extrapolation towards lower X values
C                   = +1 : extrapolation towards higher X values
C                   = +9 : error (division by zero)
C
Comment
C------------------------------------------------------------------------------
        REAL FUNCTION INTEGTAB(N,X,Y,X1,X2,IFLAG1,IFLAG2)
        IMPLICIT NONE
C       
        INTEGER N
        REAL X(N),Y(N),X1,X2
        INTEGER IFLAG1,IFLAG2
C funciones globales
        REAL LININTERP
C local variables
        INTEGER M1,M2,N1,N2
        INTEGER K,K1,K2
        REAL YSUM,Y1,Y2
        REAL XX1,XX2,YY1,YY2
C------------------------------------------------------------------------------
C protecciones
        IF(X1.GT.X2)THEN
          WRITE(*,100) 'X1,X2= '
          WRITE(*,*) X1,X2
          WRITE(*,101) 'FATAL ERROR: X1.GT.X2 in INTEGTAB'
          STOP
        END IF
        IF(X1.EQ.X2)THEN
          INTEGTAB=0.0
          RETURN
        END IF
C determinamos entre qué entradas de la tabla se ubican X1 y X2
        Y1=LININTERP(N,X,Y,X1,IFLAG1,M1,M2)
        Y2=LININTERP(N,X,Y,X2,IFLAG2,N1,N2)
        IF((IFLAG1.NE.0).OR.(IFLAG2.NE.0))THEN  !no permitimos nada raro, sorry
          INTEGTAB=0.0
          RETURN
        END IF
C sumamos el área de todos los rectángulos abarcados en el intervalo [X1,X2]
        YSUM=0.0
        K1=M1
        K2=N2-1
        IF(X2.EQ.X(N1)) K2=K2-1 !evita un efecto de borde por la forma en que
                                !trabaja BINSEARCH
        DO K=K1,K2 !bucle en número de rectángulos
          IF(K.EQ.K1)THEN !.......borde izquierdo del primer rectángulo a sumar
            XX1=X1
            YY1=Y1
          ELSE !..........borde izquierdo de un rectángulo que no es el primero
            XX1=X(K)
            YY1=Y(K)
          END IF
          IF(K.EQ.K2)THEN !.........borde derecho del último rectángulo a sumar
            XX2=X2
            YY2=Y2
          ELSE !.............borde derecho de un rectángulo que no es el último
            XX2=X(K+1)
            YY2=Y(K+1)
          END IF
          YSUM=YSUM+0.5*(YY1+YY2)*(XX2-XX1)/(X2-X1)
        END DO
        INTEGTAB=YSUM
C       
100     FORMAT(A,$)
101     FORMAT(A)
        END
