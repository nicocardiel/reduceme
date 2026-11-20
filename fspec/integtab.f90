!------------------------------------------------------------------------------
! Version 20-November-2026                                   File: integtab.f90
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
! REAL FUNCTION INTEGTAB(N,X,Y,X1,X2,IFLAG1,IFLAG2)
!
! Input: N,X,Y,X1,X2
! Output: INTEGTAB(function), IFLAG1,IFLAG2
!
! Performs the integration of a function given in a tabular form, between the
! limits X1 and X2. Note that the X matrix must be sorted in ascending order.
!
! INTEGER N -> input number of data in X and Y
! REAL    X(N) -> data matrix 
! REAL    Y(N) -> data matrix 
! REAL    X1 -> first limit of the integral
! REAL    X2 -> second limit of the integral
! INTEGER IFLAG1 -> = 0 : interpolation
!                   = -1 : extrapolation towards lower X values
!                   = +1 : extrapolation towards higher X values
!                   = +9 : error (division by zero)
! INTEGER IFLAG2 -> = 0 : interpolation
!                   = -1 : extrapolation towards lower X values
!                   = +1 : extrapolation towards higher X values
!                   = +9 : error (division by zero)
!
!Comment
!------------------------------------------------------------------------------
        REAL FUNCTION INTEGTAB(N,X,Y,X1,X2,IFLAG1,IFLAG2)
        IMPLICIT NONE
!       
        INTEGER N
        REAL X(N),Y(N),X1,X2
        INTEGER IFLAG1,IFLAG2
! funciones globales
! local variables
        INTEGER M1,M2,N1,N2
        INTEGER K,K1,K2
        REAL YSUM,Y1,Y2
        REAL XX1,XX2,YY1,YY2
!------------------------------------------------------------------------------
! protecciones
        IF(X1.GT.X2)THEN
          WRITE(*,100,ADVANCE='NO') 'X1,X2= '
          WRITE(*,*) X1,X2
          WRITE(*,101) 'FATAL ERROR: X1.GT.X2 in INTEGTAB'
          STOP
        END IF
        IF(X1.EQ.X2)THEN
          INTEGTAB=0.0
          RETURN
        END IF
! determinamos entre qué entradas de la tabla se ubican X1 y X2
        Y1=LININTERP(N,X,Y,X1,IFLAG1,M1,M2)
        Y2=LININTERP(N,X,Y,X2,IFLAG2,N1,N2)
        IF((IFLAG1.NE.0).OR.(IFLAG2.NE.0))THEN  !no permitimos nada raro, sorry
          INTEGTAB=0.0
          RETURN
        END IF
! sumamos el área de todos los rectángulos abarcados en el intervalo [X1,X2]
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
!       
100     FORMAT(A)
101     FORMAT(A)
        END
