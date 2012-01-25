C------------------------------------------------------------------------------
C Version 23-May-1995                                           File:rebining.f
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
C SUBROUTINE REBINING(X,Y,N,XX,YY,M,XINI,XINC)
C
C Input: X,Y,N,XINI,XINC
C Output: XX,YY,M
C
C Rebin a data table X(1:N),Y(1:N). If X(J)-X(J-1) .GT. XINC the routine
C performs a linear interpolation.
C
C REAL    X(N) -> ordered input array (not necesarilly equally-spaced)
C REAL    Y(N) -> input array
C INTEGER N -> no. of points in input array
C REAL    XX(M) -> equally-spaced output array
C REAL    YY(M) -> output array
C INTEGER M -> no. of points in output array
C REAL    XINI -> equal to XX(1)
C REAL    XINC -> equal to XX(J)-XX(J-1) for all J
C
Comment
C------------------------------------------------------------------------------
        SUBROUTINE REBINING(X,Y,N,XX,YY,M,XINI,XINC)
        IMPLICIT NONE
        INTEGER N
        REAL X(N),Y(N)
        INTEGER M
        REAL XX(M),YY(M)
        REAL XINI,XINC
C
        INTEGER I,J,K,L
        INTEGER IOUT,JOUT
        REAL FTP                                        !declaracion de funcion
        REAL F1,F2,X0,FF
        LOGICAL LOUTRANGE
C------------------------------------------------------------------------------
        I=1                           !valor inicial de busqueda para BINSEARCH
        J=1                           !valor inicial de busqueda para BINSEARCH
        IOUT=0             !determinamos si estamos fuera de rango al principio
        JOUT=0                 !determinamos si estamos fuera de rango al final
        DO K=1,M
c---------
          !remember that BINSEARCH return the integer N0 such that X0 is
          !between X(N0) and X(N0+1)
          LOUTRANGE=.FALSE.
          XX(K)=XINI+REAL(K-1)*XINC
          CALL BINSEARCH(X,N,XX(K)-XINC/2.,I)        !buscamos limite izquierdo
          IF((I.EQ.0).OR.(I.EQ.N))THEN
            LOUTRANGE=.TRUE.
            IF(I.EQ.0) IOUT=K
          END IF
          CALL BINSEARCH(X,N,XX(K)+XINC/2.,J)          !buscamos limite derecho
          IF((J.EQ.0).OR.(J.EQ.N))THEN
            LOUTRANGE=.TRUE.
            IF((J.EQ.N).AND.(JOUT.EQ.0)) JOUT=K
          END IF
c.........
          IF(LOUTRANGE)THEN                             !si esta fuera de rango
ccc         WRITE(*,110)'WARNING: out of range in step',K
            YY(K)=0.
c.........
          ELSE                                       !si no esta fuera de rango
            IF(I.EQ.J)THEN                    !interpolacion normal tradicional
              YY(K)=FTP(X(I),X(I+1),Y(I),Y(I+1),XX(K))
            ELSE          !hay al menos un punto en [XX(K)-XINC/2,XX(K)+XINC/2]
              X0=XX(K)-XINC/2.
              F1=FTP(X(I),X(I+1),Y(I),Y(I+1),X0)      !funcion en extremo izdo.
              FF=(X(I+1)-X0)/XINC
              YY(K)=FTP(X0,X(I+1),F1,Y(I+1),(X0+X(I+1))/2.)*FF
              X0=XX(K)+XINC/2.
              F2=FTP(X(J),X(J+1),Y(J),Y(J+1),X0)      !funcion en extremo dcho.
              FF=(X0-X(J))/XINC
              YY(K)=YY(K)+FTP(X(J),X0,Y(J),F2,(X(J)+X0)/2.)*FF
              IF(J.GT.I+1)THEN           !contribucion de los puntos interiores
                DO L=I+1,J-1
                  FF=(X(L+1)-X(L))/XINC
                  YY(K)=YY(K)+0.5*(Y(L)+Y(L+1))*FF
                END DO
              END IF
            END IF
c.........
          END IF
c---------
        END DO
C
        IF(IOUT.NE.0)THEN
          WRITE(*,120)'WARNING: Out of range from...',1,'  to...',IOUT
        END IF
        IF(JOUT.NE.0)THEN
          WRITE(*,120)'WARNING: Out of range from...',JOUT,'  to...',M
        END IF
C
ccc110     FORMAT(A,I5)
120     FORMAT(A,I5,A,I5)
        END
C
C******************************************************************************
C Interpola linealmente entre (X1,Y1) y (X2,Y2) para una abcisa igual a X
        REAL FUNCTION FTP(X1,X2,Y1,Y2,X)
        IMPLICIT NONE
        REAL X1,X2,Y1,Y2,X
C
        REAL DX
C------------------------------------------------------------------------------
        IF(X2.EQ.X1)THEN
ccc       WRITE(*,'(A)')'FATAL ERROR: division by zero in function FTP.'
ccc       STOP
          FTP=(Y1+Y2)/2.
          RETURN
        END IF
C
        IF((X1.LE.X).AND.(X.LE.X2))THEN
          DX=(X-X1)/(X2-X1)
          FTP=(1.-DX)*Y1+DX*Y2
        ELSE
          WRITE(*,100) 'X1: '
          WRITE(*,*) X1
          WRITE(*,100) 'X2: '
          WRITE(*,*) X2
          WRITE(*,100) 'Y1: '
          WRITE(*,*) Y1
          WRITE(*,100) 'Y2: '
          WRITE(*,*) Y2
          WRITE(*,100) 'X : '
          WRITE(*,*) X
          WRITE(*,'(A)')'FATAL ERROR: out of range in function FTP.'
          STOP
        END IF
100     FORMAT(A,$)
        END
