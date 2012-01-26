C------------------------------------------------------------------------------
C Version 24-July-1998                                         File: downhill.f
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
C SUBROUTINE DOWNHILL_TOTAL(N,X0,DX0,YFUNK,A,B,G,YRMSTOL,XF,DXF,NEVAL)
C
C Input N,X0,DX0,YFUNK,A,B,G,YRMSTOL
C Output XF,NEVAL
C
C Minimization of the function YFUNK of N variables using the downhill 
C simplex method, as explained by Nelder and Mead (1965, Computer Journal, 7,
C pags. 308-313). The routine returns when the stopping criterion is reached 
C (the r.m.s. of the YFUNK values computed with all the vertices of the simplex
C is .LT. YRMSTOL), or when the number of function evaluations
C is too large (NEVAL.GT.NEVALMAX).
C
C INTEGER N -> number of variables
C REAL    X0(N) -> starting point (initial solution)
C REAL    DX0(N) -> N characteristic length scales, employed to derive N 
C                  aditional starting points which, together with X0, form 
C                  the (N+1) vertices of the simplex
C REAL    YFUNK -> function to be minimized
C REAL    A -> reflection coefficient (ALPHA); tipically, A=1.0
C REAL    B -> contraction coefficient (BETA); tipically, B=0.5
C REAL    G -> expansion coefficient (GAMMA); tipically, G=2.0
C REAL    YRMSTOL -> stopping criterion
C REAL    XF(N) -> final solution
C REAL    DXF(N) -> rms of XF evaluated from the final different points of the 
C                   simplex
C INTEGER NEVAL -> number of evaluations of YFUNK employed by DOWNHIL to reach 
C                  the solution; the routine returns NEVAL=-1 if something goes
C                  wrong
C
Comment
C------------------------------------------------------------------------------
	SUBROUTINE DOWNHILL_TOTAL(N,X0,DX0,YFUNK,A,B,G,YRMSTOL,
     +   XF,DXF,NEVAL)
	IMPLICIT NONE
C
	INTEGER NMAX
	PARAMETER (NMAX=20)                        !maximum number of variables
	INTEGER NEVALMAX
	PARAMETER (NEVALMAX=5000)  !maximum number of YFUNK evaluations allowed
C
	INTEGER N
	REAL X0(N),DX0(N)
	REAL YFUNK
	EXTERNAL YFUNK
	REAL A,B,G
	REAL YRMSTOL
	REAL XF(N)
	REAL DXF(N)
	INTEGER NEVAL
C local variables
	INTEGER I,J
	INTEGER I1                   !index of point with lowest function value
	INTEGER I2             !index of point with next-highest function value
	INTEGER I3                  !index of point with highest function value
	REAL P(NMAX,NMAX+1)                            !vertices of the simplex
	REAL PCEN(NMAX)               !(simplex without highest point) centroid
	REAL PSTAR1(NMAX),YSTAR1    !reflected point and function at this point
	REAL PSTAR2(NMAX),YSTAR2     !expanded point and function at this point
	REAL PSTAR3(NMAX),YSTAR3   !contracted point and function at this point
	REAL Y(NMAX+1)          !function values at the vertices of the simplex
	REAL X(NMAX)              !coordinates of a single point of the simplex
	REAL CYRMSTOL                                          !current YRMSTOL
	REAL FN,FN1                                   !float of 1/N and 1/(N+1)
	REAL YMEAN                  !mean value of the N+1 function evaluations
	REAL XMEAN      !mean single X value of the N+1 vertices of the simplex
C
	INTEGER QN,QNEVAL
	REAL QX0(NMAX),QDX0(NMAX)
	REAL QXF(NMAX),QDXF(NMAX)
	REAL QYRMSTOL
C------------------------------------------------------------------------------
C establecemos algunas protecciones y si algo falla retornamos NEVAL=-1
	NEVAL=-1
C
	DO J=1,N
	  XF(J)=X0(J)
	  DXF(J)=0.
	END DO
C
	IF(N.GT.NMAX)THEN
	  WRITE(*,101)'ERROR in subroutine DOWNHILL:'
	  WRITE(*,100)'N, NMAX: '
	  WRITE(*,*)N,NMAX
	  WRITE(*,101)'>>> N.GT.NMAX'
	  WRITE(*,100)'Press <CR> to continue...'
	  READ(*,*)
	  RETURN
	END IF
C
	DO J=1,N
	  IF(DX0(J).EQ.0.0)THEN
	    WRITE(*,101)'ERROR in subroutine DOWNHILL:'
	    DO I=1,N
	      WRITE(*,100)'Length scale (i,value): '
	      WRITE(*,*)I,DX0(I)
	    END DO
	    WRITE(*,101)'>>> characteristic length scale.EQ.0.0'
	    WRITE(*,100)'Press <CR> to continue...'
	    READ(*,*)
	    RETURN
	  END IF
	END DO
C------------------------------------------------------------------------------
	open(77,form='formatted')
	QN=3
	QX0(1)=A
	QX0(2)=B
	QX0(3)=G
	QDX0(1)=0.1
	QDX0(2)=0.1
	QDX0(3)=0.1
	WRITE(*,100)'Enter YRMSTOL for alfa,beta,gamma: '
	READ(*,*)QYRMSTOL
C calculamos los N+1 vertices del simplex
	DO I=1,QN+1                                          !numero de vertice
	  IF(I.EQ.QN+1)THEN                                !el valor inicial X0
	    DO J=1,QN
	      P(J,I)=QX0(J)
	    END DO
	  ELSE                !el valor inicial, modificando solo una dimension
	    DO J=1,QN
	      IF(J.EQ.I)THEN
	        P(J,I)=QX0(J)+QDX0(J)
	      ELSE
	        P(J,I)=QX0(J)
	      END IF
	    END DO
	  END IF
	END DO
C evaluamos la funcion en cada uno de los vertices
	DO I=1,QN+1
	  DO J=1,QN
	    X(J)=P(J,I)
	  END DO
	  CALL DOWNHILL(N,X0,DX0,YFUNK,X(1),X(2),X(3),YRMSTOL,
     +     XF,DXF,NEVAL)
	  Y(I)=REAL(NEVAL)
	  write(77,*)0,(p(j,i),j=1,3),y(i),xf(1)
	END DO
C numero de evaluaciones de la funcion (no contamos las realizadas para
C inicializar el simplex)
	QNEVAL=0
C real de 1/N y de 1/(N+1)
	FN=1./REAL(QN)
	FN1=1./REAL(QN+1)
C------------------------------------------------------------------------------
C calculamos el numero de vertice con el valor minimo (I1), el maximo (I3) y 
C el siguiente al maximo (I2)
10	IF(Y(1).LE.Y(2))THEN
	  I3=2
	  I2=1
	ELSE
	  I3=1
	  I2=2
	END IF
	I1=I2
	IF(QN.GT.1)THEN
	  DO I=3,QN+1
	    IF(Y(I).LT.Y(I1))THEN
	      I1=I
	    ELSE
	      IF(Y(I).GT.Y(I3))THEN
	        I2=I3
	        I3=I
	      ELSE
	        IF(Y(I).GT.Y(I2)) I2=I
	      END IF
	    END IF
	  END DO
	END IF
C------------------------------------------------------------------------------
	write(77,*)qneval,(p(j,i1),j=1,3),y(i3),xf(1)
	write(77,*)qneval,(p(j,i2),j=1,3),y(i3),xf(1)
	write(77,*)qneval,(p(j,i3),j=1,3),y(i3),xf(1)
C comprobamos si se satisface la condicion de salida
	YMEAN=0.
	DO I=1,QN+1
	  YMEAN=YMEAN+Y(I)
	END DO
	YMEAN=YMEAN*FN1
	CYRMSTOL=0.
	DO I=1,QN+1
	  CYRMSTOL=CYRMSTOL+(Y(I)-YMEAN)*(Y(I)-YMEAN)
	END DO
	CYRMSTOL=SQRT(CYRMSTOL*FN)
	IF((CYRMSTOL.LT.QYRMSTOL).OR.(QNEVAL.GE.NEVALMAX))THEN            !EXIT
	  DO J=1,N  !retornamos como solucion el vertice con el valor minimo I1
	    QXF(J)=P(J,I1)
	    XMEAN=0.
	    DO I=1,N+1
	      XMEAN=XMEAN+P(J,I)
	    END DO
	    XMEAN=XMEAN*FN1
	    QDXF(J)=0.
	    DO I=1,N+1
	      QDXF(J)=QDXF(J)+(P(J,I)-XMEAN)*(P(J,I)-XMEAN)
	    END DO
	    QDXF(J)=SQRT(QDXF(J)*FN)
	  END DO
	print*,'* valores optimos, NEVAL:',y(i1)
	print*,'alfa,  (rms): ',qxf(1),qdxf(1)
	print*,'beta,  (rms): ',qxf(2),qdxf(2)
	print*,'gamma, (rms): ',qxf(3),qdxf(3)
	CALL DOWNHILL(N,X0,DX0,YFUNK,A,B,G,YRMSTOL,
     +   XF,DXF,NEVAL)
	print*,'* valores iniciales, NEVAL:',NEVAL
	print*,'alfa,  (rms): ',a
	print*,'beta,  (rms): ',b
	print*,'gamma, (rms): ',g
	close(77)
	  RETURN
	END IF
C------------------------------------------------------------------------------
C calculamos el centroide de los vertices del simplex sin el punto I3
	DO J=1,QN
	  PCEN(J)=0.
	END DO
C sumamos vertices
	DO I=1,QN+1
	  IF(I.NE.I3)THEN
	    DO J=1,QN
	      PCEN(J)=PCEN(J)+P(J,I)
	    END DO
	  END IF
	END DO
C normalizamos
	DO J=1,QN
	  PCEN(J)=PCEN(J)*FN
	END DO
C------------------------------------------------------------------------------
C ponemos en marcha la "maquina"
	DO J=1,QN
	  PSTAR1(J)=(1.+A)*PCEN(J)-A*P(J,I3)                        !reflection
	END DO
	write(77,'(a)')'alfa'
	CALL DOWNHILL(N,X0,DX0,YFUNK,PSTAR1(1),PSTAR1(2),PSTAR1(3),
     +   YRMSTOL,XF,DXF,NEVAL)
	YSTAR1=REAL(NEVAL)
	QNEVAL=QNEVAL+1
	IF(YSTAR1.LE.Y(I1))THEN              !reflection produced a new minimum
	  DO J=1,QN
	    PSTAR2(J)=G*PSTAR1(J)+(1.-G)*PCEN(J)   !try an additional expansion
	  END DO
	  write(77,'(a)')'gamma'
	  CALL DOWNHILL(N,X0,DX0,YFUNK,PSTAR2(1),PSTAR2(2),PSTAR2(3),
     +     YRMSTOL,XF,DXF,NEVAL)
	  YSTAR2=REAL(NEVAL)
	  QNEVAL=QNEVAL+1
	  IF(YSTAR2.LE.Y(I1))THEN               !additional expansion succeeded
	    DO J=1,QN                          !replace highest point by PSTAR2
	      P(J,I3)=PSTAR2(J)
	    END DO
	    Y(I3)=YSTAR2
	  ELSE                            !additional expansion did not succeed
	    DO J=1,QN                          !replace highest point by PSTAR1
	      P(J,I3)=PSTAR1(J)
	    END DO
	    Y(I3)=YSTAR1
	  END IF
	ELSE                          !reflection did not produce a new minimum
	  IF(YSTAR1.GE.Y(I2))THEN
	    IF(YSTAR1.LT.Y(I3))THEN     !reflected point is better than highest
	      DO J=1,QN                        !replace highest point by PSTAR1
	        P(J,I3)=PSTAR1(J)
	      END DO
	      Y(I3)=YSTAR1
	    END IF
	    DO J=1,QN
	      PSTAR3(J)=B*P(J,I3)+(1.-B)*PCEN(J)             !try a contraction
	    END DO
	    write(77,'(a)')'beta'
	    CALL DOWNHILL(N,X0,DX0,YFUNK,PSTAR3(1),PSTAR3(2),PSTAR3(3),
     +       YRMSTOL,XF,DXF,NEVAL)
	    YSTAR3=REAL(NEVAL)
	    QNEVAL=QNEVAL+1
	    IF(YSTAR3.GE.Y(I3))THEN                !contraction did not succeed
	      DO I=1,QN+1        !multiple contraction towards the lowest point
	        IF(I.NE.I1)THEN
	          DO J=1,QN
	            X(J)=0.5*(P(J,I)+P(J,I1))
	            P(J,I)=X(J)
	          END DO
	          CALL DOWNHILL(N,X0,DX0,YFUNK,X(1),X(2),X(3),
     +             YRMSTOL,XF,DXF,NEVAL)
	          Y(I)=REAL(NEVAL)
	          QNEVAL=QNEVAL+1
	        END IF
	      END DO
	      write(77,'(a)')'multiple contraction'
	    ELSE                                     !the contraction succeeded
	      DO J=1,QN                        !replace highest point by PSTAR3
	        P(J,I3)=PSTAR3(J)
	      END DO
	      Y(I3)=YSTAR3
	    END IF
	  ELSE    !initial reflected point lies between lowest and next-highest
	    DO J=1,QN                          !replace highest point by PSTAR1
	      P(J,I3)=PSTAR1(J)
	    END DO
	    Y(I3)=YSTAR1
	  END IF
	END IF
	GOTO 10
C------------------------------------------------------------------------------
100	FORMAT(A,$)
101	FORMAT(A)
	END
